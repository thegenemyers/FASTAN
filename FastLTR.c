/*******************************************************************************************
 *
 *  Search an assembly for satellitic repeats
 *
 *  Author:   Gene Myers
 *  Creation: Jan 2024
 *  Last Mod: July 2025
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <dirent.h>

#include "GDB.h"
#include "ANO.h"
#include "align.h"
#include "alncode.h"

#undef   PROLOG
#undef   SORT1
#undef   SORT2
#undef   SHOW_SEEDS
#undef   SHOW_SEARCH
#undef   SHOW_ALIGNMENTS

#define DIAG_MAX       20000
#define DIAG_MIN        3000
#define BLOCK_OVERLAP  20000  //  Must be less than 32000

#define TSPACE     100
#define VERSION  "0.5"
#define MBUF_LEN   200  //  Must be even

static char *Usage = "[-vM] [-T(8)] <source:path>[<fa_extn>|<1_extn>] <target>[.1aln]";

static int NTHREADS;
static int VERBOSE;
static int MAKE_MASK;

static char dna[4] = { 'a', 'c', 'g', 't' };

static void Print_Seq(uint8 *seq, int len)
{ int j;

  for (j = 0; j < len; j++)
    printf("%c",dna[seq[j]]);
}

static char *emer(int x, int unit)
{ static char mer[9];
  int i;

  mer[unit] = '\0';
  for (i = unit-1; i >= 0; i--)
    { mer[i] = dna[x&0x3];
      x >>= 2;
    }
  return (mer);
}

/*******************************************************************************************
 *
 *  DIAGONAL HITS DETECTOR
 *
 ********************************************************************************************/

  //  Thread bundle

typedef struct
   { int         tid;
     OneFile    *ofile;
     GDB        _gdb, *gdb;

     Work_Data  *work;     //  alignment machinery
     Align_Spec *spec;
     Overlap     _over, *over;
     Alignment   _align, *align;

     uint8      *buffer;   //  contig buffer

     void       *block;    //  memory block for diagonal analyzer (524MB)

     int         tmax;     //  i64 trace vector for 1-file
     int64      *trace;

     OneFile    *mfile;
     int64      *mlist;
     char       *mstring;

     int         mint;
     char       *mstr;
     int         mscf;
     int64       moff;
   } S_Bundle;

  //  Return average diagonal of trace points + the 2 end points

static int ave_tp_diag(Path *path)
{ int     tlen;
  uint16 *trace;
  int64   ave;
  int     ab, bb;
  int     i;

  tlen  = path->tlen-2;
  trace = path->trace;

  ab = path->abpos;
  bb = path->bbpos;
  ave = ab-bb;
  ab = (ab/TSPACE)*TSPACE;
  for (i = 1; i < tlen; i += 2)
    { ab = ab + TSPACE;
      bb = bb + trace[i];
      ave += (ab-bb); 
    }
  ave += path->aepos - path->bepos;
  return ((int) ((ave/(tlen/2+2.))+.5));
}

typedef struct
  { uint16  diag;
    uint16  ibeg;
  } Seed;

typedef struct
  { uint16 diag;
    uint16 count;
  } Chord;

static int CSORT(const void *l, const void *r)
{ Chord *x = (Chord *) l;
  Chord *y = (Chord *) r;

  return (y->count - x->count);
}

static int spectrum_block(uint8 *seq, int off, int len, S_Bundle *bundle)
{ Alignment  *align = bundle->align;
  Overlap    *over  = bundle->over;
  Align_Spec *spec  = bundle->spec;
  Work_Data  *work  = bundle->work;
  OneFile    *ofile = bundle->ofile;
  int64      *t64   = bundle->trace;
  int         tmax  = bundle->tmax;

  int64      *mlist   = bundle->mlist;
  char       *mstring = bundle->mstring;
  OneFile    *mfile   = bundle->mfile;
  int         mint    = bundle->mint;
  char       *mstr    = bundle->mstr;
  int         mscf    = bundle->mscf;
  int64       moff    = bundle->moff;

  int     i, p, x, c;
  int     d, e, f;
  uint16  kmer;
  uint16 *index;  // 0x08000
  uint16 *count;  // 0x10000
  uint16 *diags;  // DIAG_MAX < 0x08000
  Seed   *post;   // 0x08000
  Seed   *hits;   // 0x08000
  uint8  *s7;
  int     l7;
  double  freq[4];

  count = (uint16 *) bundle->block;
  index = count + 0x10000;
  diags = index + 0x08000; 
  post  = (Seed *) (diags +  0x08000);
  hits  = post + 0x08000;

  (void) off;

#ifdef PROLOG
  printf("\nPANEL %d-%d\n",off,off+0x8000);
  fflush(stdout);
#endif

  for (i = 0; i < 0x10000; i++)   //  Init counters
    count[i] = 0;

  s7 = seq+7;
  l7 = len-7;

  { int    fqi[4];
    uint16 u;

    for (i = 0; i < 4; i++)
      fqi[i] = 0;

    kmer = seq[0];                 //  count # of each 8-mer
    fqi[kmer] = 1;
    for (i = 1; i < 7; i++)
      { u = seq[i];
        kmer = (kmer << 2) | u;
        fqi[u] += 1;
      }
    for (i = 0; i < l7; i++)
      { u = s7[i];
        kmer = (kmer << 2) | u;
        count[kmer] += 1;
        fqi[u] += 1;
      }

    for (i = 0; i < 4; i++)
      freq[i] = (1.*fqi[i]) / len;
  }

  p = 0;                         //  turn counts into ptrs
  for (i = 0; i < 0x10000; i++)
    { x = count[i];
      count[i] = p;
      p += x;
    }

  kmer = seq[0];                 //  place positions in index in order of 8-mer
  for (i = 1; i < 7; i++)
    kmer = (kmer << 2) | seq[i];
  for (i = 0; i < l7; i++)
    { kmer = (kmer << 2) | s7[i];
      index[count[kmer]++] = i;
    }
  index[l7] = 0;

  index[0] |= 0x8000;           //  mark bucket ends and reset count
  for (i = 0; i < 0xffff; i++)
    { index[count[i]] |= 0x8000;
      count[i] = 0;
    }
  count[0xffff] = 0;

#ifdef SORT1
  for (i = 0; i < l7; i++)
    { f = index[i];
      p = index[i] & 0x7fff;
      printf("%c %5d: ",p==f?' ':'+',p);
      Print_Seq(seq+p,8);
      printf("\n");
    }
#endif

  e = index[0] & 0x7fff;
  for (i = 1; i < l7; i++)     //  count ibeg's of all same-kmer adjacent position pairs
    { f = index[i];            //    that are within diag < DIAG_MAX (8Kbp)
      if (f < 0x8000)
        { d = f-e;
          if (d < DIAG_MAX)
            count[e] += 1;
          e = f;
        }
      else
        e = f & 0x7fff;
    }

  p = 0;                          //  turn counts into ptrs
  for (i = 0; i < 0x08000; i++)
    { x = count[i];
      count[i] = p;
      p += x;
    }

  for (i = 0; i < DIAG_MAX; i++)  //  init diagonal tube counters
    diags[i] = 0;

  e = index[0] & 0x7fff;
  for (i = 1; i < l7; i++)       //   place seed pairs in post sorted on ibeg using count
    { f = index[i];              //     ptrs.  Also count diagonal tubes for next sort.
      if (f < 0x8000)
        { d = f-e;
          if (d < DIAG_MAX)
            { c = count[e]++;
              post[c].ibeg = e;
              post[c].diag = d;
              diags[d] += 1;
            }
          e = f;
        }
      else
        e = f & 0x7fff;
    }

  p = 0;                        //  turn diag counts into sort ptrs
  for (i = 0; i < DIAG_MAX; i++)
    { x = diags[i];
      diags[i] = p;
      p += x;
    }

#ifdef SORT2
  printf("Sorted on Anti\n");
  for (c = 0; c < p; c++)
    printf(" %5d %5d\n",post[c].diag,post[c].ibeg); 
#endif

  for (i = 0; i < p; i++)       //  place ibeg/diag pairs in hits in order of diag then ibeg
    { c = post[i].diag;
      hits[diags[c]++] = post[i];
    }

#ifdef SHOW_SEEDS
  p = 0;
  for (i = 1; i < DIAG_MAX; i++)
    { f = diags[i];
      if (p >= e)
        continue;
      printf("Diagonal %d : %d\n",i,f);
      for ( ; p < f; p++)
        { d = hits[p].diag;
          e = hits[p].ibeg;
          printf("   %4d : %5d  ",d,e);
          Print_Seq(seq+e,8);
          printf("\n");
        }
    }
#endif

  { int   ncnt;
    int   outhit, end, beg;
    Chord *hist = (Chord *) post;
    int   unit, wide, anti, last;
    Path *path;

    ncnt = 0;
    p = diags[1];
    for (i = 2; i < DIAG_MAX; i++)
      { f = diags[i];
        if (f-p > 1 && f-p > (i>>6))
          { hist[ncnt].count = f-p;
            hist[ncnt].diag  = i;
            ncnt += 1;
          }
        p = f;
      }

    for (i = 0; i < 0x10000; i++)   //  Init counters for model subroutine
      count[i] = 0;

    qsort(hist,ncnt,sizeof(Chord),CSORT);

#ifdef SHOW_SEARCH
    printf("Histo: %d\n",ncnt);
    for (i = 0; i < ncnt && i < 100; i++)
      printf(" %4d: %5d\n",hist[i].diag,hist[i].count);
#endif

    outhit = 0;
    for (i = 0; i < ncnt; i++)
      { d = hist[i].diag;
        if (d < DIAG_MIN || d > DIAG_MAX)
          continue;

        last = -1;
#ifdef SHOW_SEARCH
        printf(" %4d: %5d\n",d,diags[d]-diags[d-1]);
#endif
        for (x = diags[d-1]+1; x < diags[d]; x++)
          { p = hits[x].ibeg;
#ifdef SHOW_SEARCH
            printf("  p = %d (%d %d) %d\n",p,last,hits[x-1].ibeg,off);
#endif
            if (p < last || p - hits[x-1].ibeg > d || seq[p] >= 4 || seq[p+7] >= 4)
              continue;

            wide = .2*d;
            if (wide < 1)
              wide = 1;
            anti = 2*(off + p) + d;
            Local_Alignment(align,work,spec,d,d,anti,wide,wide);

            path = align->path;
#ifdef SHOW_SEARCH
            printf("    %d (%d)  %d-%d-%d\n",path->aepos-path->bbpos,2*d,d-wide,d,d+wide);
            printf(" Hit spans %d-%d (unit = %d)\n",path->bbpos,path->aepos,unit);
#endif

            end = path->aepos - off;
            beg = path->bbpos - off;

            if (beg > p || end <= p)
              continue;
            if (end > last)
              last = end;
            if (path->bepos >= path->abpos)
              continue;

            unit = ave_tp_diag(path);

            if (over->path.tlen > tmax)
              { tmax = bundle->tmax = 1.2*over->path.tlen + 1000;
                t64  = bundle->trace = realloc(t64,sizeof(int64)*tmax);
              }
            Write_Aln_Overlap(ofile,over);
            Compress_TraceTo8(over,0);
            Write_Aln_Trace(ofile,over->path.trace,over->path.tlen,t64,unit);

            if (MAKE_MASK)
              { if (mint >= MBUF_LEN)
                  { oneInt(mfile,0) = mscf;
                    oneWriteLine(mfile,'M',mint,mlist);
                    oneWriteLine(mfile,'L',mstr-mstring,mstring);
printf("list =");
for (i = 0; i < mint; i++)
  printf(" %lld",mlist[i]);
printf("\n");
*mstr = '\0';
printf("string = '%s'\n",mstring);
                    mint = 0;
                    mstr = mstring;
                  }
                if (mint == 0)
                  mstr += sprintf(mstr,"%d",unit);
                else
                  mstr += sprintf(mstr,"\t%d",unit);
                mlist[mint++] = path->bbpos + moff;
                mlist[mint++] = path->aepos + moff;
              }

#if defined(SHOW_ALIGNMENTS) || defined(DO_CUT)
            Decompress_TraceTo16(over);
            Compute_Trace_PTS(align,work,100,GREEDIEST,d-wide,d+wide);
#endif

#ifdef SHOW_ALIGNMENTS
#ifndef SHOW_SEARCH
            if (last < 0)
              printf(" %4d: %5d\n",d,diags[d]-diags[d-1]);
            printf("\n");
#endif
            printf(" Hit spans %d-%d (unit = %d)\n",path->bbpos,path->aepos,unit);
            Print_Reference(stdout,align,work,8,100,10,0,10,0);
#endif

            if (path->bepos < path->abpos)
              { for (f = (path->abpos-off); f < end; f++)
                  seq[f] = 4;
                end = path->bepos-off;
                for (f = beg; f < end; f++)
                  seq[f] = 4;

#ifdef SHOW_ALIGNMENTS
                printf("  Near Tandem len = %d gap = %d\n",
                       path->aepos-path->bbpos,path->abpos-path->bepos);
#endif
              }
            else
              { for (f = beg; f < end; f++)
                  seq[f] = 4;
              }

            if (path->aepos > outhit)
              outhit = path->aepos;
          }
      }

    bundle->mint = mint;
    bundle->mstr = mstr;

    return (outhit);
  }
}


/*******************************************************************************************
 *
 *  THREADS: ONE PER CONTIG
 *
 ********************************************************************************************/

static pthread_mutex_t TMUTEX;
static pthread_cond_t  TCOND;

//  Tstack[0..Tavail-1] is a stack of available threads at any moment.
//  It is always manipulated inside the mutex TMUTEX

static int *Tstack;
static int  Tavail;

//  Thread to process a contig of a GDB.  Set up thread's personal data structures and
//    then process the contig in 32Kbp blocks overlaping by 8Kbp.  The one exception to
//    this is if a very long satellite whose alignment extends beyond the current 32Kbp
//    block, in which case the next block begins at the end of this alignment.

static void *compress_thread(void *args)
{ S_Bundle *bundle = (S_Bundle *) args;
  uint8    *buffer = bundle->buffer;
  GDB      *gdb    = bundle->gdb;
  int       last, clen;
  int       i, p;
 
  i = bundle->over->aread;

  Get_Contig(gdb,i,NUMERIC,(char *) buffer);
#ifdef PROLOG
  printf("CONTIG %d\n",i+1);
#endif

  clen = gdb->contigs[i].clen;
  bundle->align->aseq  = bundle->align->bseq = (char *) buffer;
  bundle->align->alen  = bundle->align->blen = clen;
  bundle->over->bread  = i;

  bundle->mint = 0;
  bundle->mstr = bundle->mstring;
  bundle->mscf = gdb->contigs[i].scaf;
  bundle->moff = gdb->contigs[i].sbeg;

  last = -1;
  if (clen < 0x8000)
    spectrum_block(buffer,0,clen,bundle);
  else
    for (p = 0; p+0x2000 <= clen; p += BLOCK_OVERLAP)
      { if (p+0x8000 > clen)
          spectrum_block(buffer+p,p,clen-p,bundle);
        else
          { last = spectrum_block(buffer+p,p,0x8000,bundle);
            if (last >= p+BLOCK_OVERLAP)
              p = last-BLOCK_OVERLAP;
          }
      }

  if (bundle->mint > 0)
    { oneInt(bundle->mfile,0) = bundle->mscf;
      oneWriteLine(bundle->mfile,'M',bundle->mint,bundle->mlist);
      oneWriteLine(bundle->mfile,'L',bundle->mstr-bundle->mstring,bundle->mstring);
    }
 
  pthread_mutex_lock(&TMUTEX);   //  Put this thread back on the avail stack
    Tstack[Tavail++] = bundle->tid;
  pthread_mutex_unlock(&TMUTEX);

  pthread_cond_signal(&TCOND);   //  Signal a thread is available

  return (NULL);
}


/*******************************************************************************************
 *
 *  MAIN
 *
 ********************************************************************************************/

int main(int argc, char *argv[])
{ FILE     **units;
  char      *spath;
  GDB       _gdb, *gdb = &_gdb;
  OneFile   *Ofile;
  OneFile   *Mfile;
  OneSchema *anoSchema;

  (void) Print_Seq;
  (void) emer;

  //   Process command line

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("FasTAN")

    NTHREADS = 8;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vM")
            break;
          case 'T':
            ARG_NON_NEGATIVE(NTHREADS,"number of threads to use");
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE   = flags['v'];
    MAKE_MASK = flags['M'];

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        fprintf(stderr,"      -M: Make a .1ano mask of the hits found.\n");
        exit (1);
      }
  }

  //  Get GDB or make a temporary if a fasta

  { char *cpath, *APATH, *AROOT;
    int   ftype;

    ftype = Get_GDB_Paths(argv[1],NULL,&spath,&cpath,0);

    free(cpath);

    if (MAKE_MASK && ftype != IS_GDB)
      { fprintf(stderr,"%s: The source must be a GDB when masking (-m) is on\n",Prog_Name);
        exit (1);
      }
  
    units = Get_GDB(gdb,spath,".",NTHREADS,NULL);

    //  Open 1aln file for threaded writing

    APATH = PathTo(argv[2]);
    AROOT = Root(argv[2],".1aln");
    cpath = getcwd(NULL,0);

    Ofile = open_Aln_Write(Catenate(APATH,"/",AROOT,".1aln"),NTHREADS,Prog_Name,VERSION,
                           Command_Line,TSPACE,spath,NULL,cpath);

    if (MAKE_MASK)
      { anoSchema = make_ANO_Schema();
        Mfile = oneFileOpenWriteNew(Catenate(APATH,"/",AROOT,".1ano"),anoSchema,"ano",1,NTHREADS);

        oneAddProvenance(Mfile,Prog_Name,VERSION,Command_Line);

        oneAddReference(Mfile,gdb->srcpath,1);

        Write_Skeleton(Mfile,gdb);
      }
    free(cpath);
    free(AROOT);
    free(APATH);

    Write_Skeleton(Ofile,gdb);
  }

  if (VERBOSE)
    { fprintf(stderr,"\n  Database loaded, begin scan of %d contigs\n\n",gdb->ncontig);
      fflush(stderr);
    }

  StartTime();

  { int       i, tid;
    int       done, launch;
    pthread_t threads[NTHREADS];
    S_Bundle  parm[NTHREADS];
    int       tstack[NTHREADS];

    for (i = 0; i < NTHREADS; i++)
      { parm[i].tid   = i;
        parm[i].ofile = Ofile + i;
        parm[i].gdb   = gdb;
        parm[i]._gdb  = _gdb;
        parm[i]._gdb.seqs = units[i];
        parm[i].work  = New_Work_Data();
        if (i == 0)
          parm[i].spec = New_Align_Spec(.7,TSPACE,gdb->freq,0);
        else
          parm[i].spec = parm[i-1].spec;
        parm[i].over  = &parm[i]._over;
        parm[i].align = &parm[i]._align;
        parm[i].align->path  = &(parm[i]._over.path);
        parm[i].align->flags = 0;
        parm[i].over->flags  = 0;
        parm[i].block   = malloc(9*0x10000);   // 576KB
        parm[i].buffer  = ((uint8 *) malloc(gdb->maxctg + 4)) + 1;
        if (MAKE_MASK)
          { parm[i].mfile   = Mfile + i;
            parm[i].mlist   = malloc(sizeof(int64)*MBUF_LEN);
            parm[i].mstring = malloc(11*MBUF_LEN);
          }
        parm[i].tmax    = 10000;
        parm[i].trace   = malloc(sizeof(int64)*10000);

        if (parm[i].block == NULL || parm[i].buffer == NULL || parm[i].trace == NULL)
          { fprintf(stderr,"%s: Not enough memory\n",Prog_Name);
            exit (1);
          }
        if (MAKE_MASK && (parm[i].mlist == NULL || parm[i].mstring == NULL))
          { fprintf(stderr,"%s: Not enough memory\n",Prog_Name);
            exit (1);
          }
      }

    Tstack = tstack;
    for (i = 0; i < NTHREADS; i++)
      Tstack[i] = i;
    Tavail = NTHREADS;

    pthread_mutex_init(&TMUTEX,NULL);
    pthread_cond_init(&TCOND,NULL);

    done   = -NTHREADS;
    launch = 0;
    for (i = 0; i < gdb->ncontig; i++)
      { pthread_mutex_lock(&TMUTEX);

        if (Tavail <= 0)                       //  all threads are busy, wait
          pthread_cond_wait(&TCOND,&TMUTEX);

        tid = Tstack[--Tavail];                //  thread tid is available

        pthread_mutex_unlock(&TMUTEX);

        done   += 1;
        launch += 1;

        // Launching job for contig i on thread tid

        parm[tid].over->aread = i;

        if (VERBOSE)
          { if (done >= 0)
              fprintf(stderr,"\r  Launched %3d  Finished %3d",launch,done);
            else
              fprintf(stderr,"\r  Launched %3d  Finished   0",launch);
            fflush(stdout);
          }

        pthread_create(threads+tid,NULL,compress_thread,parm+tid);
      }

#ifndef DEBUG_THREADS
    pthread_mutex_lock(&TMUTEX);   //  Wait for all the jobs to complete
    while (Tavail < NTHREADS)
      { pthread_cond_wait(&TCOND,&TMUTEX);
        done += 1;
        if (VERBOSE)
          { fprintf(stderr,"\r  Launched %3d  Finished %3d",gdb->ncontig,done);
            fflush(stdout);
          }
      }
    pthread_mutex_unlock(&TMUTEX);
#endif

   if (VERBOSE)
     { fprintf(stderr,"\n");
       fflush(stderr);
     }

    for (i = 0; i < NTHREADS; i++)
      { free(parm[i].trace);
        if (MAKE_MASK)
          { free(parm[i].mstring);
            free(parm[i].mlist);
          }
        free(parm[i].buffer-1);
        free(parm[i].block);
        if (i == 0)
          Free_Align_Spec(parm[i].spec);
        Free_Work_Data(parm[i].work);
      }

    if (MAKE_MASK)
      { oneFileClose(Mfile);
        oneSchemaDestroy(anoSchema);
      }
    oneFileClose(Ofile);

    Close_GDB(gdb);

    if (VERBOSE)
      { TimeTo(stderr,0,1);
        TimeTo(stderr,1,0);
      }

    free(spath);

    Catenate(NULL,NULL,NULL,NULL);
    Numbered_Suffix(NULL,0,NULL);
    free(Prog_Name);

    exit (0);
  }
}
