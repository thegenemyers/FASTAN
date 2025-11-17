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
#include "align.h"
#include "alncode.h"

#define VERSION "0.1"
#define TSPACE   100
#define DIAG_MAX 8000

static char *Usage = "[-vm] [-T(8)] <source:path>[<fa_extn>|<1_extn>] <target>[.1aln]";

static int NTHREADS;
static int VERBOSE;
static int MASK_DB;

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
 *  DYNAMIC MASK LISTS
 *
 ********************************************************************************************/

typedef struct _Mask_Block
  { struct _Mask_Block *link;
    int64               beg[1000]; 
    int64               end[1000]; 
  } Mask_Block;

typedef struct _Mask_List
  { Mask_Block *start;
    Mask_Block *current;
    int         ctop;
    int         nblk;
  } Mask_List;

static void Init_Masks(Mask_List *masks)
{ masks->current = Malloc(sizeof(Mask_Block),"Allocating first Mask_Block");
  if (masks->current == NULL)
    exit (1);
  masks->start = masks->current;
  masks->current->link = NULL;
  masks->ctop = 0;
  masks->nblk = 0;
}

static void Add_Mask(int beg, int end, Mask_List *masks)
{ int         ctop;
  Mask_Block *curr, *next;
 
  ctop = masks->ctop;
  curr = masks->current;
  if (ctop >= 1000)
    { next = curr->link;
      if (next == NULL)
        { curr->link = next = Malloc(sizeof(Mask_Block),"Allocating another Mask_Block");
          if (next == NULL)
            exit (1);
          next->link = NULL;
        }
      masks->current = curr = next;
      ctop  = 0;
      masks->nblk += 1;
    }
  curr->beg[ctop] = beg;
  curr->end[ctop] = end;
  masks->ctop = ctop+1;
}

static void Free_Masks(Mask_List *masks)
{ Mask_Block *m, *n;

  for (m = masks->start; m != NULL; m = n)
    { n = m->link;
      free(m);
    }
}


/*******************************************************************************************
 *
 *  SPECIAL VERSION OF GDB WRITE
 *
 ********************************************************************************************/

//  If the genome/sequence is masked then it must be upper-case (i.e. u-line is present)
    
static char *gdbSchemaText =
  "1 3 def 1 0                         schema for genome skeleton\n"
  ".\n"
  "P 3 gdb                             GDB\n"
  "D f 4 4 REAL 4 REAL 4 REAL 4 REAL   global: base frequency vector\n"
  "D u 0                               global: upper case when displayed\n"
  "O S 1 6 STRING                      id for a scaffold\n"
  "D G 1 3 INT                         gap of given length\n"
  "D C 1 3 INT                         contig of given length\n"
  "D M 1 8 INT_LIST                    mask pair list for a contig\n"
;         
          
static OneSchema *make_GDB_Schema()
{ return (oneSchemaCreateFromText(gdbSchemaText)); }


// Write the given gdb to the file 'tpath'.  The GDB must have seqstate EXTERNAL and tpath
//   must be consistent with the name of the .bps file.

extern bool addProvenance(OneFile *of, OneProvenance *from, int n) ; // backdoor - clean up some day

int Write_GDB_Special(GDB *gdb, char *tpath)
{ OneSchema *schema;                   
  OneFile   *of;                       
  bool       binary;                   
  int64      spos, len;
  char      *head;
  int        s, c;
  GDB_MASK  *m, *n;

  if (strcmp(tpath+strlen(tpath)-4,".gdb") == 0)
    binary = false;
  else 
    binary = true;

  schema = make_GDB_Schema();
  if (schema == NULL)
    { EPRINTF("Failed to create GDB schema (Write_GDB)");
      EXIT(1);
    }

  of = oneFileOpenWriteNew(tpath,schema,"gdb",binary,1);
  if (of == NULL)
    { EPRINTF("Failed to open GDB file %s (Write_GDB)",tpath);
      oneSchemaDestroy(schema);
      EXIT(1);
    }

  addProvenance(of,gdb->prov,gdb->nprov);
  oneAddProvenance(of,Prog_Name,"0.1",Command_Line);

  oneAddReference(of,gdb->srcpath,1);

  oneReal(of,0) = gdb->freq[0];
  oneReal(of,1) = gdb->freq[1];
  oneReal(of,2) = gdb->freq[2];
  oneReal(of,3) = gdb->freq[3];
  oneWriteLine(of,'f',0,0);
  oneWriteLine(of,'u',0,0);

  for (s = 0; s < gdb->nscaff; s++)
    { head = gdb->headers + gdb->scaffolds[s].hoff;
      oneWriteLine(of,'S',strlen(head),head);

     spos = 0;
      for (c = gdb->scaffolds[s].fctg; c < gdb->scaffolds[s].ectg; c++)
        { if (gdb->contigs[c].sbeg > spos)
            { oneInt(of,0) = gdb->contigs[c].sbeg - spos;
              oneWriteLine(of,'G',0,0);
            }
          len = gdb->contigs[c].clen;
          oneInt(of,0) = len;
          oneWriteLine(of,'C',0,0);
          spos = gdb->contigs[c].sbeg + len;
          
          m = n = (GDB_MASK *) (gdb->contigs[c].moff);
          while (n->beg >= 0)
            n += 1;
          len = n-m;
          if (len > 0)
            oneWriteLine(of,'M',len,(int64 *) m);
        }
      if (gdb->scaffolds[s].slen > spos)
        { oneInt(of,0) = gdb->scaffolds[s].slen - spos;
          oneWriteLine(of,'G',0,0);
        }
    }

  oneFileClose(of);
  oneSchemaDestroy(schema);
  return (0);
}


/*******************************************************************************************
 *
 *  SEED CHAIN DETECTOR
 *
 ********************************************************************************************/

#undef   PROLOG
#undef   SORT1
#undef   SORT2
#undef   SHOW_SEEDS
#undef   SHOW_SEARCH
#undef   SHOW_ALIGNMENTS
 
typedef struct
   { int         tid;
     OneFile    *ofile;
     Work_Data  *work;
     Align_Spec *spec;
     GDB         _gdb, *gdb;
     Overlap     _over, *over;
     Alignment   _align, *align;
     void       *block;
     uint8      *buffer;
     int         tmax;
     int64      *trace;
     Mask_List   masks;
   } S_Bundle;

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
  Mask_List  *masks = &bundle->masks;

  int     i, p, x, c;
  int     d, e, f;
  uint16  kmer;
  uint16 *index;  // 0x08000
  uint16 *count;  // 0x10000
  uint16 *diags;  // DIAG_MAX < 0x08000
  Seed   *post;   // 0x08000
  Seed   *hits;   // 0x08000
  uint8 *s7;
  int    l7;

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

  kmer = seq[0];                 //  count # of each 8-mer
  for (i = 1; i < 7; i++)
    kmer = (kmer << 2) | seq[i];
  for (i = 0; i < l7; i++)
    { kmer = (kmer << 2) | s7[i];
      count[kmer] += 1;
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

  for (i = 0; i < DIAG_MAX; i++)        //  init diagonal tube counters
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
    int   wide, anti, last;
    Path *bpath;

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
        last = -1;
#ifdef SHOW_SEARCH
        printf(" %4d: %5d\n",d,diags[d]-diags[d-1]);
#endif
        for (x = diags[d-1]+1; x < diags[d]; x++)
          { p = hits[x].ibeg;
#ifdef SHOW_SEARCH
            printf("  p = %d (%d %d %d\n",p,last,hits[x-1].ibeg,mark[p+1]);
#endif
            if (p < last || p - hits[x-1].ibeg > d || seq[p+1] >= 4)
              continue;
            wide = .2*d;
            if (wide < 1)
              wide = 1;
            anti = 2*(off + p) + d;
            bpath = Local_Alignment(align,work,spec,d,d,anti,wide,wide);
#ifdef SHOW_SEARCH
            if (bpath == NULL)
              printf("    NULL\n");
            else
              printf("    %d (%d)  %d-%d-%d\n",bpath->bepos-bpath->abpos,2*d,d-wide,d,d+wide);
#endif
            if (bpath == NULL)
              continue;

            end = bpath->bepos - off;
            beg = bpath->abpos - off;

            if (beg > p || end <= p)
              continue;
            if (end > last)
              last = end;
            if (end-beg < 1.95*d)
              continue;

            if (over->path.tlen > tmax)
              { tmax = bundle->tmax = 1.2*over->path.tlen + 1000;
                t64  = bundle->trace = realloc(t64,sizeof(int64)*tmax);
              }
            Write_Aln_Overlap(ofile,over);
            Compress_TraceTo8(over,0);
            Write_Aln_Trace(ofile,over->path.trace,over->path.tlen,t64);
            oneInt(ofile,0) = d;
            oneWriteLine(ofile,'U',0,0);

            if (MASK_DB)
              Add_Mask(bpath->abpos,bpath->bepos,masks);

#ifdef SHOW_ALIGNMENTS
#ifndef SHOW_SEARCH
            if (last < 0)
              printf(" %4d: %5d\n",d,diags[d]-diags[d-1]);
            printf("\n");
#endif
            Decompress_TraceTo16(over);
            printf(" Hit spans %d-%d\n",bpath->abpos,bpath->bepos);
            Compute_Trace_PTS(align,work,100,GREEDIEST,d-wide,d+wide);
            Print_Alignment(stdout,align,work,8,100,10,0,10,0);
#endif

            if (bpath->aepos < bpath->bbpos)
              { for (f = (bpath->bbpos-off)+1; f <= end; f++)
                  seq[f] = 4;
                end = bpath->aepos-off;
                for (f = beg+1; f <= end; f++)
                  seq[f] = 4;
              }
            else
              { for (f = beg+1; f <= end; f++)
                  seq[f] = 4;
              }

            if (bpath->bepos > outhit)
              outhit = bpath->bepos;
          }
      }

    return (outhit);
  }
}

static pthread_mutex_t TMUTEX;
static pthread_cond_t  TCOND;

//  Tstack[0..Tavail-1] is a stack of available threads at any moment.
//  It is always manipulated inside the mutex TMUTEX

static int *Tstack;
static int  Tavail;

//  for 1st k-mer byte range [beg,end), find each group of equal k-mers, then
//    overwrite to the bottom of the range.  K-mer payloads contain the prefix mask, count,
//    and then the contig/position location.

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

  bundle->masks.ctop = 0;
  bundle->masks.current = bundle->masks.start;

  last = -1;
  if (clen < 0x8000)
    spectrum_block(buffer,0,clen,bundle);
  else
    for (p = 0; p+0x2000 <= clen; p += 0x6000)
      { if (p+0x8000 > clen)
          spectrum_block(buffer+p,p,clen-p,bundle);
        else
          { last = spectrum_block(buffer+p,p,0x8000,bundle);
            if (last >= p+0x6000)
              p = last-0x6000;
          }
      }

  if (MASK_DB)
    { Mask_Block *m;
      GDB_MASK   *cm;
      int         c, t, ctop;

      cm = Malloc(sizeof(GDB_MASK)*(1000*bundle->masks.nblk + bundle->masks.ctop + 1),
                         "Allocating contig mask array");
      gdb->contigs[i].moff = (int64) cm;

      t = 0;
      for (m = bundle->masks.start; 1; m = m->link)
        { if (m == bundle->masks.current)
            ctop = bundle->masks.ctop;
          else
            ctop = 1000;
          for (c = 0; c < ctop; c++, t++)
            { cm[t].beg = m->beg[c];
              cm[t].end = m->end[c];
              // printf("M %d %d\n",m->beg[c],m->end[c]);
            }
          if (m == bundle->masks.current)
            break;
        }
      cm[t].beg = -1;
    }
 
  pthread_mutex_lock(&TMUTEX);   //  Put this thread back on the avail stack
    Tstack[Tavail++] = bundle->tid;
  pthread_mutex_unlock(&TMUTEX);

  pthread_cond_signal(&TCOND);   //  Signal a thread is available

  return (NULL);
}

int main(int argc, char *argv[])
{ FILE   **units;
  char    *spath;
  GDB     _gdb, *gdb = &_gdb;
  OneFile *Ofile;

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
            ARG_FLAGS("vm")
            break;
          case 'T':
            ARG_NON_NEGATIVE(NTHREADS,"number of threads to use");
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    MASK_DB = flags['m'];

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        fprintf(stderr,"      -m: Mask the data base with hits found.\n");
        exit (1);
      }
  }

  //  Get GDB or make a temporary if a fasta

  { char *cpath, *APATH, *AROOT;
    int   ftype;

    ftype = Get_GDB_Paths(argv[1],NULL,&spath,&cpath,0);

    free(cpath);

    if (MASK_DB && ftype != IS_GDB)
      { fprintf(stderr,"%s: The source must be a GDB when masking (-m) is on\n",Prog_Name);
        exit (1);
      }
  
    units = Get_GDB(gdb,spath,".",NTHREADS);

    //  Open 1aln file for threaded writing

    APATH = PathTo(argv[2]);
    AROOT = Root(argv[2],".1aln");
    cpath = getcwd(NULL,0);

    Ofile = open_Aln_Write(Catenate(APATH,"/",AROOT,".1aln"),NTHREADS,Prog_Name,VERSION,
                           Command_Line,TSPACE,spath,NULL,cpath);
    free(cpath);
    free(AROOT);
    free(APATH);

    // oneAddProvenance(Ofile,Prog_Name,"v0.5",Command_Line);

    Write_Aln_Skeleton(Ofile,gdb);
  }

  if (VERBOSE)
    { fprintf(stderr," Database built, begin scan\n");
      fflush(stderr);
    }

  StartTime();

  { int       i, tid;
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
        parm[i].block  = malloc(9*0x10000);   // 576KB
        parm[i].buffer = ((uint8 *) malloc(gdb->maxctg + 4)) + 1;
        parm[i].tmax   = 10000;
        parm[i].trace  = malloc(sizeof(int64)*10000);
        if (parm[i].block == NULL || parm[i].buffer == NULL || parm[i].trace == NULL)
          { fprintf(stderr,"%s: Not enough memory\n",Prog_Name);
            exit (1);
          }
        Init_Masks(&parm[i].masks);
      }

    Tstack = tstack;
    for (i = 0; i < NTHREADS; i++)
      Tstack[i] = i;
    Tavail = NTHREADS;

    pthread_mutex_init(&TMUTEX,NULL);
    pthread_cond_init(&TCOND,NULL);

    for (i = 0; i < gdb->ncontig; i++)
      { pthread_mutex_lock(&TMUTEX);

        if (Tavail <= 0)                       //  all threads are busy, wait
          pthread_cond_wait(&TCOND,&TMUTEX);

        tid = Tstack[--Tavail];                //  thread tid is available

        pthread_mutex_unlock(&TMUTEX);

        // Launching job for contig i on thread tid

// printf("Launching %d\n",i); fflush(stdout);

        parm[tid].over->aread = i;

        pthread_create(threads+tid,NULL,compress_thread,parm+tid);
      }

#ifndef DEBUG_THREADS
    pthread_mutex_lock(&TMUTEX);   //  Wait for all the jobs to complete
    while (Tavail < NTHREADS)
      pthread_cond_wait(&TCOND,&TMUTEX);
    pthread_mutex_unlock(&TMUTEX);
#endif

    for (i = 0; i < NTHREADS; i++)
      { free(parm[i].trace);
        free(parm[i].buffer-1);
        free(parm[i].block);
        if (i == 0)
          Free_Align_Spec(parm[i].spec);
        Free_Work_Data(parm[i].work);
        Free_Masks(&parm[i].masks);
      }

    oneFileClose(Ofile);

    if (MASK_DB)
      Write_GDB_Special(gdb,spath);

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
