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

#define TIR_MIN    300   //  Minimum length of interior segment
#define TIR_MAX  16000   //  Includes TIRs, must be less than 24K, time depends critically
                         //     on this as it determines the amount of overlap between blocks
#define TSPACE     100
#define VERSION  "0.5"

static char *Usage = "[-v] [-T(8)] [-o<target>] <source:path>[<fa_extn>|<1_extn>|.1gdb]";

static int NTHREADS;
static int VERBOSE;

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
 *  ANTI-DIAGONAL HITS DETECTOR
 *
 ********************************************************************************************/

typedef struct
  { uint16  anti;
    uint16  ibeg;
  } Seed;

typedef struct
  { uint16 anti;
    uint16 count;
  } Chord;

  //  Thread bundle

typedef struct
   { int         tid;
     OneFile    *ofile;
     GDB        _gdb, *gdb;

     Work_Data  *work;     //  alignment machinery
     Align_Spec *spec;
     Overlap     _over, *over;
     Alignment   _align, *align;

     uint8      *forseq;   //  contig buffer
     uint8      *revseq;   //  contig buffer

     void       *block;    //  memory block for anti-diagonal analyzer (524MB)
     int         pmax;
     Seed       *posts;

     int         tmax;     //  i64 trace vector for 1-file
     int64      *trace;
   } S_Bundle;

static int CSORT(const void *l, const void *r)
{ Chord *x = (Chord *) l;
  Chord *y = (Chord *) r;

  return (y->count - x->count);
}

static uint16 comph[256];
static uint16 compl[256];

static void spectrum_block(uint8 *seq, int off, int len, S_Bundle *bundle)
{ Alignment  *align = bundle->align;
  Overlap    *over  = bundle->over;
  Align_Spec *spec  = bundle->spec;
  Work_Data  *work  = bundle->work;
  OneFile    *ofile = bundle->ofile;
  int64      *t64   = bundle->trace;
  int         tmax  = bundle->tmax;
  Seed       *post  = bundle->posts;
  int         pmax  = bundle->pmax;

  (void) align;
  (void) over;
  (void) spec;
  (void) work;
  (void) ofile;
  (void) t64;
  (void) tmax;
  (void) off;

  uint16  kmer, hmer;
  uint16 *index;  // 0x08000
  uint16 *count;  // 0x10001
  int    *antis;  // 0x10000

  Seed   *hits;   // dynamic
  int     npairs;

  uint8  *s7;
  int     l7;
  int     i, p, x;

  antis = (int *) bundle->block;
  index = (uint16 *) (antis + 0x10000);
  count = index + 0x08000; 

#ifdef PROLOG
  printf("\nPANEL %d-%d\n",off,off+0x8000);
  fflush(stdout);
#endif

  for (i = 0; i < 0x10000; i++)   //  Init counters
    antis[i] = count[i] = 0;

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

  for (i = 0xffff; i > 0; i--)   //  reset count ptrs
    count[i] = count[i-1];
  count[0] = 0;
  count[0x10000] = l7;

#ifdef SORT1
  x = 0;
  for (i = 0; i < l7; i++)
    { p = index[i];
      printf("%5d: ",p);
      Print_Seq(seq+p,8);
      if (count[x] == i)
        { printf(" +");
          while (count[x] == i)
            x += 1;
        }
      printf("\n");
    }
#endif

  npairs = 0;
  kmer = seq[0];            //  Record # of anti-pairs at least TIR_MIN apart in sorted order
  for (i = 1; i < 7; i++)   //    of 1st pair pos, prepare to radix sort on anti-diagonal
    kmer = (kmer << 2) | seq[i];
  for (i = 0; i < l7; i++)
    { kmer = (kmer << 2) | s7[i];
      hmer = comph[kmer&0xff] | compl[kmer>>8];
      for (p = count[hmer+1]-1; p >= count[hmer]; p--)
        { x = index[p];
          if (x < i+TIR_MIN)
            break;
          if (npairs >= pmax)
            { pmax = 1.2*npairs + 1000;
              post = realloc(post,2*sizeof(Seed)*pmax);
              if (post == NULL)
                { fprintf(stderr,"%s: Out of memory realloc'ing seeds\n",Prog_Name);
                  exit (1);
                }
            }
          x += i;
          post[npairs].ibeg = i;
          post[npairs].anti = x;
          antis[x] += 1;
          npairs += 1;
        }
    }

  p = 0;                          //  turn counts into ptrs
  for (i = 0; i < 0x10000; i++)
    { x = antis[i];
      antis[i] = p;
      p += x;
    }

  hits = post + npairs;

#ifdef SORT2
  printf("\nSorted on 1st pos\n");
  for (i = 0; i < p; i++)
    printf(" %5d %5d\n",post[i].ibeg,post[i].anti);
#endif

  for (i = 0; i < p; i++)       //  place ibeg/anti pairs in hits in order of anti then ibeg
    { x = post[i].anti;
      hits[antis[x]++] = post[i];
    }

#ifdef SHOW_SEEDS
  printf("\nSorted on Anti, then 1st pos\n");
  for (i = 0; i < npairs; i++)
    printf("   %4d : %5d\n",hits[i].anti,hits[i].ibeg);
#endif

  { int    ncnt;
    Chord *hist = (Chord *) post;
    int    b, c, cover;

    ncnt = 0;
    p = antis[0];
    for (i = 1; i < 0x10000; i++)
      { x = antis[i];
        c = -8;
        cover = 0;
        for ( ; p < x; p++)
          { b = hits[p].ibeg;
            c = b-c;
            if (c >= 8)
              cover += 8;
            else
              cover += c;
            c = b;
          }
        if (cover >= 24)
          { hist[ncnt].count = cover;
            hist[ncnt].anti  = i;
            ncnt += 1;
          }
        p = x;
      }

    qsort(hist,ncnt,sizeof(Chord),CSORT);

    printf("\nBuckets:\n");
    for (i = 0; i < ncnt; i++) 
      { x = hist[i].anti;
        printf("%5d (%d):\n",x,hist[i].count);
        for (p = antis[x-1]; p < antis[x]; p++)
          printf(" %5d [%5d]\n",hits[p].ibeg,hits[p].anti);
      }
  }

  bundle->pmax  = pmax;
  bundle->posts = post;
  return;
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
  uint8    *forseq = bundle->forseq;
  uint8    *revseq = bundle->revseq;
  GDB      *gdb    = bundle->gdb;
  int       ctg, clen;
  int       p, q;
 
  ctg = bundle->over->aread;

  clen = gdb->contigs[ctg].clen;
  Get_Contig(gdb,ctg,NUMERIC,(char *) forseq);

  q = clen-1;
  for (p = 0; p < clen; p++)
    revseq[p] = 3-forseq[q--];
  revseq[-1] = revseq[clen] = 4;
  
#ifdef PROLOG
  printf("CONTIG %d\n",ctg+1);
#endif

  bundle->align->aseq = (char *) forseq;
  bundle->align->bseq = (char *) revseq;
  bundle->align->alen = bundle->align->blen = clen;
  bundle->over->bread = ctg;

  if (clen < 0x8000)
    spectrum_block(forseq,0,clen,bundle);
  else
    for (p = 0; p < clen; p += TIR_MAX)
      { if (p+0x8000 >= clen)
          { spectrum_block(forseq+p,p,clen-p,bundle);
            break;
          }
        else
          spectrum_block(forseq+p,p,0x8000,bundle);
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
  char      *TRGT_PATH;
  OneFile   *Ofile;

  (void) Print_Seq;
  (void) emer;

  //   Process command line

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("FastTIR")

    NTHREADS  = 8;
    TRGT_PATH = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'o':
            TRGT_PATH = argv[i]+2;
            break;
          case 'T':
            ARG_NON_NEGATIVE(NTHREADS,"number of threads to use");
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE   = flags['v'];

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -o: Root path of .1aln/.1ano file (default root path of input).\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        exit (1);
      }
  }

  //  Get GDB or make a temporary if a fasta

  { char *cpath, *APATH, *AROOT;

    Get_GDB_Paths(argv[1],NULL,&spath,&cpath,0);

    free(cpath);

    units = Get_GDB(gdb,spath,".",NTHREADS,NULL);

    //  Get target root path

    if (TRGT_PATH == NULL)
      { AROOT = Root(spath,NULL);
        APATH = PathTo(spath);
      }
    else
      { if (strcmp(TRGT_PATH + (strlen(TRGT_PATH)-4), ".1aln") ||
            strcmp(TRGT_PATH + (strlen(TRGT_PATH)-4), ".1ano"))
          AROOT = Root(TRGT_PATH,NULL);
        else
          AROOT = Root(TRGT_PATH,"");
        APATH = PathTo(TRGT_PATH);
      }

    cpath = getcwd(NULL,0);

    //  Open 1aln file for threaded writing

    Ofile = open_Aln_Write(Catenate(APATH,"/",AROOT,".1aln"),NTHREADS,Prog_Name,VERSION,
                           Command_Line,TSPACE,spath,NULL,cpath);

    Write_Skeleton(Ofile,gdb);

    free(cpath);
    free(AROOT);
    free(APATH);
  }

  if (VERBOSE)
    { fprintf(stderr,"\n  Database loaded, begin scan of %d contigs\n\n",gdb->ncontig);
      fflush(stderr);
    }

  StartTime();

  //  Establish mapping tables

  { int i, i1, i2, i3, i4;

    i = 0;
    for (i1 = 3; i1 >= 0; i1 -= 1)
     for (i2 = 12; i2 >= 0; i2 -= 4)
      for (i3 = 48; i3 >= 0; i3 -= 16)
       for (i4 = 192; i4 >= 0; i4 -= 64)
         { compl[i] =  i4 | i3 | i2 | i1;
           comph[i] = (compl[i] << 8);
           i += 1;
         }
  }

  { int       i, tid;
    int       done, launch;
    pthread_t threads[NTHREADS];
    S_Bundle  parm[NTHREADS];
    int       tstack[NTHREADS];

    for (i = 0; i < NTHREADS; i++)
      { parm[i].tid   = i;
        parm[i].ofile = Ofile + i;
        parm[i].gdb   = &parm[i]._gdb;
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
        parm[i].tmax    = 10000;
        parm[i].trace   = malloc(sizeof(int64)*10000);

        parm[i].block   = malloc(7*0x10000+2);   // 448KB

        parm[i].forseq  = ((uint8 *) malloc(gdb->maxctg + 4)) + 1;
        parm[i].revseq  = ((uint8 *) malloc(gdb->maxctg + 4)) + 1;

        parm[i].pmax    = 0x40000;
        parm[i].posts   = malloc(2*sizeof(Seed)*0x40000);

        if (parm[i].block == NULL || parm[i].forseq == NULL || parm[i].revseq == NULL
                   || parm[i].trace == NULL || parm[i].posts == NULL)
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
        pthread_detach(threads[tid]);
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
      { free(parm[i].revseq-1);
        free(parm[i].forseq-1);
        free(parm[i].block);
        free(parm[i].trace);
        if (i == 0)
          Free_Align_Spec(parm[i].spec);
        else
          fclose(units[i]);
        Free_Work_Data(parm[i].work);
      }

    oneFileClose(Ofile);

    if (NTHREADS > 1)
      free(units);
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
