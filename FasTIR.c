/*******************************************************************************************
 *
 *  Search an assembly for Terminal Inverted Repeats
 *
 *  Author:   Chenxi Zhou (Adapted from Gene Myers' FastLTR.c)
 *  Creation: May 2026
 *  Last Mod: May 2026
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

#undef   PROLOG
#undef   SORT1
#undef   SORT2
#undef   SHOW_SEEDS
#undef   SHOW_SEARCH
#undef   SHOW_ALIGNMENTS

#define BLOCK_OVERLAP   20000  //  Must be less than 32000

#define MINARM             20
#define MINGAP             50

#define CHAIN_MAXSKIP      20   //  max bpos/apos step between chained seeds
#define CHAIN_INDEL_TOL     4   //  max gap change within a chain
#define CHAIN_MIN_SCORE    40   //  min chain score (~5 seeds) to launch alignment
#define MAX_KMER_OCC      100   //  skip over-represented k-mers

#define TSPACE     100
#define VERSION  "0.5"

static char *Usage = "[-v] [-T(8)] <source:path>[<fa_extn>|<1_extn>] <target>[.1aln]";

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
     uint8      *aseq;     //  will point to contig buffer + off
     uint8      *bseq;     //  will fill dynamically 32kb blocks

     void       *block;    //  memory block for diagonal analyzer (128kb)
     void       *seeds;    //  memory block for seeds (dynamically allocated)
     int64       ssize;    //  size of seeds block in # seeds

     int         tmax;     //  i64 trace vector for 1-file
     int64      *trace;
   } S_Bundle;

typedef struct
  { uint16  apos;
    uint16  bpos;
    int     prev;
    int     score;
  } Seed;

typedef struct
  { int which;
    int score;
  } Sord;

static int CSCORE(const void *a, const void *b)
{ const Sord *sa = (const Sord *) a;
  const Sord *sb = (const Sord *) b;
  return (sb->score - sa->score);
}

static int spectrum_block(S_Bundle *bundle, int len, int off)
{ Alignment  *align = bundle->align;
  Overlap    *over  = bundle->over;
  Align_Spec *spec  = bundle->spec;
  Work_Data  *work  = bundle->work;
  OneFile    *ofile = bundle->ofile;
  GDB        *gdb   = bundle->gdb;
  int64      *t64   = bundle->trace;
  int         tmax  = bundle->tmax;

  uint8      *aseq  = (uint8 *) align->aseq;
  uint8      *bseq  = (uint8 *) align->bseq;
  int         clen  = gdb->contigs[over->aread].clen;

  int         i, j, p, x, c, f;
  uint16      kmer;
  uint32     *pairs;  // 0x10000
  uint16     *count;  // 0x10000
  uint16     *index;  // 0x08000
  Seed       *hits;   // 0x10000
  Sord       *sord;   // 0x10000
  uint8      *a7, *b7;
  int         l7;
  
  if (len < 8) return (0);

  pairs = (uint32 *) bundle->block;
  count = (uint16 *) (pairs + 0x10000);
  index = (uint16 *) (count + 0x10000);
  hits  = (Seed *)   bundle->seeds;
  sord  = (Sord *)   (hits + bundle->ssize);

#ifdef PROLOG
  printf("\nPANEL %d-%d\n",off,off+0x8000);
  fflush(stdout);
#endif

  for (i = 0; i < 0x10000; i++)   //  Init counters
    count[i] = 0;

  l7 = len-7;
  a7 = aseq+7;
  b7 = bseq+7;

  { kmer = a7[-7];               //  count # of each 8-mer in aseq
    for (i = -6; i < 0; i++)
      kmer = (kmer << 2) | a7[i];
    for (i = 0; i < l7; i++)
      { kmer = (kmer << 2) | a7[i];
        count[kmer] += 1;
      }
  }

  p = 0;                         //  turn counts into ptrs
  for (i = 0; i < 0x10000; i++)
    { x = count[i];
      count[i] = p;
      p += x;
    }

  { kmer = a7[-7];               //  place aseq positions in index sorted by 8-mer
    for (i = -6; i < 0; i++)
      kmer = (kmer << 2) | a7[i];
    for (i = 0; i < l7; i++)
      { kmer = (kmer << 2) | a7[i];
        index[count[kmer]++] = i;
      }
  }

  for (i = 0; i < 0x10000; i++)   //  init seed pairs counters
    pairs[i] = 0;
  
  { int b, e, m;                //  fist pass to count seed pairs
    kmer = b7[-7];
    for (i = -6; i < 0; i++)
      kmer = (kmer << 2) | b7[i];
    for (i = 0; i < l7; i++)
      { kmer = (kmer << 2) | b7[i];
        b = (kmer == 0) ? 0 : count[kmer-1];
        e = count[kmer];
        if (e - b > MAX_KMER_OCC) continue;
        for (f = b, m = len - 16 - i; f < e; f++)
          { j = index[f];
            if (j < i && j <= m)
              pairs[j] += 1;
          }
      }
  }

  p = 0;                          //  turn counts into ibeg sort ptrs
  for (i = 0; i < 0x8000; i++)
    { x = pairs[i];
      pairs[i] = p;
      p += x;
    }

  if (p > (int)bundle->ssize)    //  ensure enough space for p seeds
    { int64 ssnew = (int64)(1.2*p) + 0x8000;
      bundle->seeds = realloc(bundle->seeds, ssnew*20); // sizeof(Seed) + sizeof(Sord)
      if (bundle->seeds == NULL)
        { fprintf(stderr,"%s: Not enough memory\n",Prog_Name);
          exit (1);
        }
      bundle->ssize = ssnew;
      hits = (Seed *) bundle->seeds;
      sord = (Sord *) (hits + ssnew);
    }

  { int b, e, m;                //  second pass to collect seed pairs
    kmer = b7[-7];
    for (i = -6; i < 0; i++)
      kmer = (kmer << 2) | b7[i];
    for (i = 0; i < l7; i++)
      { kmer = (kmer << 2) | b7[i];
        b = (kmer == 0) ? 0 : count[kmer-1];
        e = count[kmer];
        if (e - b > MAX_KMER_OCC) continue;
        for (f = b, m = len - 16 - i; f < e; f++)
          { j = index[f];
            if (j < i && j <= m)
              { c = pairs[j]++;
                hits[c].apos = j;
                hits[c].bpos = i;
              }
          }
      }
  }

#ifdef SHOW_SEEDS
  { int a, b;
    for (i = 0; i < p; i++)
      { a  = hits[i].apos;
        b  = hits[i].bpos;
        printf("  gap=%6d a=%6d b=%6d  aseq:",b-a,a,b);
        Print_Seq(aseq + a, 8);
        printf("  bseq:");
        Print_Seq(bseq + b, 8);
        printf("\n");
        fflush(stdout);
      }
  }
#endif

  { int   outhit, nchain, score;
    int   s, da, db, rt, bp, bs, dgap;
    int   diag, wide, anti;
    int   abpos, aepos, bbpos, bepos, boff;
    uint8 saved;
    Path *path;

    //  Chaining DP: hits[0..p-1] already sorted by apos (aseq position)

    for (i = 0; i < p; i++)
      { bs = 8; bp = -1;
        for (j = i-1; j >= 0; j--)
          { da = (int)hits[i].apos - (int)hits[j].apos;
            if (da > CHAIN_MAXSKIP) break;
            if (da <= 0) continue;
            db = (int)hits[i].bpos - (int)hits[j].bpos;
            if (db <= 0 || db > CHAIN_MAXSKIP) continue;
            dgap = ((int)hits[i].apos - (int)hits[i].bpos)
                 - ((int)hits[j].apos - (int)hits[j].bpos);
            if (dgap < -CHAIN_INDEL_TOL || 
                dgap >  CHAIN_INDEL_TOL) 
              continue;
            score = hits[j].score + 8;
            if (score > bs)
              { bs = score;
                bp = j;
              }
          }
        hits[i].score = bs;
        hits[i].prev  = bp;
      }

    //  Sort seeds by score descending so best chain heads come first
    for (i = 0; i < p; i++)
      { sord[i].which = i;
        sord[i].score = hits[i].score;
      }
    qsort(sord, p, sizeof(Sord), CSCORE);

    //  Greedily claim chains in score order
    nchain = 0;
    for (i = 0; i < p; i++)
      { rt = sord[i].which;
        if (hits[rt].score < CHAIN_MIN_SCORE) break;
        score = 8;
        bp = rt;
        while ((s = hits[bp].prev) >= 0)
          { if (hits[s].score < 0)
              { //  already claimed by better chain
                hits[bp].prev = -1; //  break chain here
                break; 
              }
            score += 8;
            bp = s;
          }
        hits[rt].score = score;  //  save head score
        if (score < CHAIN_MIN_SCORE)
          continue;
        s = rt;
        while ((s = hits[s].prev) >= 0)
          hits[s].score = -1;  //  claim chain
        sord[nchain++] = (Sord) {rt, score};
      }

    //  Resort selected chains by score descending
    qsort(sord, nchain, sizeof(Sord), CSCORE);

    //  Align from each valid chain head in score order
    // work in the original coordinate space of the contig
    boff = clen - len - off;
    align->aseq -= off;
    align->bseq -= boff;
    saved = aseq[len]; aseq[len] = 4;
    outhit = 0;
    for (i = 0; i < nchain; i++)
      { rt = sord[i].which;

        //  Root = earliest seed in chain (lowest bpos = index bp)
        bp = rt;
        while (hits[bp].prev >= 0) bp = hits[bp].prev;

#ifdef SHOW_SEARCH
        printf("  search chain: score=%d  head=(a=%d,b=%d,gap=%d)  root=(a=%d,b=%d,gap=%d)\n",
               hits[bp].score,
               (int)hits[rt].apos, (int)hits[rt].bpos,
               (int)hits[rt].apos - (int)hits[rt].bpos,
               (int)hits[bp].apos, (int)hits[bp].bpos,
               (int)hits[bp].apos - (int)hits[bp].bpos);
        fflush(stdout);
#endif

        if (aseq[(int)hits[bp].apos] >= 4) continue;  //  already masked by another TIR

        diag = (int)hits[bp].apos - (int)hits[bp].bpos + off - boff;
        anti = (int)hits[bp].apos + (int)hits[bp].bpos + off + boff;
        wide = (int)hits[rt].apos - (int)hits[bp].apos;
        wide = (int)(wide*0.2) + CHAIN_INDEL_TOL;
        
        Local_Alignment(align, work, spec, diag, diag, anti, wide, wide);

        path = align->path;
        //  sequence chunk coordinates for TIR selection and masking
        abpos = path->abpos - off;
        aepos = path->aepos - off;
        bbpos = len - (path->bepos - boff);
        bepos = len - (path->bbpos - boff);

        if (aepos - abpos + bepos - bbpos >= 2*MINARM && bbpos - aepos >= MINGAP)
          { if (over->path.tlen > tmax)
              { tmax = bundle->tmax = (int)(1.2*over->path.tlen) + 1000;
                t64  = bundle->trace = (int64 *) realloc(t64, sizeof(int64)*tmax);
              }
            Write_Aln_Overlap(ofile, over);
            Compress_TraceTo8(over, 0);
            Write_Aln_Trace(ofile, over->path.trace, over->path.tlen, t64, 0);

#ifdef SHOW_ALIGNMENTS
            printf("\nTIR: diag=%d anti=%d wide=%d\n", diag, anti, wide);
            printf("  left  arm: aseq[%d..%d)\n", abpos+off, aepos+off);
            printf("  right arm: aseq[%d..%d)\n", bbpos+off, bepos+off);
            align->alen = align->blen = clen; //  need to be in original coordinates for Compute_Trace_PTS
            Decompress_TraceTo16(over);
            Compute_Trace_PTS(align, work, TSPACE, GREEDIEST, 1, -1);
            Print_Reference(stdout, align, work, 8, 100, 10, 0, 10, 0);
            align->alen = align->blen = len;  // restore for next round
            fflush(stdout);
#endif

            for (x = abpos; x < aepos; x++)
              aseq[x] = 4;
            for (x = bbpos; x < bepos; x++)
              aseq[x] = 4;
            if (outhit < bepos)
              outhit = bepos;
          }
      }
    //  Restore coordinates
    aseq[len] = saved;
    align->aseq += off;
    align->bseq += clen - off - len;
    
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
  char     *buffer = (char *) bundle->buffer;
  char     *bseq   = (char *) bundle->bseq;
  GDB      *gdb    = bundle->gdb;
  int       last, clen, blen;
  int       i, p;
 
  i = bundle->over->aread;

  Get_Contig(gdb,i,NUMERIC,(char *) buffer);
#ifdef PROLOG
  printf("CONTIG %d\n",i+1);
#endif

  clen = gdb->contigs[i].clen;
  
  last = -1;
  if (clen < 0x8000) {
    bundle->align->aseq = (char *) buffer;
    bundle->align->alen = bundle->align->blen = clen;
    memcpy(bseq,buffer,clen);
    Complement_Seq(bseq,clen);
    ((uint8 *) bseq)[clen] = 4;
    spectrum_block(bundle, clen, 0);
  } else
    for (p = 0; p+0x2000 <= clen; p += BLOCK_OVERLAP)
      { if (p+0x8000 > clen)
          { blen = clen-p;
            bundle->align->aseq = (char *) buffer+p;
            bundle->align->alen = bundle->align->blen = blen;
            memcpy(bseq,buffer+p,blen);
            Complement_Seq(bseq,blen);
            ((uint8 *) bseq)[blen] = 4;
            ((uint8 *) buffer)[p-1] = 4;
            spectrum_block(bundle, blen, p);
            break;
          }
        else
          { blen = 0x8000;
            bundle->align->aseq = (char *) buffer+p;
            bundle->align->alen = bundle->align->blen = blen;
            memcpy(bseq,buffer+p,blen);
            Complement_Seq(bseq,blen);
            ((uint8 *) bseq)[blen] = 4;
            ((uint8 *) buffer)[p-1] = 4;
            last = spectrum_block(bundle, blen, p);
            if (last > blen - BLOCK_OVERLAP)
              p += last - (blen - BLOCK_OVERLAP);
          }
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

  (void) Print_Seq;
  (void) emer;

  //   Process command line

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("FastTIR")

    NTHREADS = 8;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'T':
            ARG_NON_NEGATIVE(NTHREADS,"number of threads to use");
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE   = flags['v'];

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        exit (1);
      }
  }

  //  Get GDB or make a temporary if a fasta

  { char *cpath, *APATH, *AROOT;
    
    (void) Get_GDB_Paths(argv[1],NULL,&spath,&cpath,0);

    free(cpath);
  
    units = Get_GDB(gdb,spath,".",NTHREADS,NULL);

    //  Open 1aln file for threaded writing

    APATH = PathTo(argv[2]);
    AROOT = Root(argv[2],".1aln");
    cpath = getcwd(NULL,0);

    Ofile = open_Aln_Write(Catenate(APATH,"/",AROOT,".1aln"),NTHREADS,Prog_Name,VERSION,
                           Command_Line,TSPACE,spath,NULL,cpath);

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
        parm[i].gdb  = &parm[i]._gdb;
        parm[i]._gdb  = _gdb;
        parm[i]._gdb.seqs = units[i];
        parm[i].work  = New_Work_Data();
        if (i == 0)
          parm[i].spec = New_Align_Spec(.7,TSPACE,gdb->freq,0);
        else
          parm[i].spec = parm[i-1].spec;
        parm[i].block = malloc(8*0x10000);   // 512kb
        parm[i].buffer = ((uint8 *) malloc(gdb->maxctg + 4)) + 1;
        parm[i].buffer[-1] = 4;              // sentinel for Local_Alignment reverse wave
        parm[i].seeds = malloc(20*0x10000);  // initial seeds block
        parm[i].ssize = 0x10000;             //  max # seeds in block
        parm[i].bseq   = ((uint8 *) malloc(0x8000 + 4)) + 1;
        parm[i].bseq[-1] = 4;               // sentinel for Local_Alignment reverse wave
        parm[i].over   = &parm[i]._over;
        parm[i].over->flags  = COMP_FLAG;  // bseq is reverse complement
        parm[i].align = &parm[i]._align;
        parm[i].align->path  = &(parm[i]._over.path);
        parm[i].align->flags = COMP_FLAG;  // bseq is reverse complement
        parm[i].align->bseq  = (char *) parm[i].bseq;
        parm[i].tmax    = 10000;
        parm[i].trace   = malloc(sizeof(int64)*10000);

        if (parm[i].block == NULL || parm[i].trace == NULL || 
            parm[i].buffer == NULL || parm[i].bseq == NULL)
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
        parm[tid].over->bread = i;

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
      { free(parm[i].trace);
        free(parm[i].buffer-1);
        free(parm[i].bseq-1);
        free(parm[i].block);
        if (i == 0)
          Free_Align_Spec(parm[i].spec);
        Free_Work_Data(parm[i].work);
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
