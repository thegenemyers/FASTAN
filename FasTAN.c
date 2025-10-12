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

#define TSPACE   100
#define VERSION "0.1"

static char *Usage = "[-vm] [-T(8)] <source:path>[<fa_extn>|<1_extn>] <target>[.1aln]";

static int NTHREADS;
static int VERBOSE;
static int MODELS;

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
   } S_Bundle;

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
 *  SEED CHAIN DETECTOR
 *
 ********************************************************************************************/

#define DIAG_MAX 8000

#undef   PROLOG
#undef   SORT1
#undef   SORT2
#undef   SHOW_SEEDS
#undef   SHOW_CHAINS
#undef   SHOW_SEARCH
#undef   SHOW_ALIGNMENTS
#undef   COMPUTE_MODEL

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

#ifdef COMPUTE_MODEL

typedef struct
  { int    edge[4];
    int    count;
    uint16 kmer;
    uint16 mark;
    uint16 depth;
  } Dnode;

static int compute_model(uint8 *seq, int len, uint16 *alive, int unit)
{ int     i, x, c, e;
  uint16  kmer, klen, umask;
  uint8  *s7, *s8;
  int     l7;
  int     clen;

  int    stop, node;
  Dnode *stack;
  int    maxnode, maxcount;
  double density;
  int    dcut;
  int   *mark;

  if (unit > 8)
    klen = 8;
  else
    klen = unit;

  s7 = seq+(klen-1);
  s8 = seq+klen;
  l7 = len-(klen-1);
  if (klen >= 8)
    umask = 0xffff; 
  else
    umask = (1 << (2*klen)) - 1;

  stop = 0;         //  count # of each 8-mer
  kmer = seq[0];
  for (i = 1; i < klen-1; i++)
    { x = seq[i];
      kmer = (kmer << 2) | x;
    }
  for (i = 0; i < l7; i++)
    { x = s7[i];
      kmer = ((kmer << 2) | x) & umask;
      node = alive[kmer];
      if (node == 0)
        { alive[kmer] = 0x7fff;
          stop += 1;
        }
    }

  stack = malloc(sizeof(Dnode)*stop);
  mark  = malloc(sizeof(int)*stop);
  if (stack == NULL || mark == NULL)
    { printf("MALLOC\n");
      exit (1);
    }

  maxcount = 0;
  stop = 0;         //  count # of each 8-mer
  kmer = seq[0];
  for (i = 1; i < klen-1; i++)
    { x = seq[i];
      kmer = (kmer << 2) | x;
    }
  for (i = 0; i < l7; i++)
    { x = s7[i];
      kmer = ((kmer << 2) | x) & umask;
      node = alive[kmer];
      if (node == 0x7fff)
        { alive[kmer] = node = stop++;
          stack[node].edge[0] = stack[node].edge[1] = stack[node].edge[2] = stack[node].edge[3] = 0;
          stack[node].kmer    = kmer;
          stack[node].count   = 0;
          stack[node].mark    = 0;
        }
      stack[node].count += 1;
      stack[node].edge[s8[i]] += 1;
      if (stack[node].count > maxcount)
        { maxcount = stack[node].count;
          maxnode  = node;
        }
    }

  density = (1.*len)/unit;
  printf("\nDB Graph(%d): %d ave = %d/%d = %.1f\n",klen,stop,len,unit,density);
  printf("  max = %d (%d)\n",maxnode,maxcount);
  dcut = .05*density;
  if (dcut < 1)
    dcut = 1;

  { int n, m, p, y, e, x, s;
    int mtop;
    int emax, enod, etod, echr;
    int deep;

    clen = 0;
    mtop = 0;
    stack[maxnode].mark = 2;
    stack[maxnode].depth = 0;
    mark[mtop++] = maxnode;
    for (i = 0; i < 10; i++)
      { emax = 0;
        for (s = 0; s < mtop; s++)
          { n = mark[s];
            p = (stack[n].kmer << 2) & umask;
            for (x = 0; x < 4; x++)
              { if (stack[n].edge[x] <= emax) 
                  continue;
                m = alive[p | x];
                if (stack[m].mark)
                  continue;
                emax = stack[n].edge[x];
                enod = n;
                etod = m;
                echr = x;
              }
          }
        if (emax <= dcut)
          break;

        printf("\n  %d(%d) : %d\n",enod,stack[enod].count,stack[enod].depth);
        deep = stack[enod].depth+1;
        printf("      -%c(%d)-> %d(%d) : %d\n",dna[echr],emax,etod,stack[etod].count,deep);
        n = etod;
        while (stack[n].mark == 0)
          { stack[n].mark = 1;
            stack[n].depth = deep;
            y = 0;
            e = stack[n].edge[y];
            for (x = 1; x < 4; x++)
              if (stack[n].edge[x] > e)
                { y = x;
                  e = stack[n].edge[y];
                }
            p = (stack[n].kmer << 2) & umask;
            n = alive[p | y];
            deep += 1;
            printf("      -%c(%d)-> %d(%d) : %d\n",dna[y],e,n,stack[n].count,deep);
          }
        if (stack[n].depth == deep || (n == maxnode && stack[n].depth == clen))
          printf("   Detour\n");
        if (stack[n].mark == 1)
          printf("   New cycle @ %d\n",n);
        if (clen == 0)
          { clen = deep - stack[n].depth;
            printf("   Main unit len = %d\n",clen);
          }

        printf("\n   %c",dna[echr]);
        n = etod;
        while (stack[n].mark == 1)
          { stack[n].mark = 2;
            mark[mtop++]  = n;
            y = 0;
            e = stack[n].edge[y];
            p = (stack[n].kmer << 2) & umask;
            for (x = 1; x < 4; x++)
              if (stack[n].edge[x] > e)
                { y = x;
                  e = stack[n].edge[y];
                }
            n = alive[p | y];
            printf("%c",dna[y]);
          }
        printf("\n");
      }
  }

  for (i = 0; i < stop; i++)
    { kmer = stack[i].kmer;
      c    = stack[i].count;
      if (c <= dcut)
        continue;
      printf(" %3d: %s(%d)\n",i,emer(kmer,klen),c);
      for (x = 0; x < 4; x++)
        { e = stack[i].edge[x];
          if (e <= dcut)
            continue;
          printf("      -%c(%d)-> %d\n",dna[x],e,alive[((kmer << 2) | x) & umask]);
        }
    }

  for (i = 0; i < stop; i++)
    alive[stack[i].kmer] = 0;

  free(mark);
  free(stack);

  if (clen < .9*unit)
    printf("    Bad call\n");

  return (clen);
}

#endif

static int spectrum_block(uint8 *seq, int off, int len, S_Bundle *bundle)
{ Alignment  *align = bundle->align;
  Overlap    *over  = bundle->over;
  Align_Spec *spec  = bundle->spec;
  Work_Data  *work  = bundle->work;
  OneFile    *ofile = bundle->ofile;
  int64      *t64   = bundle->trace;
  int         tmax  = bundle->tmax;

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
            printf("  p = %d\n",p);
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

            Decompress_TraceTo16(over);

#ifdef SHOW_ALIGNMENTS
#ifndef SHOW_SEARCH
            if (last < 0)
              printf(" %4d: %5d\n",d,diags[d]-diags[d-1]);
            printf("\n");
#endif
            printf(" Hit spans %d-%d\n",bpath->abpos,bpath->bepos);
            Compute_Trace_PTS(Align,Work,100,GREEDIEST,d-wide,d+wide);
            Print_Alignment(stdout,Align,Work,8,100,10,0,10,0);
#endif

            if (bpath->aepos < bpath->bbpos)
              { for (f = (bpath->bbpos-off)+1; f <= end; f++)
                  seq[f] = 4;
                end = bpath->aepos-off;
                for (f = beg+1; f <= end; f++)
                  seq[f] = 4;

#ifdef COMPUTE_MODEL
                printf("  Near Tandem %d-%d\n",bpath->bbpos-bpath->aepos,bpath->aepos-bpath->abpos);
#endif
              }
            else
              { 
#ifdef COMPUTE_MODEL
                compute_model(seq+beg,end-beg,count,d);
#endif

                for (f = beg+1; f <= end; f++)
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
 
  pthread_mutex_lock(&TMUTEX);   //  Put this thread back on the avail stack
    Tstack[Tavail++] = bundle->tid;
  pthread_mutex_unlock(&TMUTEX);

  pthread_cond_signal(&TCOND);   //  Signal a thread is available

  return (NULL);
}

int main(int argc, char *argv[])
{ FILE   **units;
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
    MODELS  = flags['m'];
    MODELS  = 0;              //  Not yet realized

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        fprintf(stderr,"      -m: Compute models of each hit (not yet implemented).\n");
        exit (1);
      }
  }

  //  Get GDB or make a temporary if a fasta

  { char *tpath, *spath;
    char *cpath, *APATH, *AROOT;

    Get_GDB_Paths(argv[1],NULL,&spath,&tpath,0);
  
    units = Get_GDB(gdb,spath,".",NTHREADS);

    free(tpath);

    //  Open 1aln file for threaded writing

    APATH = PathTo(argv[2]);
    AROOT = Root(argv[2],".1aln");
    cpath = getcwd(NULL,0);

    Ofile = open_Aln_Write(Catenate(APATH,"/",AROOT,".1aln"),NTHREADS,Prog_Name,VERSION,
                           Command_Line,TSPACE,spath,NULL,cpath);
    free(cpath);
    free(AROOT);
    free(APATH);
    free(spath);

    Write_Aln_Skeleton(Ofile,gdb);
  }

  if (VERBOSE)
    { fprintf(stderr," Database built, begin scan\n");
      fflush(stderr);
    }

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

printf("Launching %d\n",i); fflush(stdout);

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
      }

    oneFileClose(Ofile);

    Close_GDB(gdb);

    if (VERBOSE)
      { TimeTo(stderr,0,1);
        TimeTo(stderr,1,0);
      }

    Catenate(NULL,NULL,NULL,NULL);
    Numbered_Suffix(NULL,0,NULL);
    free(Prog_Name);

    exit (0);
  }
}
