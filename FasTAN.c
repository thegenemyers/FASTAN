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

#undef   SHOW_BITVEC
#undef   SHOW_PARSE

#define DIAG_MAX        30000
#define DIAG_MIN            2
#define BLOCK_OVERLAP  0x6000

#define TSPACE     100
#define VERSION  "0.5"
#define MBUF_LEN   200  //  Must be even

static char *Usage = "[-vamp] [-T(8)] [-o<target>] <source:path>[<fa_extn>|<1_extn>|.1gdb]";

static int NTHREADS;
static int VERBOSE;
static int MAKE_ALN;
static int MAKE_ANO;
static int PARSE;

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
 *  CUT SEQUENCE INTO UNITS USING BITVECTOR MODEL
 *
 ********************************************************************************************/

#define START(P,M,s)		\
{ P = 0xffffffffffffffffllu;	\
  M = 0;			\
  s = patlen;			\
}

#define START_BLK(p,m)		\
{ p = 0xffffffffffffffffllu;	\
  m = 0;			\
}

//  (P,M) is the delta-state of the current column
//  v is the 0/1 match vector for the symbol to advance over
//  s is the score of the element in the top row before on entry, and after on completion
//  c is whether or not to inject a new start (0), or not (1)

#define ADVANCE(P,M,v,c,s)		\
{ uint64 V, H, Q, N;			\
					\
  V  = v | (c < 0);			\
  H  = (((V & P) + P) ^ P) | V;		\
  V |= M;				\
					\
  Q = M | ~ (H | P);			\
  N = P & H;				\
					\
  H = (Q << 1) | (c > 0);		\
  P = (N << 1) | ~ (V | H) | (c < 0);	\
  M = H & V;				\
					\
  if (Q & bit)				\
    s += 1;				\
  else if (N & bit)			\
    s -= 1;				\
}

//  Here (p,m) is the delta-state, and c is the delta between columns (-1,0,1)

#define ADVANCE_BLK(p,m,v,c) 	 	\
{ uint64 V, H, Q, N, P, M, x;      	\
					\
  P  = p;                       	\
  M  = m;				\
  x  = (c < 0);				\
  V  = v | x;				\
  H  = (((V & P) + P) ^ P) | V;		\
  V |= M;				\
					\
  Q = M | ~ (H | P);			\
  N = P & H;				\
					\
  H = (Q << 1) | (c > 0);		\
  p = (N << 1) | ~ (V | H) | x;		\
  m = H & V;				\
					\
  if (Q & 0x8000000000000000llu)	\
    c = 1;				\
  else if (N & 0x8000000000000000llu)	\
    c = -1;				\
  else					\
    c = 0;				\
}

#define PRINT_BLOCK(p,m) 		\
{ uint64 P, M;			      	\
  int d;				\
					\
  P  = p;                       	\
  M  = m;				\
  for (d = 0; d < 64; d++)		\
    if (P & (((long) 1) << d))		\
      printf("+");			\
    else if (M & (((long) 1) << d))	\
      printf("-");			\
    else				\
      printf("0");			\
  printf(" ");				\
}

#define PRINT_COLUMN(P,M,s,x)	 	\
{ int d;				\
					\
  for (d = 0; d < patrem; d++)		\
    if (P & (((long) 1) << d))		\
      printf("+");			\
    else if (M & (((long) 1) << d))	\
      printf("-");			\
    else				\
      printf("0");			\
  printf(" [%d]",s);			\
  if (x >= 0)				\
    printf(" [%d]",x);			\
  printf("\n");				\
  fflush(stdout);			\
}

static int64 *bitvec_partition(uint8 *seq, int *len, uint8 *pat, int patlen, int unit)
{ int64 *partition;
  int    parlen;

  if (patlen > unit) patlen = unit;

  int nblocks = (patlen>>6);
  int patrem  = (patlen & 0x3f);
  int seqlen  = *len;
  int var     = (unit-2)/3+1;
  uint64 bit  = (0x1llu << (patrem-1));

  parlen    = 0;
  partition = malloc(sizeof(int64)*(1.2*(seqlen/unit)+10));
  if (partition == NULL)
    { fprintf(stderr,"%s: Out of memory partitioning tandem array\n",Prog_Name);
      exit (1);
    }

#ifdef SHOW_BITVEC
  printf("\nBit vect partition: %d = %d blocks + %d rem, unit = %d, var = %d\n",
         patlen,nblocks,patrem,unit,var);
#endif

  { uint64 revt[4][nblocks+1];
    uint64 fort[4][nblocks+1];
    uint64 p[nblocks+1];
    uint64 m[nblocks+1];
    uint64 a[nblocks+1];
    uint64 b[nblocks+1];
  
    //  setup pattern vectors
  
    { uint8 *rat;
      int    i, j, k;
      int    s, n;
      uint64 P, M, Q;
  
      rat = pat+(patlen-1);
      for (j = n = 0; j < patlen; j += 64, n++)
        { for (s = 0; s < 4; s++)
            { P = 0;
              Q = 0;
              M = 1;
              k = j+64;
              for (i = j; i < k; i++)
                { if (i >= patlen)
                    break;
                  if (s == pat[i])
                    P |= M;
                  if (s == rat[-i])
                    Q |= M;
                  M <<= 1;
                }
              for ( ; i < k; i++)
                { P |= M;
                  Q |= M;
                  M <<= 1;
                }
              fort[s][n] = P;
              revt[s][n] = Q;
            }
        }
  
#ifdef SHOW_BITVEC
      { int     d;
        uint64 *v;
  
        for (d = 0; d < 2; d++)
          { printf("     ");
            if (d == 0)
              { for (i = 0; i < patlen; i++)
                  printf("%c",dna[pat[i]]);
              }
            else
              { for (i = 0; i < patlen; i++)
                  printf("%c",dna[rat[-i]]);
              }
            printf("\n");
            for (s = 0; s < 4; s++)
              { if (d == 0)
                  v = fort[s];
                else
                  v = revt[s];
                printf("  %c: ",dna[s]);
                for (i = 0; i < patlen; i++)
                  { if (*v & (0x1llu<<(i&0x3f)))
                      printf("1");
                    else
                      printf("0");
                    if ((i&0x3f) == 63)
                      v += 1;
                  }
                printf("\n");
              }
          }
      }
#endif
    }
   
    //  forward scan
  
    { int    Score[128], s;
      uint64 P, M, *v;
      int    tlow, tmid, thgh;
      int    bloc, ploc, oloc;
      int    best, pest, oest;
      int    i, n, c;
  
#ifdef SHOW_PARSE
      printf("Forward Scan\n");
#endif
  
      for (n = 0; n < nblocks; n++)
        START_BLK(p[n],m[n]);
      START(P,M,s);

      //  find first partition point
  
      thgh = unit + var;
      pest = patlen;
      best = patlen;
      bloc = 0;
      for (i = 0; i < thgh; )
        { c = 0;
          v = fort[seq[i]];
          for (n = 0; n < nblocks; n++)
            ADVANCE_BLK(p[n],m[n],*v++,c);
          ADVANCE(P,M,*v,c,s)
#ifdef SHOW_BITVEC
          printf("  %6d %c: ",i+1,dna[seq[i]]);
          for (n = 0; n < nblocks; n++)
            PRINT_BLOCK(p[n],m[n]);
          PRINT_COLUMN(P,M,s,-1);
#endif
          i += 1;
          if (s < best)
            { best = s; bloc = i; }
          else if (i-bloc >= var && best <= .3*patlen)
            break;
        }
#ifdef SHOW_BITVEC
      printf("Report %d @ %d\n",best,bloc);
      fflush(stdout);
#endif
  
      tmid = bloc + unit;
      tlow = tmid - var;
      thgh = tmid + var;

      oloc = 0;
      oest = -1;
      ploc = bloc;
      pest = best;
      best = patlen;

      //  carry on finding partition points [unit-var,unit+var] from the last one
  
      while (i < seqlen)
        { c = 0;
          v = fort[seq[i]];
          for (n = 0; n < nblocks; n++)
            ADVANCE_BLK(p[n],m[n],*v++,c);
          ADVANCE(P,M,*v,c,s)
#ifdef SHOW_BITVEC
          printf("  %6d %c: ",i+1,dna[seq[i]]);
          for (n = 0; n < nblocks; n++)
            PRINT_BLOCK(p[n],m[n]);
          PRINT_COLUMN(P,M,s,-1);
          fflush(stdout);
#endif
  
          i += 1;
          if (i >= tlow)
            { if (i == thgh)
                {
#ifdef SHOW_BITVEC
                  printf("Report %d @ %d\n",best,bloc);
                  fflush(stdout);
#endif
                  if (unit == patlen && pest != 0)
                    { uint64 A, B;
                      int    nloc, nest;
                      int    j, t, x, y;
  
                      tlow = ploc-var;
                      if (tlow <= oloc)
                        tlow = oloc+1;
                      thgh = ploc+var;
                      pest = (patlen<<3);

                      for (n = 0; n < nblocks; n++)
                        START_BLK(a[n],b[n]);
                      START(A,B,t);
                      for (j = oloc; j <= thgh; )
                        { c = 1;
                          v = fort[seq[j]];
                          for (n = 0; n < nblocks; n++)
                            ADVANCE_BLK(a[n],b[n],*v++,c);
                          ADVANCE(A,B,*v,c,t);
                          j += 1;
                          Score[j-oloc] = t;
                        }
                          
                      for (n = 0; n < nblocks; n++)
                        START_BLK(a[n],b[n]);
                      START(A,B,t);
                      for (j = bloc-1; j >= tlow; j--)
                        { c = 1;
                          v = revt[seq[j]];
                          for (n = 0; n < nblocks; n++)
                            ADVANCE_BLK(a[n],b[n],*v++,c);
                          ADVANCE(A,B,*v,c,t);
                          if (j > ploc+var)
                            continue;
                          y = Score[j-oloc];
                          x = t+y;
#ifdef SHOW_BITVEC
                          printf("  %6d %c: ",j,dna[seq[j]]);
                          for (n = 0; n < nblocks; n++)
                            PRINT_BLOCK(a[n],b[n]);
                          PRINT_COLUMN(A,B,t,y);
#endif
                          if (x < pest)
                            { pest = x;
                              nloc = j;
                              nest = t;
                            }
                          else if (x == pest && t < nest) 
                            { nloc = j;
                              nest = t;
                            }
                        }
#ifdef SHOW_BITVEC
                      printf("Best is %d @%d delta = %d\n",pest,nloc,nloc-ploc);
                      fflush(stdout);
#endif
                      oest = pest-nest;
                      ploc = nloc;
                      pest = nest;
                    }
#ifdef SHOW_PARSE
                  if (oloc == 0)
                    printf("       ");
                  else
                    printf("%5d: ",parlen-1);
                  Print_Seq(seq+oloc,ploc-oloc);
                  if (oloc == 0)
                    printf(" ? [%d..%d] tail\n",oloc,ploc);
                  else
                    printf(" %d [%d..%d]\n",oest,oloc,ploc);
                  fflush(stdout);
#else
                  (void) oest;
#endif
                  if (oloc < ploc)
                    partition[parlen++] = ploc;

                  tmid = bloc + unit;
                  tlow = tmid - var;
                  thgh = tmid + var;

                  oest = pest;
                  oloc = ploc;
                  pest = best;
                  ploc = bloc;
                  best = patlen;
                }
              else
                { if (s < best)
                    { best = s; bloc = i; }
                  else if (s == best && abs(bloc-tmid) >= abs(i-tmid))
                    { best = s; bloc = i; } 
                }
            }
        }
#ifdef SHOW_PARSE
      printf("%5d: ",parlen-1);
      Print_Seq(seq+oloc,ploc-oloc);
      printf(" %d [%d..%d] \n",oest,oloc,ploc);
#endif
      partition[parlen++] = ploc;
      if (best < .3*patlen && bloc < seqlen-10)    //  the last minimum probably a partition point
        { oloc = ploc;                             //  if <30% difference and a bit away from end
          ploc = bloc;                             //   of sequence
#ifdef SHOW_PARSE
          printf("%5d: ",parlen-1);
          Print_Seq(seq+oloc,ploc-oloc);
          printf(" %d [%d..%d] \n",pest,oloc,ploc);
#endif
          partition[parlen++] = bloc;              //  not harmful to include it in consensus
                                                   //   evaluation even if truncated.
        }
#ifdef SHOW_PARSE
      printf("       ");
      Print_Seq(seq+ploc,seqlen-ploc);
      printf(" ? [%d..%d] tail\n",ploc,seqlen);
      fflush(stdout);
#endif
    }
  }

  *len = parlen;
  return (partition);
}

  //  return 1st position of x-mer with the highest sum of the x-7 of 8-mer counts within it

static int best_seed(uint8 *seq, int len, uint16 *alive, int xmer, int unit)
{ int     i, x, y, r;
  uint16  kmer;
  int     best, kloc;
  int     count[xmer], wrap, score;
  int     wmax;

  wrap = xmer-7;
  wmax = len/unit+1;

  kmer = seq[0];
  for (i = 1; i < 7; i++)
    kmer = (kmer << 2) | seq[i];
  for (i = 7; i < len; i++)
    { kmer = ((kmer << 2) | seq[i]) & 0xffff;
      x = alive[kmer];
      if (x < wmax)
        alive[kmer] = x+1;
    }

  kmer = seq[0];
  for (i = 1; i < 7; i++)
    kmer = (kmer << 2) | seq[i];
  r = 0;
  score = 0;
  for (i = 7; i < xmer; i++)
    { kmer = ((kmer << 2) | seq[i]) & 0xffff;
      x = alive[kmer];
      count[r++] = x;
      score += x;
    }
  r = 0;
  best = score;
  kloc = xmer-1;
  for (i = xmer; i < len; i++)
    { kmer = ((kmer << 2) | seq[i]) & 0xffff;
      y = count[r];
      x = alive[kmer];
      count[r] = x;
      score += x-y;
      if (score > best)
        { best = score;
          kloc = i; 
        }
      if (++r == wrap)
        r = 0;
    }
  kloc -= xmer-1;

  if (len < 0x8000)
    { for (i = 0; i < 0x10000; i++)
        alive[i] = 0;
    }
  else
    { kmer = seq[0];
      for (i = 1; i < 7; i++)
        kmer = (kmer << 2) | seq[i];
      for (i = 7; i < len; i++)
        { kmer = ((kmer << 2) | seq[i]) & 0xffff;
          alive[kmer] = 0;
        }
    }

#ifdef SHOW_BITVEC
  printf("Best %d-mer is at %d with 8-mer score %d (",xmer,kloc,best);
  Print_Seq(seq+kloc,xmer);
  printf(")\n");
#endif

  return (kloc);
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

  OneFile    *mfile   = bundle->mfile;
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
    p = diags[DIAG_MIN-1];
    for (i = DIAG_MIN; i < DIAG_MAX; i++)
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
            printf("  p = %d (%d %d) %d\n",p,last,hits[x-1].ibeg,off);
            fflush(stdout);
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
            fflush(stdout);
#endif

            end = path->aepos - off;
            beg = path->bbpos - off;

            if (beg > p || end <= p)
              continue;
            if (end > last)
              last = end;
            if (end-beg < 1.95*d)
              continue;

            unit = ave_tp_diag(path);

            if (over->path.tlen > tmax)
              { tmax = bundle->tmax = 1.2*over->path.tlen + 1000;
                t64  = bundle->trace = realloc(t64,sizeof(int64)*tmax);
              }

            if (MAKE_ALN)
              { Write_Aln_Overlap(ofile,over);
                Compress_TraceTo8(over,0);
                Write_Aln_Trace(ofile,over->path.trace,over->path.tlen,t64,unit);
              }

            if (MAKE_ANO)
              { char mstring[20];
                int  mlen;

                oneInt(mfile,0) = mscf;
                oneInt(mfile,1) = path->bbpos + moff;
                oneInt(mfile,2) = path->aepos + moff;
                oneWriteLine(mfile,'M',0,NULL);
                mlen = sprintf(mstring,"%d",unit);
                oneWriteLine(mfile,'L',mlen,mstring);
                oneInt(mfile,0) = (100.*path->diffs) / (path->aepos-path->abpos);
                oneWriteLine(mfile,'X',0,NULL);
              }

#if defined(SHOW_ALIGNMENTS) || defined(DO_CUT)
            if (MAKE_ALN)
              Decompress_TraceTo16(over);
            Compute_Trace_PTS(align,work,100,GREEDIEST,d-wide,d+wide);
#endif

#ifdef SHOW_ALIGNMENTS
#ifndef SHOW_SEARCH
            if (last < 0)
              printf(" %4d: %5d\n",d,diags[d]-diags[d-1]);
            printf("\n");
            fflush(stdout);
#endif
            printf(" Hit spans %d-%d (unit = %d)\n",path->bbpos,path->aepos,unit);
            Print_Reference(stdout,align,work,8,100,10,0,10,0);
            fflush(stdout);
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
                fflush(stdout);
#endif
              }
            else
              { if (PARSE)
                  { int    xmer, kloc, len;
                    int64 *part;

                    if (unit < 8)
                      xmer = 8;
                    else if (unit <= 64)
                      xmer = unit;
                    else
                      xmer = (unit-64)*.07 + 64;  // sync seed length is 7% of unit length

                    kloc = best_seed(seq+beg,end-beg,count,xmer,unit);
                    len  = end-beg;
                    part = bitvec_partition(seq+beg,&len,seq+beg+kloc,xmer,unit);

                    oneWriteLine(mfile,'P',len,part);

                    free(part);
                  }

                for (f = beg; f < end; f++)
                  seq[f] = 4;
              }

            if (path->aepos > outhit)
              outhit = path->aepos;
          }
      }

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

  bundle->mscf = gdb->contigs[i].scaf;
  bundle->moff = gdb->contigs[i].sbeg;

  last = -1;
  if (clen <= 0x8000)
    spectrum_block(buffer,0,clen,bundle);
  else
    for (p = 0; p+7 < clen; p += BLOCK_OVERLAP)
      { if (p+0x8000 >= clen)
          { spectrum_block(buffer+p,p,clen-p,bundle);
            break;
          }
        else
          { last = spectrum_block(buffer+p,p,0x8000,bundle);
            if (last >= p+BLOCK_OVERLAP)
              p = last-BLOCK_OVERLAP;
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
  char      *TRGT_PATH;
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

    NTHREADS  = 8;
    TRGT_PATH = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vamp")
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

    VERBOSE  = flags['v'];
    MAKE_ALN = flags['a'];
    MAKE_ANO = flags['m'];
    PARSE    = flags['p'];

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -a: Make a .1aln of all the hits found.\n");
        fprintf(stderr,"      -m: Make a .1ano mask of the hits found.\n");
        fprintf(stderr,"      -p: Parse hits into units.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -o: Root path of .1aln/.1ano file (default root path of input).\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        exit (1);
      }

    if (! (MAKE_ALN || MAKE_ANO))
      { fprintf(stderr,"%s: Must set at least one of -a or -m options\n",Prog_Name);
        exit (1);
      }
    if (PARSE && !MAKE_ANO)
      { fprintf(stderr,"%s: Must set at least one of -m if set -p option\n",Prog_Name);
        exit (1);
      }
  }

  //  Get GDB or make a temporary if a fasta

  { char *cpath, *APATH, *AROOT;
    int   ftype;

    ftype = Get_GDB_Paths(argv[1],NULL,&spath,&cpath,0);

    free(cpath);

    if (MAKE_ANO && ftype != IS_GDB)
      { fprintf(stderr,"%s: The source must be a GDB when masking (-m) is on\n",Prog_Name);
        exit (1);
      }
  
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

    //  Open output files

    if (MAKE_ALN)
      { Ofile = open_Aln_Write(Catenate(APATH,"/",AROOT,".1aln"),NTHREADS,Prog_Name,VERSION,
                               Command_Line,TSPACE,spath,NULL,cpath);

        Write_Skeleton(Ofile,gdb);
      }

    if (MAKE_ANO)
      { anoSchema = make_ANO_Schema();
        Mfile = oneFileOpenWriteNew(Catenate(APATH,"/",AROOT,".1ano"),anoSchema,"ano",1,NTHREADS);

        oneAddProvenance(Mfile,Prog_Name,VERSION,Command_Line);

        oneAddReference(Mfile,gdb->srcpath,1);

        Write_Skeleton(Mfile,gdb);
      }

    free(cpath);
    free(AROOT);
    free(APATH);
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
        if (MAKE_ALN)
          parm[i].ofile = Ofile + i;
        else
          parm[i].ofile = NULL;

        parm[i].gdb  = &parm[i]._gdb;
        parm[i]._gdb = _gdb;
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
        if (MAKE_ANO)
          parm[i].mfile = Mfile + i;
        else
          parm[i].mfile = NULL;
        parm[i].tmax    = 10000;
        parm[i].trace   = malloc(sizeof(int64)*10000);

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
      { free(parm[i].trace);
        free(parm[i].buffer-1);
        free(parm[i].block);
        if (i == 0)
          Free_Align_Spec(parm[i].spec);
        else
          fclose(units[i]);
        Free_Work_Data(parm[i].work);
      }

    if (MAKE_ANO)
      { oneFileClose(Mfile);
        oneSchemaDestroy(anoSchema);
      }
    if (MAKE_ALN)
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
