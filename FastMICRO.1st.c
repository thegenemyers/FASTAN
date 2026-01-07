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

#define  ID_MIN  .7
#define  ID_REG  12  //  ceiling(ID_MIN * 16)

#define MIN_SPAN 8

#define TSPACE     100
#define VERSION  "0.5"
#define DIAG_MAX 30000
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

static int GoodX[0x10000];
static int GoodR[0x10000];

static void suffix_positive(int val, int total, int good)
{ if (total == 16)
    GoodX[val] = 1;
  else
    { suffix_positive((0x1<<total)|val,total+1,good+1);
      if (good > ID_MIN*(total+1))
        suffix_positive(val,total+1,good);
    }
}

static void good_region(int val, int bits, int good)
{ if (good < ID_REG)
    { if ((16-bits)+good >= ID_REG)
        { val <<= 1;
          bits += 1;
          good_region(val,bits,good);
          good_region(val+1,bits,good+1);
        }
    }
  else
    { int i;

      val <<= (16-bits);
      for (i = (1<<(16-bits)) - 1; i >= 0; i--)
        GoodR[val+i] = 1;
    }
}

static void Compute_Terminators()
{ int i;

  for (i = 0; i < 0x10000; i++)
    GoodX[i] = GoodR[i] = 0;

  suffix_positive(0,0,0);

  good_region(0,0,0);

  // for (i = 0xffff; i >= 0; i--)
    // printf(" %04x: %d %d\n",i,GoodX[i],GoodR[i]);
}

static int extend_right(uint8 * seq, int slen, uint8 *unit, int ulen)
{ int k, d;
  int y, z;
  int b, c;
  int j;

  int    _F[1001], *F = _F+500;
  uint16 _B[1001], *B = _B+500;

  j = 0;
  for (z = 0; z < slen; z++)
    { if (seq[z] != unit[j])
        break;
      z += 1;
      j += 1;
      if (j >= ulen)
        j = 0;
    }
  F[0] = z;
  B[0] = 0xffff;
  printf("Wave  0: %4d\n         %4x\n         %4d\n",F[0],B[0],GoodX[B[0]]);

  for (d = 1; d < 10; d++)
    { k = -d;
      y = F[k-1] = F[d+1] = F[k] = F[d] = -1;
      c = 0;
      for ( ; k <= d; k++)
        { z = y;
          y = F[k];
          b = c;
          c = B[k];
          if (z < y)
            { z = y;
              b = c;
            }
          if (z < F[k+1])
            { z = F[k+1];
              b = B[k+1];
            }
          else
            z += 1;
          b <<= 1;
          j = (z-k) % ulen;
          while (z < slen)
            { if (seq[z] != unit[j])
                break;
              z += 1;
              j += 1;
              b = (b<<1) | 1;
              if (j >= ulen)
                j = 0;
            }
          F[k] = z;
          B[k] = b;
        }
      printf("Wave %2d:",d);
      for (k = -d; k <= d; k++)
        printf(" %4d",F[k]);
      printf("\n        ");
      for (k = -d; k <= d; k++)
        printf(" %04x",B[k]);
      printf("\n        ");
      for (k = -d; k <= d; k++)
        printf(" %4d",GoodX[B[k]]);
      printf("\n");
      fflush(stdout);
    }
  return (0);
}

/*
static int extend_right2(uint8 *seq, int len, int a, int b)
{ int valab, valba;
  int x, y, z;
  int i, val, score;
  int ans;

  valab = (a<<2) | b;
  valba = (b<<2) | a;

  score = 0;

  ans = 0;
  y = seq[-2];
  z = seq[-1];
  for (i = 0; i < len; i++)
    { x = y;
      y = z;
      z = seq[i];
      val = (z<<2)|y;
      if (z == a || z == b)
        { if (z == x && z != y)
            score += 1;
          if (val == valab || val == valba)
            { if (score >= 0)
                ans = i;
              printf(" @%d: %d\n",i,score);
            }
        }
      else
        score -= 2;
      if (score < -10)
        break;
    }
  return (ans);
}

static int extend_right3(uint8 *seq, int len, int a, int b, int c)
{ int valabc, valbca, valcab;
  int u, x, y, z;
  int i, val, score;
  int ans;

  valabc = (a<<4) | (b<<2) | c;
  valbca = (b<<4) | (c<<2) | a;
  valcab = (c<<4) | (a<<2) | b;

  score = 0;

  ans = 0;
  x = seq[-3];
  y = seq[-2];
  z = seq[-1];
  for (i = 0; i < len; i++)
    { u = x;
      x = y;
      y = z;
      z = seq[i];
      val = (z<<4)|(y<<2)|x;
      if (z == a || z == b || z == c)
        { if (z == u && (z != y || z != x))
            score += 1;
          if (val == valabc || val == valbca || val == valcab)
            { if (score >= 0)
                ans = i;
              printf(" @%d: %d\n",i,score);
            }
        }
      else
        score -= 2;
      if (score < -10)
        break;
    }
  return (ans);
}
*/

static void micro_contig(uint8 *seq, int len, S_Bundle *bundle)
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

  int lastr, lastl;
  int a, b, c;
  int i, j;
  uint8 unit[4];

/*
  lastr = -1;
  b = seq[0];
  for (i = 0; i < len-MIN_SPAN; i++)
    { a = b;
      b = seq[i+1];
      if (a != b)
        { for (j = i+2; j < len; j += 2)
            if (a != seq[j] || b != seq[j+1])
              break;
          if (j-i >= MIN_SPAN)
            { printf("Hit @%d: ",i);
              Print_Seq(seq+i,j-i);
              printf("\n");
              if (lastr >= 0)
                { printf("Extend from %d: ",lastr);
                  Print_Seq(seq+lastl,lastr-lastl);
                  printf(".");
                  Print_Seq(seq+lastr,j-lastr);
                  printf("\n");
                  lastr += extend_right(seq+lastr,j-lastr,unit,2);
                }
              if (lastr > i+1)
                lastr = j; 
              else
                { printf("\nReport %d-%d %c%c:   \n",lastl,lastr,dna[unit[0]],dna[unit[1]]);
                  Print_Seq(seq+lastl,lastr-lastl);
                  printf("\n\n");
                  unit[0] = a;
                  unit[1] = b;
                  // i -= extend_left(seq+lastr,i-lastr,unit,2);
                  lastl = i;
                  lastr = j;
                }
	      i = j;
            }
        }
    }
  if (lastr >= 0)
    extend_right(seq+lastr,len-lastr,unit,2);
*/


  lastr = -1;
  b = seq[0];
  c = seq[1];
  for (i = 0; i < len-MIN_SPAN; i++)
    { a = b;
      b = c;
      c = seq[i+2];
      if (a != b || a != c)
        { for (j = i+3; j < len; j += 3)
            if (a != seq[j] || b != seq[j+1] || c != seq[j+2])
              break;
          if (j-i >= MIN_SPAN)
            { printf("Hit @%d: ",i);
              Print_Seq(seq+i,j-i);
              printf("\n");
              if (lastr >= 0)
                { printf("Extend from %d: ",lastr);
                  Print_Seq(seq+lastl,lastr-lastl);
                  printf(".");
                  Print_Seq(seq+lastr,j-lastr);
                  printf("\n");
                  lastr += extend_right(seq+lastr,j-lastr,unit,3);
                }
              if (lastr > i+1)
                lastr = j; 
              else
                { if (lastr >= 0)
                    { printf("\nReport %d-%d %c%c%c:   \n",lastl,lastr,dna[unit[0]],dna[unit[1]],dna[unit[2]]);
                      Print_Seq(seq+lastl,lastr-lastl);
                      printf("\n\n");
                    }
                  unit[0] = a;
                  unit[1] = b;
                  unit[2] = c;
                  // i -= extend_left3(seq+lastr,i-lastr,unit,3);
                  lastl = i;
                  lastr = j;
                }
	      i = j;
            }
        }
    }
  if (lastr >= 0)
    extend_right(seq+lastr,len-lastr,unit,3);

  // Write_Aln_Overlap(ofile,over);
  // Compress_TraceTo8(over,0);
  // Write_Aln_Trace(ofile,over->path.trace,over->path.tlen,t64,unit);

/*
  if (MAKE_MASK)
    { if (mint >= MBUF_LEN)
        { oneInt(mfile,0) = mscf;
          oneWriteLine(mfile,'M',mint,mlist);
          oneWriteLine(mfile,'L',mstr-mstring,mstring);
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
*/

  bundle->mint = mint;
  bundle->mstr = mstr;
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
  int       i, clen;
 
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

  micro_contig(buffer,clen,bundle);

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

    Compute_Terminators();

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

        break;
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
