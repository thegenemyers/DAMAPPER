/*******************************************************************************************
 *
 *  Fast local alignment filter for long, noisy reads based on "dumbing down" of my RECOMB 2005
 *     filter with Jens Stoye, and a "smarting up" of the k-mer matching by turning it into
 *     a threaded sort and merge paradigm using a super cache coherent radix sort.  Local
 *     alignment is accomplised with dynamically-banded O(nd) algorithm that terminates when
 *     it fails to find a e-matching patch for a significant distance, and polishes the match
 *     to the last e-prefix-positive 32-mer.
 *
 *  Author :  Gene Myers
 *  First  :  June 2013
 *  Current:  June 1, 2014
 *
 ********************************************************************************************/

//  A complete threaded code for the filter

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "DB.h"
#include "align.h"
#include "map.h"

#define THREAD    pthread_t

#define MAX_BIAS  2    //  In -b mode, don't consider tuples with specificity
                       //     <= 4 ^ -(kmer-MAX_BIAS)

#define MAXGRAM 10000  //  Cap on k-mer count histogram (in count_thread, merge_thread)

#define BINWIDTH  64   //  minimum width of a diagonal bucket
#define HITMIN    3    //  minimum # of k-mers in a potential seed

#define MAX_GAP   1000 //  maximum separation between two k-mers in a chain
#define MIN_PIECE  300 //  two seeds are distinct iff one covers uniquely >= MIN_PIECE bp's

#define MATCH_CHUNK    100     //  Max expected number of hits between two reads
#define TRACE_CHUNK  20000     //  Max expected trace points in hits between two reads

#define CHAIN_OFF      500.    //  Two hits with Agap and Bgap bases between them in A and B
#define CHAIN_OVL      400.    //     are chainable if they overlap by less than CHAIN_OVL in
#define CHAIN_PLAY     1.4     //     both A and B, and
                               //  1/CHAIN_PLAY <= (Agap+CHAIN_OFF)/(Bgap+CHAIN_OFF) <= CHAIN_PLAY

#define DIFF_SCORE     2.3     //  Score of a LA = length - DIFF_SCORE*diffs
#define TIE_SCORE       50     //  Score is tied if less than this different
#define TIE_GAP        500     //  Two chain edges are tied if the gap difference is less than this

#undef  TEST_LSORT
#undef  TEST_KSORT
#undef  TEST_PAIRS
#undef  TEST_CSORT
#define    HOW_MANY   3000   //  Print first HOW_MANY items for each of the TEST options above

#undef  TEST_CANDID
#undef  TEST_CHAIN
#undef  TEST_CONTAIN
#undef  TEST_ALIGN
#undef  TEST_SELECT

#if defined(TEST_CHAIN) || defined(TEST_CANDID)
#define NOTHREAD
#endif

#if defined(TEST_CONTAIN) || defined(TEST_ALIGN) || defined(TEST_SELECT)
#define NOTHREAD
#endif

typedef struct
  { uint64 p1;   //  The lower half
    uint64 p2;
  } Double;

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__

typedef struct
  { uint64 code;
    int    rpos;
    int    read;
  } KmerPos;

typedef struct
  { int    diag;
    int    apos;
    int    bread;
    int    aread;
  } SeedPair;

#else

typedef struct
  { uint64 code;
    int    read;
    int    rpos;
  } KmerPos;

typedef struct
  { int    apos;
    int    diag;
    int    aread;
    int    bread;
  } SeedPair;

#endif

/*******************************************************************************************
 *
 *  PARAMETER SETUP
 *
 ********************************************************************************************/

static int Kmer;
static int Suppress;

static int    Kshift;         //  2*Kmer
static uint64 Kmask;          //  4^Kmer-1
static int    TooFrequent;    //  (Suppress != 0) ? Suppress : INT32_MAX

static int    NTHREADS;       //  Adjusted downward to nearest power of 2
static int    NSHIFT;         //  NTHREADS = 1 << NSHIFT

int Set_Filter_Params(int kmer, int suppress, int nthreads)
{ if (kmer <= 1)
    return (1);

  Kmer     = kmer;
  Suppress = suppress;

  Kshift = 2*Kmer;
  if (Kmer == 32)
    Kmask = 0xffffffffffffffffllu;
  else
    Kmask = (0x1llu << Kshift) - 1;

  if (Suppress == 0)
    TooFrequent = INT32_MAX;
  else
    TooFrequent = Suppress;

  NTHREADS = 1;
  NSHIFT   = 0;
  while (2*NTHREADS <= nthreads)
    { NTHREADS *= 2;
      NSHIFT   += 1;
    }

  return (0);
}


/*******************************************************************************************
 *
 *  LEXICOGRAPHIC SORT
 *
 ********************************************************************************************/

#define BMER      4
#define BSHIFT    8             //  = 2*BMER
#define BPOWR   256             //  = 2^BSHIFT
#define BMASK  0xffllu          //  = BPOWR-1

static uint32  NMASK;           //  = NTHREADS-1
static uint64  QMASK;           //  = BMASK << NSHIFT
static int     LEX_shift;
static int64  *LEX_zthresh;
static int     LEX_last;
static int     LEX_next;
static Double *LEX_src;
static Double *LEX_trg;

typedef struct
  { int64  beg;
    int64  end;
    int64  tptr[BPOWR];
    int64 *sptr;
    uint64 table[BPOWR];
  } Lex_Arg;

static void *lex_thread(void *arg)
{ Lex_Arg    *data  = (Lex_Arg *) arg;
  int64      *sptr  = data->sptr;
  int64      *tptr  = data->tptr;
  uint64     *table = data->table;
  int         shift = LEX_shift;   //  Must be a multiple of 8 in [0,120]
  int        qshift = (LEX_next - LEX_shift) - NSHIFT;
  Double     *src   = LEX_src;
  int64      *thrsh = LEX_zthresh;
  Double     *trg   = LEX_trg;
  int64       i, n, x;
  uint64      c, b;
  int         d, q;

  n = data->end;
  if (shift >= 64)
    { shift -= 64;
      if (LEX_last)
        for (i = data->beg; i < n; i++)
          { c = src[i].p2;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
          }
      else
        for (i = data->beg; i < n; i++)
          { c = src[i].p2;
            b = (c >> shift);
            d = (b & BMASK);
            x = tptr[d]++;
            trg[x] = src[i];
            q = table[d];
            if (q >= NTHREADS)
              { q &= NMASK;
                if (x >= thrsh[q])
                  { table[d] -= NMASK;
                    q += 1;
                  }
              }
            sptr[((b >> qshift) & QMASK) | q] += 1;
          }
    }

  else if ( ! LEX_last && LEX_next >= 64)   //  && LEX_shift < 64

    { qshift = (LEX_next - 64) - NSHIFT;
      if (qshift < 0)
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            b = (c >> shift);
            d = (b & BMASK);
            x = tptr[d]++;
            trg[x] = src[i];
            q = table[d];
            if (q >= NTHREADS)
              { q &= NMASK;
                if (x >= thrsh[q])
                  { table[d] -= NMASK;
                    q += 1;
                  }
              }
            sptr[((src[i].p2 << NSHIFT) & QMASK) | q] += 1;
          }
      else
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            b = (c >> shift);
            d = (b & BMASK);
            x = tptr[d]++;
            trg[x] = src[i];
            q = table[d];
            if (q >= NTHREADS)
              { q &= NMASK;
                if (x >= thrsh[q])
                  { table[d] -= NMASK;
                    q += 1;
                  }
              }
            sptr[((src[i].p2 >> qshift) & QMASK) | q] += 1;
          }
    }

  else // LEX_last || LEX_next < 64
    if (LEX_last)
      if (shift == 0)
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            x = tptr[c&BMASK]++;
            trg[x] = src[i];
          }
      else
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
          }
    else
      if (shift == 0)
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            d = (c & BMASK);
            x = tptr[d]++;
            trg[x] = src[i];
            q = table[d];
            if (q >= NTHREADS)
              { q &= NMASK;
                if (x >= thrsh[q])
                  { table[d] -= NMASK;
                    q += 1;
                  }
              }
            sptr[((c >> qshift) & QMASK) | q] += 1;
          }
      else
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            b = (c >> shift);
            d = (b & BMASK);
            x = tptr[d]++;
	    trg[x] = src[i];
            q = table[d];
            if (q >= NTHREADS)
              { q &= NMASK;
                if (x >= thrsh[q])
                  { table[d] -= NMASK;
                    q += 1;
                  }
              }
	    sptr[((b >> qshift) & QMASK) | q] += 1;
          }

  return (NULL);
}

static Double *lex_sort(int bytes[16], Double *src, Double *trg, Lex_Arg *parmx)
{ THREAD  threads[NTHREADS];
  int64   zthresh[NTHREADS];

  int64   len, x, y;
  Double *xch;
  int     i, j, k, z;
  int     b, c, fb;

  len       = parmx[NTHREADS-1].end;
  LEX_src   = src;
  LEX_trg   = trg;
  QMASK     = (BMASK << NSHIFT);
  NMASK     = NTHREADS-1;

  z = (len-1)/NTHREADS + 1;
  for (i = 0; i < NTHREADS; i++)
    zthresh[i] = (i+1)*z;
  LEX_zthresh = zthresh;

  for (c = 0; c < 16; c++)
    if (bytes[c])
      break;
  fb = c;
  for (b = c; b < 16; b = c)
    { for (c = b+1; c < 16; c++)
        if (bytes[c])
          break;
      LEX_last  = (c >= 16);
      LEX_shift = (b << 3);
      LEX_next  = (c << 3);

      if (b == fb)
        { for (i = 0; i < NTHREADS; i++)
            for (z = 0; z < NTHREADS*BPOWR; z++)
              parmx[i].sptr[z] = 0;
        }
      else
        { x = 0;
          for (i = 0; i < NTHREADS; i++)
            { parmx[i].beg = x;
              x = LEX_zthresh[i];
              if (x > len)
                x = len;
              parmx[i].end = x;
              for (j = 0; j < BPOWR; j++)
                parmx[i].tptr[j] = 0;
            }
          parmx[NTHREADS-1].end = len;

          for (j = 0; j < BPOWR; j++)
            { k = (j << NSHIFT);
              for (z = 0; z < NTHREADS; z++)
                for (i = 0; i < NTHREADS; i++)
                  { parmx[i].tptr[j] += parmx[z].sptr[k+i];
                    parmx[z].sptr[k+i] = 0;
                  }
            }
	}

      k = 0;
      x = 0;
      for (j = 0; j < BPOWR; j++)
        for (i = 0; i < NTHREADS; i++)
          { y = parmx[i].tptr[j];
            parmx[i].tptr[j] = x;
            if (x >= LEX_zthresh[k])
              k += 1;
            parmx[i].table[j] = k;
            while (x+y > LEX_zthresh[k])
              { parmx[i].table[j] += NTHREADS;
                k += 1;
              }
            x += y;
          }

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,lex_thread,parmx+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);

      xch     = LEX_src;
      LEX_src = LEX_trg;
      LEX_trg = xch;

#ifdef TEST_LSORT
      printf("\nLSORT %d\n",LEX_shift);
      if (LEX_shift >= 64)
        { x = (1llu << ((LEX_shift-64)+BSHIFT))-1;
          for (i = 0; i < len; i++)
            if (i > 0 && (LEX_src[i].p1 < LEX_src[i].p1 ||
                               (LEX_src[i].p1 == LEX_src[i].p1 && 
                               (LEX_src[i].p2 & x) < (LEX_src[i-1].p2 & x))))
              { printf("%6d: %8llx %8llx %8llx %8llx : %4llx",
                       i,LEX_src[i].p2>>32,(LEX_src[i].p2)&0xffffffffll,LEX_src[i].p1>>32,
                       LEX_src[i].p1&0xffffffffll,LEX_src[i].p2&x);
                printf(" OO");
                printf("\n");
              }
            else
              { printf("%6d: %8llx %8llx %8llx %8llx : %4llx",
                       i,LEX_src[i].p2>>32,(LEX_src[i].p2)&0xffffffffll,LEX_src[i].p1>>32,
                       LEX_src[i].p1&0xffffffffll,LEX_src[i].p2&x);
                printf("\n");
              }
        }
      else
        { x = (1llu << (LEX_shift+BSHIFT))-1;
          for (i = 0; i < len; i++)
            if (i > 0 && ((LEX_src[i].p1 & x) >> LEX_shift) < ((LEX_src[i-1].p1 & x) >> LEX_shift))
              { printf("%6d: %8llx %8llx %8llx %8llx : %4llx",
                       i,LEX_src[i].p2>>32,(LEX_src[i].p2)&0xffffffffll,LEX_src[i].p1>>32,
                       LEX_src[i].p1&0xffffffffll,(LEX_src[i].p1&x)>>LEX_shift);
                printf(" OO");
                printf("\n");
              }
            else
              { printf("%6d: %8llx %8llx %8llx %8llx : %4llx",
                       i,LEX_src[i].p2>>32,(LEX_src[i].p2)&0xffffffffll,LEX_src[i].p1>>32,
                       LEX_src[i].p1&0xffffffffll,(LEX_src[i].p1&x)>>LEX_shift);
                printf("\n");
              }
        }
#endif
    }

  return (LEX_src);
}


/*******************************************************************************************
 *
 *  INDEX BUILD
 *
 ********************************************************************************************/

static int *NormShift = NULL;
static int  LogNorm, LogThresh;
static int  LogBase[4];

static DAZZ_DB    *TA_block;
static KmerPos    *TA_list;
static DAZZ_TRACK *TA_track;

typedef struct
  { int64 *kptr;
    int    fill;
    int    tbeg;
    int    tend;
    int    rbeg;
    int    rend;
  } Tuple_Arg;

static void *tuple_thread(void *arg)
{ Tuple_Arg  *data  = (Tuple_Arg *) arg;
  int64      *kptr  = data->kptr;
  DAZZ_READ  *reads = TA_block->reads;
  KmerPos    *list  = TA_list;
  int         i, j;
  int         n, m;
  char       *s;

  n = data->tbeg;
  i = data->rbeg;
  m = data->tend;
  j = data->rend;
  s = ((char *) TA_block->bases) + reads[i].boff;

  if (TA_track != NULL)
    { int64     *anno1 = ((int64 *) (TA_track->anno)) + 1;
      int       *point = (int *) (TA_track->data);
      int        p, q;
      int        io, jo;
      int64      a, b, f; 

      io = n + i*Kmer - reads[i].boff;
      jo = TA_block->maxlen;
      q  = 0;
      f  = anno1[i-1];
      for (; i <= j; i++)
        { if (i == j)
            { jo = m + j*Kmer - reads[j].boff;
              if (jo == 0)
                continue;
              jo += Kmer-1;
            }
          b = f;
          f = anno1[i];
          for (a = b; a <= f; a += 2)
            { if (a == b)
                p = 0;
              else
                p = point[a-1];
              if (a == f)
                q = reads[i].rlen;
              else
                q = point[a];
              if (p < io)
                p = io;
              if (q > jo)
                q = jo;
              if (p+Kmer <= q)
                { int    x;
                  uint64 c;

                  c = 0;
                  for (x = 1; x < Kmer; x++)
                    c = (c << 2) | s[p++];
                  while (p < q)
                    { x = s[p];
                      c = ((c << 2) | x) & Kmask;
                      list[n].read = i;
                      list[n].rpos = p++;
                      list[n].code = c;
                      n += 1;
                      kptr[c & BMASK] += 1;
                    }
                }
            }
          s += (q+1);
          io = 0;
        }

      kptr[BMASK] += (data->fill = m-n);
      while (n < m)
        { list[n].code = 0xffffffffffffffffllu;
          list[n].read = 0xffffffff;
          list[n].rpos = 0xffffffff;
          n += 1;
        }
    }

  else
    { int    p, q;
      int    io;

      io = n + i*Kmer - reads[i].boff;
      for (; i <= j; i++)
        { p = io;
          if (i == j)
            q = m + j*Kmer + (Kmer-1) - reads[j].boff;
          else
            q = reads[i].rlen;
          if (p+Kmer <= q)
            { int    x;
              uint64 c;

              c = 0;
              for (x = 1; x < Kmer; x++)
                c = (c << 2) | s[p++];
              while (p < q)
                { x = s[p];
                  c = ((c << 2) | x) & Kmask;
                  list[n].read = i;
                  list[n].rpos = p++;
                  list[n].code = c;
                  n += 1;
                  kptr[c & BMASK] += 1;
                }
            }
          s += (q+1);
          io = 0;
        }
    }

  return (NULL);
}

static void *biased_tuple_thread(void *arg)
{ Tuple_Arg  *data  = (Tuple_Arg *) arg;
  DAZZ_READ  *reads = TA_block->reads;
  int64      *kptr  = data->kptr;
  KmerPos    *list  = TA_list;
  int         i, j;
  int         m, n;
  char       *s;

  n = data->tbeg;
  i = data->rbeg;
  m = data->tend;
  j = data->rend;
  s = ((char *) TA_block->bases) + reads[i].boff;

  if (TA_track != NULL)
    { int64     *anno1 = ((int64 *) (TA_track->anno)) + 1;
      int       *point = (int *) (TA_track->data);
      int        p, q;
      int        io, jo;
      int64      a, b, f; 
      char      *s1;

      io = n + i*Kmer - reads[i].boff;
      jo = TA_block->maxlen;
      q = 0;
      f = anno1[i-1];
      for (; i <= j; i++)
        { if (i == j)
            { jo = m + j*Kmer - reads[j].boff;
              if (jo == 0)
                continue;
              jo += Kmer-1;
            }
          b = f;
          f = anno1[i];
          s1 = s+1;
          for (a = b; a <= f; a += 2)
            { if (a == b)
                p = 0;
              else
                p = point[a-1];
              if (a == f)
                q = reads[i].rlen;
              else
                q = point[a];
              if (p < io)
                p = io;
              if (q > jo)
                q = jo;
              if (p+Kmer <= q)
                { int    x, k, z;
                  uint64 c, d;

                  c = 0;
                  z = 0;
                  k = 1;
                  while (p < q)
                    { x = s[p];
                      z += LogBase[x];
                      c  = ((c << 2) | x);
                      while (z < LogNorm && k < Kmer)
                        { if (++p >= q)
                            break;
                          k += 1;
                          x  = s[p];
                          z += LogBase[x];
                          c  = ((c << 2) | x);
                        }
                      while (1)
                        { int sc = z-LogBase[(int) s1[p-k]];
                          if (sc < LogNorm) break;
                          z  = sc;
                          k -= 1;
                        }
                      if (z > LogThresh)
                        { d = ((c << NormShift[k]) & Kmask);
                          list[n].read = i;
                          list[n].rpos = p;
                          list[n].code = d;
                          n += 1;
                          kptr[d & BMASK] += 1;
                        }
                      p += 1;
                      z -= LogBase[(int) s[p-k]];
                    }
                }
            }
          s += (q+1);
          io = 0;
	}
    }

  else
    { int    p, q;
      int    io;
      char  *s1;

      io = n + i*Kmer - reads[i].boff;
      for (; i <= j; i++)
        { s1 = s+1;
          p  = io;
          if (i == j)
            q = m + j*Kmer + (Kmer-1) - reads[j].boff;
          else
            q = reads[i].rlen;
          if (p+Kmer <= q)
            { int    x, k, z;
              uint64 c, d;

              c = 0;
              z = 0;
              k = 1;
              while (p < q)
                { x = s[p];
                  z += LogBase[x];
                  c  = ((c << 2) | x);
                  while (z < LogNorm && k < Kmer)
                    { if (++p >= q)
                        break;
                      k += 1;
                      x  = s[p];
                      z += LogBase[x];
                      c  = ((c << 2) | x);
                    }
                  while (1)
                    { int sc = z-LogBase[(int) s1[p-k]];
                      if (sc < LogNorm) break;
                      z  = sc;
                      k -= 1;
                    }
                  if (z > LogThresh)
                    { d = ((c << NormShift[k]) & Kmask);
                      list[n].read = i;
                      list[n].rpos = p;
                      list[n].code = d;
                      n += 1;
                      kptr[d & BMASK] += 1;
                    }
                  p += 1;
                  z -= LogBase[(int) s[p-k]];
                }
            }
          s += (q+1);
          io = 0;
        }
    }

  kptr[BMASK] += (data->fill = m-n);
  while (n < m)
    { list[n].code = 0xffffffffffffffffllu;
      list[n].read = 0xffffffff;
      list[n].rpos = 0xffffffff;
      n += 1;
    }

  return (NULL);
}

static KmerPos *FR_src;
static KmerPos *FR_trg;

typedef struct
  { int  beg;
    int  end;
    int  kept;
  } Comp_Arg;

static void *compsize_thread(void *arg)
{ Comp_Arg   *data  = (Comp_Arg *) arg;
  int         end   = data->end;
  KmerPos    *src   = FR_src;
  int         n, i, c, p;
  uint64      h, g;

  i = data->beg;
  h = src[i].code;
  n = 0;
  while (i < end)
    { p = i++;
      while ((g = src[i].code) == h)
        i += 1;
      if ((c = (i-p)) < TooFrequent)
        n += c;
      h = g;
    }

  data->kept = n;
  return (NULL);
}

static void *compress_thread(void *arg)
{ Comp_Arg   *data  = (Comp_Arg *) arg;
  int         end   = data->end;
  KmerPos    *src   = FR_src;
  KmerPos    *trg   = FR_trg;
  int         n, i, p;
  uint64      h, g;

  i = data->beg;
  h = src[i].code;
  n = data->kept;
  while (i < end)
    { p = i++;
      while ((g = src[i].code) == h)
        i += 1;
      if (i-p < TooFrequent)
        { while (p < i)
            trg[n++] = src[p++];
        }
      h = g;
    }

  return (NULL);
}

static int find_read(int x, DAZZ_READ *a, int n)
{ int64 l, r, m;

  // largest k s.t. a[k].boff - k*Kmer <= x (or n if does not exist)

  l = 0;
  r = n;
  while (l < r)
    { m = ((l+r+1) >> 1);
      if (a[m].boff - m*Kmer > x)
        r = m-1;
      else
        l = m;
    }
  return (l);
}

void *Sort_Kmers(DAZZ_DB *block, int *len)
{ THREAD    threads[NTHREADS];
  Tuple_Arg parmt[NTHREADS];
  Comp_Arg  parmf[NTHREADS];
  Lex_Arg   parmx[NTHREADS];
  int       mersort[16];

  KmerPos  *src, *trg, *rez;
  int       kmers, nreads;
  int       i, j, x, z;
  uint64    h;

  for (i = 0; i < NTHREADS; i++)
    parmx[i].sptr = (int64 *) alloca(NTHREADS*BPOWR*sizeof(int64));

  for (i = 0; i < 16; i++)
    mersort[i] = 0;
  for (i = 0; i < Kshift; i += 8)
    mersort[i>>3] = 1;

  if (NormShift == NULL && BIASED)
    { double scale;

      NormShift = (int *) Malloc(sizeof(int)*(Kmer+1),"Allocating Sort_Kmers bias shift");
      if (NormShift == NULL)
        exit (1);
      for (i = 0; i <= Kmer; i++)
        NormShift[i] = Kshift - 2*i;
      LogNorm = 10000 * Kmer;
      LogThresh = 10000 * (Kmer-MAX_BIAS);

      scale = -10000. / log(4.);
      for (i = 0; i < 4; i++)
        LogBase[i] = (int) ceil( scale * log(block->freq[i]) );
    }

  nreads = block->nreads;
  kmers  = block->reads[nreads].boff - Kmer * nreads;

  if (kmers <= 0)
    goto no_mers;

  if (( (Kshift-1)/BSHIFT + (TooFrequent < INT32_MAX) ) & 0x1)
    { trg = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+2),"Allocating Sort_Kmers vectors");
      src = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+2),"Allocating Sort_Kmers vectors");
    }
  else
    { src = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+2),"Allocating Sort_Kmers vectors");
      trg = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+2),"Allocating Sort_Kmers vectors");
    }
  if (src == NULL || trg == NULL)
    exit (1);

  if (VERBOSE)
    { printf("\n   Kmer count = ");
      Print_Number((int64) kmers,0,stdout);
      printf("\n   Using %.2fGb of space\n",(1. * kmers) / 33554432);
      fflush(stdout);
    }

  TA_block = block;
  TA_list  = src;
  TA_track = block->tracks;

  x = z = 0;
  for (i = 0; i < NTHREADS; i++)
    { parmt[i].kptr = parmx[i].tptr;
      for (j = 0; j < BPOWR; j++)
        parmt[i].kptr[j] = 0;
      parmt[i].tbeg = parmx[i].beg = x;
      parmt[i].rbeg = z;
      parmt[i].tend = parmx[i].end = x = (((int64) kmers)*(i+1) >> NSHIFT);
      parmt[i].rend = z = find_read(x,block->reads,nreads);
    }

  if (BIASED)
    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,biased_tuple_thread,parmt+i);
  else
    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,tuple_thread,parmt+i);

  for (i = 0; i < NTHREADS; i++)
    pthread_join(threads[i],NULL);

  rez = (KmerPos *) lex_sort(mersort,(Double *) src,(Double *) trg,parmx);

  if (BIASED || TA_track != NULL)
    for (i = 0; i < NTHREADS; i++)
      kmers -= parmt[i].fill;

  if (TooFrequent < INT32_MAX && kmers > 0)
    { parmf[0].beg = 0;
      for (i = 1; i < NTHREADS; i++)
        { x = (((int64) i)*kmers) >> NSHIFT;
          h = rez[x-1].code;
          while (rez[x].code == h)
            x += 1;
          parmf[i-1].end = parmf[i].beg = x;
        }
      parmf[NTHREADS-1].end = kmers;

      if (rez[kmers-1].code == 0xffffffffffffffffllu)
        rez[kmers].code = 0;
      else
        rez[kmers].code = 0xffffffffffffffffllu;

      if (src == rez)
        { FR_src = src;
          FR_trg = rez = trg;
        }
      else
        { FR_src = trg;
          FR_trg = rez = src;
        }

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,compsize_thread,parmf+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);

      x = 0;
      for (i = 0; i < NTHREADS; i++)
        { z = parmf[i].kept;
          parmf[i].kept = x;
          x += z;
        }
      kmers = x;

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,compress_thread,parmf+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);
    }

  rez[kmers].code   = 0xffffffffffffffffllu;
  rez[kmers+1].code = 0;
    
  if (src != rez)
    free(src);
  else
    free(trg);

#ifdef TEST_KSORT
  { int i;

    printf("\nKMER SORT:\n");
    for (i = 0; i < HOW_MANY && i < kmers; i++)
      { KmerPos *c = rez+i;
        printf(" %5d / %5d / %10lld\n",c->read,c->rpos,c->code);
      }
    fflush(stdout);
  }
#endif

  if (VERBOSE)
    { if (TooFrequent < INT32_MAX || BIASED || TA_track != NULL)
        { printf("   Revised kmer count = ");
          Print_Number((int64) kmers,0,stdout);
          printf("\n");
        }
      printf("   Index occupies %.2fGb\n",(1. * kmers) / 67108864);
      fflush(stdout);
    }

  if (kmers <= 0)
    { free(rez);
      goto no_mers;
    }

  if (kmers > (int64) (MEM_LIMIT/(4*sizeof(KmerPos))))
    { fprintf(stderr,"Warning: Block size too big, index occupies more than 1/4 of");
      if (MEM_LIMIT == MEM_PHYSICAL)
        fprintf(stderr," physical memory (%.1fGb)\n",(1.*MEM_LIMIT)/0x40000000ll);
      else
        fprintf(stderr," desired memory allocation (%.1fGb)\n",(1.*MEM_LIMIT)/0x40000000ll);
      fflush(stderr);
    }

  *len = kmers;
  return (rez);

no_mers:
  *len = 0;
  return (NULL);
}


/*******************************************************************************************
 *
 *  MERGE SORTED K-MER LISTS
 *
 ********************************************************************************************/

static int find_tuple(uint64 x, KmerPos *a, int n)
{ int l, r, m;

  // smallest k s.t. a[k].code >= x (or n if does not exist)

  l = 0;
  r = n;
  while (l < r)
    { m = ((l+r) >> 1);
      if (a[m].code < x)
        l = m+1;
      else
        r = m;
    }
  return (l);
}

static int find_aread(int x, SeedPair *a, int64 n)
{ int64 l, r, m;

  // smallest k s.t. a[k].aread >= x (or n if does not exist)

  l = 0;
  r = n;
  while (l < r)
    { m = ((l+r) >> 1);
      if (a[m].aread < x)
        l = m+1;
      else
        r = m;
    }
  return (l);
}

  //  Determine what *will* be the size of the merged list and histogram of sizes for given cutoffs

static KmerPos  *MG_alist;
static KmerPos  *MG_blist;
static SeedPair *MG_hits;
static int       MG_comp;

typedef struct
  { int    abeg, aend;
    int    bbeg, bend;
    int64 *kptr;
    int64  nhits;
    int    limit;
    int64  hitgram[MAXGRAM];
  } Merge_Arg;

static void *count_thread(void *arg)
{ Merge_Arg  *data  = (Merge_Arg *) arg;
  KmerPos    *asort = MG_alist;
  KmerPos    *bsort = MG_blist;
  int64      *gram  = data->hitgram;
  int64       nhits = 0;
  int         aend  = data->aend;

  int64  ct;
  int    ia, ib;
  int    jb, ja;
  uint64 ca, cb;
  uint64 da, db;

  ia = data->abeg;
  ca = asort[ia].code;
  ib = data->bbeg;
  cb = bsort[ib].code;
  while (1)
    { while (cb < ca)
        cb = bsort[++ib].code;
      while (cb > ca)
        ca = asort[++ia].code;
      if (cb == ca)
        { ja = ia++;
          while ((da = asort[ia].code) == ca)
            ia += 1;
          jb = ib++;
          while ((db = bsort[ib].code) == cb)
            ib += 1;

          if (ia > aend)
            { if (ja >= aend)
                break;
              da = asort[ia = aend].code;
              db = bsort[ib = data->bend].code;
            }

          ct  = (ia-ja);
          ct *= (ib-jb);

          nhits += ct;
          ca = da;
          cb = db;

          if (ct < MAXGRAM)
            gram[ct] += 1;
        }
    }

  data->nhits = nhits;

  return (NULL);
}

  //  Produce the merged list now that the list has been allocated and
  //    the appropriate cutoff determined.

static void *merge_thread(void *arg)
{ Merge_Arg  *data  = (Merge_Arg *) arg;
  int64      *kptr  = data->kptr;
  KmerPos    *asort = MG_alist;
  KmerPos    *bsort = MG_blist;
  SeedPair   *hits  = MG_hits;
  int64       nhits = data->nhits;
  int         aend  = data->aend;
  int         limit = data->limit;

  int64  ct;
  int    ia, ib;
  int    jb, ja;
  uint64 ca, cb;
  uint64 da, db;
  int    ap;
  int    a, b;

  ia = data->abeg;
  ca = asort[ia].code;
  ib = data->bbeg;
  cb = bsort[ib].code;
  while (1)
    { while (cb < ca)
        cb = bsort[++ib].code;
      while (cb > ca)
        ca = asort[++ia].code;
      if (cb == ca)
        { if (ia >= aend) break;
          ja = ia++;
          while ((da = asort[ia].code) == ca)
            ia += 1;
          jb = ib++;
          while ((db = bsort[ib].code) == cb)
            ib += 1;

          if (ia > aend)
            { if (ja >= aend)
                break;
              da = asort[ia = aend].code;
              db = bsort[ib = data->bend].code;
            }

          ct = ib-jb;
          if ((ia-ja)*ct < limit)
            { for (a = ja; a < ia; a++)
                { ap = asort[a].rpos;
                  kptr[ap & BMASK] += ct;
                  for (b = jb; b < ib; b++)
                    { hits[nhits].bread = bsort[b].read;
                      hits[nhits].aread = asort[a].read;
                      hits[nhits].apos  = ap;
                      hits[nhits].diag  = ap - bsort[b].rpos;
                      nhits += 1;
                    }
                }
            }
          ca = da;
          cb = db;
        }
    }

  return (NULL);
}


/*******************************************************************************************
 *
 *  GLOBAL VARS SEEN BY HIT HANDLER AND COMPONENTS
 *
 ********************************************************************************************/

static DAZZ_DB    *MR_ablock;
static DAZZ_DB    *MR_bblock;
static SeedPair   *MR_hits;
static Align_Spec *MR_spec;
static int         MR_start;
static int         MR_doA;
static int         MR_doB;


/*******************************************************************************************
 *
 *  SPLAY TREE K-MER CHAIN ROUTINES
 *
 ********************************************************************************************/

typedef struct _splay
  { struct _splay *lft;
    struct _splay *rgt;
    struct _splay *nxt;
    int            bmin;
    int            bmax;
    int            bpos;
    int            apos;
    int            diag;
    struct _splay *from;
    struct _splay *orig;
    int            cost;
  } Splay;

  //  Add node o to tree x.  Upon return, o is always the root of subtree
  //  The tree is ordered on the diag-field.

static inline void update(Splay *x)
{ Splay *l, *r;

  l = x->lft;
  r = x->rgt;
  if (r == NULL)
    if (l == NULL)
       x->bmin = x->bmax = x->bpos;
    else
      { if (x->bpos < l->bmin)
          x->bmin = x->bpos;
        else
          x->bmin = l->bmin;
        if (x->bpos > l->bmax)
          x->bmax = x->bpos;
        else
          x->bmax = l->bmax;
      }
  else
    if (l == NULL)
      { if (x->bpos < r->bmin)
          x->bmin = x->bpos;
        else
          x->bmin = r->bmin;
        if (x->bpos > r->bmax)
          x->bmax = x->bpos;
        else
          x->bmax = r->bmax;
      }
    else
      { if (x->bpos < l->bmin)
          if (x->bpos < r->bmin)
            x->bmin = x->bpos;
          else
            x->bmin = r->bmin;
        else
          if (l->bmin < r->bmin)
            x->bmin = l->bmin;
          else
            x->bmin = r->bmin;
        if (x->bpos > l->bmax)
          if (x->bpos > r->bmax)
            x->bmax = x->bpos;
          else
            x->bmax = r->bmax;
        else
          if (l->bmax > r->bmax)
            x->bmax = l->bmax;
          else
            x->bmax = r->bmax;
      }
}

static void add(Splay *x, Splay *o)
{ Splay *u;

  if (x == NULL)
    o->lft = o->rgt = NULL;
  else if (o->diag > x->diag || (o->diag == x->diag && o->apos > x->apos))
    { u = x->lft;
      if (u == NULL)
        { o->lft = NULL;
          o->rgt = x;
        }
      else if (o->diag > u->diag || (o->diag == u->diag && o->apos > u->apos))
        { add(u->lft,o);
          x->lft = u->rgt;
          u->rgt = x;
          u->lft = o->rgt;
          o->rgt = u;
          update(x);
          update(u);
        }
      else
        { add(u->rgt,o);
          u->rgt = o->lft;
          x->lft = o->rgt;
          o->lft = u;
          o->rgt = x;
          update(x);
          update(u);
        }
    }
  else
    { u = x->rgt;
      if (u == NULL)
        { o->lft = x;
          o->rgt = NULL;
        }
      else if (o->diag > u->diag || (o->diag == u->diag && o->apos > u->apos))
        { add(u->lft,o);
          x->rgt = o->lft;
          u->lft = o->rgt;
          o->lft = x;
          o->rgt = u;
          update(x);
          update(u);
        }
      else
        { add(u->rgt,o);
          x->rgt = u->lft;
          u->lft = x;
          u->rgt = o->lft;
          o->lft = u;
          update(x);
          update(u);
        }
    }
}

static Splay *del(Splay *x, Splay *o)
{ if (o == x)
    { Splay *l, *r;
      Splay *c, *d;

      l = o->lft;
      r = o->rgt;
      if (l != NULL)
        { if (l->rgt != NULL)
            { c = l;
              while (c->rgt->rgt != NULL)
                c = c->rgt;
              d = c->rgt;
              c->rgt = d->lft;
              add(l,d);
            }
          else
            d = l;
          d->rgt = r;
          update(d);
          return (d);
        }
      else if (r != NULL)
        { if (r->lft != NULL)
            { c = r;
              while (c->lft->lft != NULL)
                c = c->lft;
              d = c->lft;
              c->lft = d->rgt;
              add(r,d);
            }
          else
            d = r;
          d->lft = l;
          update(d);
          return (d);
        }
      else
        return (NULL);
    }
  else if (o->diag > x->diag || (o->diag == x->diag && o->apos > x->apos))
    { x->lft = del(x->lft,o);
      update(x);
      return (x);
    }
  else
    { x->rgt = del(x->rgt,o);
      update(x);
      return (x);
    }
}

  //  Convert tree h into a rgt-linked list that ends at to

static Splay *_linear(Splay *h, Splay *to)
{ Splay *head;

  if (h->lft != NULL)
    head = _linear(h->lft,h);
  else
    head = h;
  if (h->rgt != NULL)
    h->rgt = _linear(h->rgt,to);
  else
    h->rgt = to;
  return (head);
}

static Splay *linearize(Splay *h, Splay *e)
{ if (h == NULL)
    return (e);
  return (_linear(h,e));
}

  //  Diagnostic: print splay tree h at indent level lev.  ref is the base point of memory
  //     holding tree so that pointers can be printed as small integers

void print_tree(Splay *ref, Splay *h, int lev)
{ if (h == NULL) return;
  print_tree(ref,h->rgt,lev+2);
  printf("  %*s %4ld: %4ld < %8d(%5d,%8d) [%8d,%8d] > %4ld (->%4ld @%4ld)\n",
         lev,"",h-ref,(h->lft==NULL)?-1:h->lft-ref,h->diag,h->apos,h->bpos,h->bmin,h->bmax,
                      (h->rgt==NULL)?-1:h->rgt-ref,(h->from==NULL)?-1:h->from-ref,
                      (h->orig==NULL)?-1:h->orig-ref);
  print_tree(ref,h->lft,lev+2);
}

  //  Compress a chain (linked by from-field) by removing all same diagonal successors
  //    not further than 100bp away.  Return the length of the compressed chain.

static int chain_length(Splay *x)
{ Splay *y;
  int    n, da;

  n = 0;
  y = x->from;
  while (y != NULL)
    { da = x->apos - y->apos;
      if (da == x->bpos - y->bpos && da < 100)
        y = x->from = y->from;
      else
        { n += 1;
          x = y;
          y = x->from;
        }
    }
  return (n);
}

static Splay *predOf(Splay *x, int thresh)
{ x = x->lft;
  if (x == NULL || x->bmax < thresh)
    return (NULL);
  while (1)
    { if (x->rgt != NULL && x->rgt->bmax >= thresh)
        x = x->rgt;
      else if (x->bpos >= thresh)
        return (x);
      else
        { x = x->lft;
          if (x == NULL)
            return (NULL);
        }
    }
}

static Splay *leftmost(Splay *x, int diag, int thresh)
{ Splay *a;

  a = NULL;
  while (x != NULL)
    if (diag > x->diag)
      x = x->lft;
    else if (diag < x->diag)
      x = x->rgt;
    else
      { if (x->lft != NULL && x->lft->bmax >= thresh)
          { if (x->bpos >= thresh)
              a = x;
            x = x->lft;
          }
        else if (x->bpos >= thresh)
          return (x);
        else
          x = x->rgt;
      }
  return (a);
}
    

static Splay *succOf(Splay *x, int thresh)
{ x = x->rgt;
  if (x == NULL || x->bmin > thresh)
    return (NULL);
  while (1)
    { if (x->lft != NULL && x->lft->bmin <= thresh)
        x = x->lft;
      else if (x->bpos <= thresh)
        return (x);
      else
        { x = x->rgt;
          if (x == NULL)
            return (NULL);
        }
    }
}

  //  Print contents of splay tree h in decreasing order, suppressing any successors
  //    that have the same orig(in).  Ignore all entries with cost < hithr.  From the
  //    top level, start with orig = NULL.

Splay *print_list(Splay *h, int hithr, Splay *orig)
{ if (h == NULL)
    return (orig);

  orig = print_list(h->rgt,hithr,orig);

  if (h->cost >= hithr) // && h->orig != orig && h != orig)
    { int    da, db;
      double aspect, density;

      if (h->from == NULL)
        { printf("  %4d: (%5d,%5d)  (%8d,%8d) %8d",h->cost,h->apos,h->apos,
                                                           h->bpos,h->bpos,h->diag);
          aspect = 1.;
          density = 1.;
        }
      else
        { printf("  %4d: (%5d,%5d)  (%8d,%8d) %8d",h->cost,h->orig->apos,h->apos,
                                                         h->orig->bpos,h->bpos,h->diag);
          da = h->apos - h->orig->apos;
          db = h->bpos - h->orig->bpos;
          aspect = ( 1. * da) / db;
          if (da > db)
            density = h->cost / (1. * da);
          else
            density = h->cost / (1. * db);
        }
      printf(" <%8d,%8d>",h->bmin,h->bmax);
      printf("  %6.1f  %6.3f",aspect,density);
      if (h->from != NULL)
        { if (h->orig->orig == h)
            printf(" **");
        }
      else
        { if (h->orig == h)
            printf(" aa");
        }
      printf("\n");
      orig = h->orig;
    }

  orig = print_list(h->lft,hithr,orig);

  return (orig);
}


/*******************************************************************************************
 *
 *  EXAMINE AND CHAIN K-MER DAZZ
 *
 ********************************************************************************************/

  //  Analyze the k-mer matches between an A-read and all the possible B-reads that it
  //    could map to, finding and reporting the best non-redundant alignments.

typedef struct
  { int     next;
    uint16  adisp[5];
    uint16  bdisp[5];
  } Jump;

typedef struct
  { int  next;
    int  score;
    int  length;
    int  chain;
    int  bread;
    int  comp;
    int  afirst;
    int  alast;
    int  bfirst;
    int  blast;
  } Candidate;

void print_candidates(int aread, int cand, Candidate *cbase, Jump *jbase, int showJ)
{ int d;

  printf("Read %5d\n",aread);
  for (d = cand; d >= 0; d = cbase[d].next)
    { int     da, db;
      double  aspect, density;

      da = cbase[d].alast - cbase[d].afirst;
      db = cbase[d].blast - cbase[d].bfirst;
      aspect = ( 1. * da) / db;
      if (da > db)
        density = cbase[d].score / (1. * da);
      else
        density = cbase[d].score / (1. * db);

      printf("    %5d = %5d%c [%5d,%5d] vs [%8d,%8d]  |%4d|  %6.1f  %6.3f\n",
             cbase[d].score,cbase[d].bread,(cbase[d].comp?'c':'n'),
             cbase[d].afirst,cbase[d].alast,
             cbase[d].bfirst,cbase[d].blast,cbase[d].length,aspect,density);
      fflush(stdout);

      if (showJ)
        { int   j;
          int   n, m;

          m = cbase[d].length;
          for (j = cbase[d].chain; j >= 0; j = jbase[j].next)
            { if (m < 5)
                for (n = 0; n < m; n++)
                  printf("   (%4d,%4d)\n",jbase[j].adisp[n],jbase[j].bdisp[n]);
              else
                for (n = 0; n < 5; n++)
                  printf("   (%4d,%4d)\n",jbase[j].adisp[n],jbase[j].bdisp[n]);
              m -= 5;
            }
        }

      fflush(stdout);
    }
}

typedef struct
  { int64       hbeg, hend;
    int         abeg, aend;
    int64       nfilt;
    int64       ncheck;

    Candidate  *cbase;
    int         cmax;
    int         cfree;

    Jump       *jbase;
    int         jmax;
    int         jfree;

    int16      *cover;
    int         comax;

    Work_Data  *work;
    FILE       *afile;
    FILE       *bfile;
  } Report_Arg;

static void *chain_thread(void *arg)
{ Report_Arg  *data   = (Report_Arg *) arg;
  SeedPair    *hits   = MR_hits;
  DAZZ_READ   *aread  = MR_ablock->reads;
  int          bstart = MR_bblock->tfirst;
#if defined(TEST_CHAIN) || defined(TEST_CANDID)
  int          astart = MR_ablock->tfirst;
#endif

  Splay       *sbase, *sfree;    //  Splay space, allocated linearly

  int          cmax, cfree;      //  Candidate space, free list
  Candidate   *cbase;

  int          jmax, jfree;      //  K-mer position space, free list
  Jump        *jbase;

  int          comax;
  int16       *cover, *cnt;

  int64   nfilt = 0;

  int     hithr;
  int     ar, nar;
  int64   nidx, eidx;

#define INIT_FREE_SPACE(max,base,free,Type,FACTOR)					\
  max = data->max;									\
  if (MR_start)										\
    { int i;										\
											\
      if (max == 0)									\
        { max  = FACTOR*((data->aend - data->abeg) + 1);				\
          base = (Type *) Malloc(sizeof(Type)*max,"Allocating free space");		\
          if (base == NULL)								\
            exit (1);									\
        }										\
      else										\
        base = data->base;								\
      for (i = 0; i < max; i++)								\
        base[i].next = i+1;								\
      base[max-1].next = -1;								\
      free = 0;										\
    }											\
  else											\
    { base = data->base;								\
      free = data->free;								\
    }

#define EXPAND_FREE_SPACE(max,base,free,Type)					\
{ int i, est;									\
										\
  est  = (1.0*max*(data->aend-data->abeg)) / ((ar-data->abeg)+1.) + 2000;	\
  base = (Type *) Realloc(base,sizeof(Type)*est,"Reallocating free space");	\
  if (base == NULL)								\
    exit (1);									\
  for (i = max; i < est; i++)							\
    base[i].next = i+1;								\
  base[est-1].next = -1;							\
  free = max;									\
  max  = est;									\
}

#define FREE_LIST(x,free,base)				\
{ int j, k;						\
							\
  j = x;						\
  if (j >= 0)						\
    { for (k = base[j].next; k >= 0; k = base[k].next)	\
        j = k;						\
      base[j].next = free;				\
      free    = x;					\
    }							\
}

#define DEL_CELL(c,d,e,head,free,base)	\
  e = base[d].next;			\
  if (c < 0)				\
    head = e;				\
  else					\
    base[c].next = e;			\
  base[d].next = free;			\
  free = d;
  

  //  In ovl and align roles of A and B are reversed, as the B sequence must be the
  //    complemented sequence !!

  hithr  = HITMIN*Kmer;
  eidx   = data->hend;

  { Double *hitd = (Double *) MR_hits;
    uint64  npair, cpair;
    int     maxk;

    maxk = 0;
    nidx = data->hbeg;
    for (npair = hitd[nidx].p2; nidx < eidx; npair = cpair)
      { int64 f = nidx++;
        while ((cpair = hitd[nidx].p2) == npair)
          nidx += 1;
        if (nidx-f > maxk)
          maxk = (int) (nidx-f);
      }

    sbase = Malloc(sizeof(Splay)*(maxk+1),"Allocating splay free space\n");
    if (sbase == NULL)
      exit (1);

    cover = data->cover;
    comax = data->comax;
    if (MR_start)
      { int i, b;

        for (i = data->abeg; i < data->aend; i++)
          aread[i].coff = -1;

        if (PROFILE)
          { b = (aread[data->aend].boff-aread[data->abeg].boff)/SPACING + 2*(data->aend-data->abeg);
            if (b > comax)
              { comax = b;
                cover = (int16 *) Realloc(cover,sizeof(int16)*comax,"Allocating coverage vectors");
               }
            bzero(cover,b*sizeof(int16));
          }
      }

    INIT_FREE_SPACE(cmax,cbase,cfree,Candidate,8)
    INIT_FREE_SPACE(jmax,jbase,jfree,Jump,40)
  }

  //  For every k-mer hit in hits[nidx..eidx) do

  nidx = data->hbeg;
  nar  = hits[nidx].aread;
  cnt  = cover;
  if (PROFILE)
    { int x;

      for (x = data->abeg; x < nar; x++)
        cnt = cnt + (aread[x].rlen-1) / SPACING + 2;
    }

  for (ar = nar; nidx < eidx; ar = nar)
    { int        br, nbr;
      int        apos, bpos;
      int        ncan;
      Splay     *chain, *expired, *queue, *qlast;
      int        cand;

      //  Diagonal bucket size and range depends on A-reads length and average error rate

      sfree   = sbase;
      chain   = NULL;
      expired = NULL;
      queue   = NULL;
      qlast   = NULL;
      cand    = aread[ar].coff;
      ncan    = 0;

      br = hits[nidx].bread;
      while (1)
        { nar = hits[nidx].aread;
          nbr = hits[nidx].bread;

          //  If last hit with a given B-read, analyze significant K-mer chains as candidates

          if (nbr != br || nar != ar)

            { Splay *h;

              chain = linearize(chain,expired);

              //  For each significant K-mer chain do:

#ifdef TEST_CANDID
              printf("Read %d vs Contig %d\n",ar+astart,br+bstart); fflush(stdout);
#endif
              for (h = chain; h != NULL; h = h->rgt)
                if (h->cost >= hithr && h->orig->orig == h)
                  { int    ab, bb;
                    int    ae, be;

                    ab = h->orig->apos - Kmer;
                    bb = h->orig->bpos - Kmer;
                    ae = h->apos;
                    be = h->bpos;
#ifdef TEST_CANDID
                    printf("  %5d: (%5d,%5d)  (%8d,%8d)\n",h->cost,ab,ae,bb,be); fflush(stdout);
#endif

                    if (PROFILE)
                      { int tb, te;
                        int cb, ce;

                        tb = ab/SPACING;
                        te = (ae-1)/SPACING+1;
                        cb = cnt[tb];
                        ce = cnt[te];
                        if (cb < 0x7fff && ce >-0xffff)
                          { cnt[tb] = (int16) (cb+1);
                            cnt[te] = (int16) (ce-1);
                          }
                      }

                    { int d, c, e;

                      //  A seed c is dominated by another d if it does not contain at least
                      //    MIN_PIECE bp's not in d and its score is less than that of d.
                      //    Remove all seeds in candidate stack dominated by c and remove c
                      //    if it is dominated by a candidate in the stack.

                      c = -1;
                      for (d = cand; d >= 0; d = e)
                        { if (cbase[d].afirst < ab+MIN_PIECE && cbase[d].alast > ae-MIN_PIECE)
                            if (ab < cbase[d].afirst+MIN_PIECE && ae > cbase[d].alast-MIN_PIECE)
                              { if (.9 * cbase[d].score >= h->cost)
                                  break;
                                else if (cbase[d].score <= .9 * h->cost)
                                  { FREE_LIST(cbase[d].chain,jfree,jbase);
                                    DEL_CELL(c,d,e,cand,cfree,cbase);
                                    ncan -= 1;
                                  }
                                else
                                  { c = d;
                                    e = cbase[d].next;
                                  }
                              }
                            else
                              { if (.9 * cbase[d].score >= h->cost)
                                  break;
                                c = d;
                                e = cbase[d].next;
                              }
                          else
                            if (ab < cbase[d].afirst+MIN_PIECE && ae > cbase[d].alast-MIN_PIECE)
                              { if (cbase[d].score <= .9 * h->cost)
                                  { FREE_LIST(cbase[d].chain,jfree,jbase);
                                    DEL_CELL(c,d,e,cand,cfree,cbase);
                                    ncan -= 1;
                                  }
                                else
                                  { c = d;
                                    e = cbase[d].next;
                                  }
                              }
                            else
                              { c = d;
                                e = cbase[d].next;
                              }
                        }

                      if (d < 0)
                        { int len;

                          if (cfree < 0)
                            EXPAND_FREE_SPACE(cmax,cbase,cfree,Candidate)
                          d = cfree;
                          cfree = cbase[d].next;
                          cbase[d].next = cand;
                          cand  = d;
                          cbase[d].bread  = br + bstart;
                          cbase[d].comp   = MG_comp;
                          cbase[d].score  = h->cost;
                          cbase[d].afirst = ab;
                          cbase[d].alast  = ae;
                          cbase[d].bfirst = bb;
                          cbase[d].blast  = be;

                          cbase[d].length = len = chain_length(h);

                          if (len == 0)
                            cbase[d].chain = -1;
                          else
                            { Splay *f, *g;
                              int    j,  k;
                              int    n;

                              if (jfree < 0)
                                EXPAND_FREE_SPACE(jmax,jbase,jfree,Jump)
                              cbase[d].chain = j = jfree;
                              g = h;
                              n = 0;
                              for (f = h->from; f != NULL; f = f->from)
                                { if (n == 5) 
                                    { k = jbase[j].next;
                                      if (k < 0)
                                        { EXPAND_FREE_SPACE(jmax,jbase,jfree,Jump)
                                          k = jbase[j].next = jfree;
                                        }
                                      j = k;
                                      n = 0;
                                    }
                                  jbase[j].adisp[n]   = (uint16) (g->apos - f->apos); 
                                  jbase[j].bdisp[n++] = (uint16) (g->bpos - f->bpos); 
                                  g = f;
                                }
                              jfree         = jbase[j].next;
                              jbase[j].next = -1; 
                            }
     
                          ncan += 1;
                        }
                    }
                  }

              if (nar != ar)
                break;
              br = nbr;

              sfree   = sbase;
              chain   = NULL;
              expired = NULL;
              queue   = NULL;
            }

          //  Add next k-mer match to its best chain extension

          { Splay *l, *r;
            int    lcost, rcost;

            apos = hits[nidx].apos + 1;
            bpos = apos - hits[nidx].diag;

            while (queue != NULL && queue->apos < apos-MAX_GAP)
              { if (queue->bmin >= 0)
                  { chain = del(chain,queue);
                    if (queue->orig->orig == queue)
                      { queue->rgt = expired;
                        expired = queue;
                      }
                  }
                queue = queue->nxt;
              }

            sfree->apos = apos;
            sfree->bpos = bpos;
            sfree->diag = apos - bpos;

            add(chain,sfree);
            chain = sfree++;
            update(chain);

            l = predOf(chain,bpos-MAX_GAP);
            if (l != NULL)
              l = leftmost(chain,l->diag,bpos-MAX_GAP);
            r = succOf(chain,bpos);
            lcost = rcost = 0;
            if (l != NULL)
              { if (apos >= l->apos + Kmer)
                  lcost = l->cost + Kmer;
                else
                  lcost = l->cost + (apos - l->apos);
              }
            if (r != NULL)
              { if (bpos >= r->bpos + Kmer)
                  rcost = r->cost + Kmer;
                else
                  rcost = r->cost + (bpos - r->bpos);
              }
            if (lcost > rcost)
              rcost = 0;
            else
              lcost = 0;
                
            if (lcost > 0)
              { chain->from = l;
                chain->cost = lcost;
                if (l->from == NULL)
                  chain->orig = l;
                else
                  chain->orig = l->orig;
                if (lcost >= chain->orig->orig->cost)
                  { chain->orig->orig = chain;
                    if (abs(l->diag-chain->diag) <= .2*(chain->apos-l->apos))
                      { del(chain,l);
                        l->bmin = -1;
                      }
                  }
              }
            else if (rcost > 0)
              { chain->from = r;
                chain->cost = rcost;
                if (r->from == NULL)
                  chain->orig = r;
                else
                  chain->orig = r->orig;
                if (rcost >= chain->orig->orig->cost)
                  { chain->orig->orig = chain;
                    if (abs(r->diag-chain->diag) <= .2*(chain->apos-r->apos))
                      { del(chain,r);
                        r->bmin = -1;
                      }
                  }
              }
            else
              { chain->from = NULL;
                chain->cost = Kmer;
                chain->orig = chain;
              }

            chain->nxt = NULL;
            if (queue != NULL)
              qlast->nxt = chain;
            else
              queue = chain;
            qlast = chain;
          }

          nidx += 1;
        }

      aread[ar].coff = cand;
      nfilt += ncan;

#ifdef TEST_CHAIN
      print_candidates(ar+astart,cand,cbase,jbase,0);
      if (PROFILE)
        { int i, c, n;

          printf(" Profile: %ld\n",cnt-cover);
          n = (aread[ar].rlen-1) / SPACING + 1;
          c = 0;
          for (i = 0; i <= n; i++)
            if (cnt[i] != 0)
              { c += cnt[i];
                printf(" %5d: %3d\n",i*SPACING,c);
              }
        }
#endif

      if (PROFILE)
        { int x;

          if (nar >= data->aend)
            nar = data->aend;
          for (x = ar; x < nar; x++)
            cnt = cnt + (aread[x].rlen-1) / SPACING + 2;
        }
    }

  //  Clean up

  free(sbase);

  data->jmax  = jmax;
  data->jbase = jbase;
  data->jfree = jfree;

  data->cmax  = cmax;
  data->cbase = cbase;
  data->cfree = cfree;

  data->comax = comax;
  data->cover = cover;

  data->nfilt  = nfilt;

  return (NULL);
}


/*******************************************************************************************
 *
 *  EXAMINE ALIGNMENT DAZZ
 *
 ********************************************************************************************/

typedef struct
  { uint64   max;
    uint64   top;
    uint16  *trace;
  } Trace_Buffer;

  //  Reverse complement s into t padded up to have length divisible by trace spacing
  //    Return the padded length.

static void complement(char *t, char *s, int len)
{ int i;

  t += len;
  *t = 4;
  for (i = 0; i < len; i++)
    *--t = (char ) (3-s[i]);
  *--t = 4;
}

  //  Determine if jpath and kpath have a trace point in common, returning
  //     the a-coordinate in *where if so.

static int Entwine(Path *jpath, Path *kpath, Trace_Buffer *tbuf, int *where)
{ int   ac, b2, y2, ae;
  int   i, j, k;
  int   num, den, min;
#ifdef SEE_ENTWINE
  int   strt = 1;
  int   iflare, oflare;
#endif

  uint16 *ktrace = tbuf->trace + (uint64) (kpath->trace);
  uint16 *jtrace = tbuf->trace + (uint64) (jpath->trace);

  min   = 10000;
  num   = 0;
  den   = 0;

#ifdef SEE_ENTWINE
  printf("\n");
#endif

  y2 = jpath->bbpos;
  j  = jpath->abpos/SPACING;

  b2 = kpath->bbpos;
  k  = kpath->abpos/SPACING;

  if (jpath->abpos == kpath->abpos)
    { min = abs(y2-b2);
      if (min == 0)
        *where = kpath->abpos;
    }

  if (j < k)
    { ac = k*SPACING;

      j = 1 + 2*(k-j);
      k = 1;

      for (i = 1; i < j; i += 2)
        y2 += jtrace[i];
    }
  else
    { ac = j*SPACING;

      k = 1 + 2*(j-k);
      j = 1;

      for (i = 1; i < k; i += 2)
        b2 += ktrace[i];
    }

  ae = jpath->aepos;
  if (ae > kpath->aepos)
    ae = kpath->aepos;

  while (1)
    { ac += SPACING;
      if (ac >= ae)
        break;
      y2 += jtrace[j];
      b2 += ktrace[k];
      j += 2;
      k += 2;

#ifdef SEE_ENTWINE
      printf("   @ %5d : %5d %5d = %4d\n",ac,y2,b2,abs(b2-y2));
#endif

      i = abs(y2-b2);
      if (i <= min)
        { min = i;
          if (i == 0)
            *where = ac;
        }
      num += i;
      den += 1;
#ifdef SEE_ENTWINE
      if (strt)
        { strt   = 0;
          iflare = i;
        }
      oflare = i;
#endif
    }

  if (jpath->aepos == kpath->aepos)
    { i = abs(jpath->bepos-kpath->bepos);
      if (i <= min)
        { min = i;
          if (i == 0)
            *where = kpath->aepos;
        }
    }

#ifdef SEE_ENTWINE
  if (den == 0)
    printf("Nothing\n");
  else
    printf("MINIM = %d AVERAGE = %d  IFLARE = %d  OFLARE = %d\n",min,num/den,iflare,oflare);
#endif

  if (den == 0)
    return (-1);
  else
    return (min);
}


  //  Produce the concatentation of path1 and path2 where they are known to meet at
  //    the trace point with coordinate ap. Place this result in a big growing buffer,
  //    that gets reset when fusion is called with path1 = NULL

static void Fusion(Path *path1, int ap, Path *path2, Trace_Buffer *tbuf)
{ int     k, k1, k2;
  int     len, diff;
  uint16 *trace;

  k1 = 2 * ((ap/SPACING) - (path1->abpos/SPACING));
  k2 = 2 * ((ap/SPACING) - (path2->abpos/SPACING));

  len = k1+(path2->tlen-k2);

  if (tbuf->top + len >= tbuf->max)
    { tbuf->max = 1.2*(tbuf->top+len) + 1000;
      tbuf->trace = (uint16 *) Realloc(tbuf->trace,sizeof(uint16)*tbuf->max,"Allocating paths");
      if (tbuf->trace == NULL)
        exit (1);
    }

  trace = tbuf->trace + tbuf->top;
  tbuf->top += len;

  diff = 0;
  len  = 0;
  if (k1 > 0)
    { uint16 *t = tbuf->trace + (uint64) (path1->trace);
      for (k = 0; k < k1; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
          diff += t[k];
        }
    }
  if (k2 < path2->tlen)
    { uint16 *t = tbuf->trace + (uint64) (path2->trace);
      for (k = k2; k < path2->tlen; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
          diff += t[k];
        }
    }

  path1->aepos = path2->aepos;
  path1->bepos = path2->bepos;
  path1->diffs = diff;
  path1->trace = (void *) (trace - tbuf->trace);
  path1->tlen  = len;
}


  //  Fuse any alignments in amatch[0..novls) that share a trace point, and eliminate
  //    any alignments whose a-span is entired contained in another.  Return the #
  //    of alignment remaining.

static int Handle_Redundancies(Overlap *amatch, int novls, 
                               Overlap *bmatch, Trace_Buffer *tbuf, int cm)
{ Path *jpath, *kpath;
  Path *jmath, *kmath;
  int   j, k, no;
  int   dist;
  int   awhen = 0, bwhen = 0;
  int   hasB;

#ifdef TEST_CONTAIN
  for (j = 0; j < novls; j++)
    printf("  %3d: [%5d,%5d] x [%5d,%5d]\n",j,amatch[j].path.abpos,amatch[j].path.aepos,
                                              amatch[j].path.bbpos,amatch[j].path.bepos);
#endif

  hasB = (bmatch != NULL);

  for (j = 1; j < novls; j++)
    { jpath = &(amatch[j].path);
      jmath = &(bmatch[j].path);
      for (k = j-1; k >= 0; k--)
        { kpath = &(amatch[k].path);
          kmath = &(bmatch[k].path);

          if (kpath->abpos < 0)
            continue;

          if (jpath->abpos < kpath->abpos)

            { if (kpath->abpos <= jpath->aepos && kpath->bbpos <= jpath->bepos)
                { dist = Entwine(jpath,kpath,tbuf,&awhen);
                  if (dist == 0)
                    { if (kpath->aepos > jpath->aepos)
                        { if (hasB)
                            { if (cm)
                                { dist = Entwine(kmath,jmath,tbuf,&bwhen);
                                  if (dist != 0)
                                    continue;
                                  Fusion(jpath,awhen,kpath,tbuf);
                                  Fusion(kmath,bwhen,jmath,tbuf);
                                  *jmath = *kmath;
#ifdef TEST_CONTAIN
                                  printf("  Really 1");
#endif
                                }
                              else
                                { dist = Entwine(jmath,kmath,tbuf,&bwhen);
                                  if (dist != 0)
                                    continue;
                                  Fusion(jpath,awhen,kpath,tbuf);
                                  Fusion(jmath,bwhen,kmath,tbuf);
#ifdef TEST_CONTAIN
                                  printf("  Really 2");
#endif
                                }
                            }
                          else
                            { Fusion(jpath,awhen,kpath,tbuf);
#ifdef TEST_CONTAIN
                              printf("  Really 3");
#endif
                            }
                          k = j;
                        }
                      kpath->abpos = -1;
#ifdef TEST_CONTAIN
                      printf("  Fuse! A %d %d\n",j,k);
#endif
                      break;
                    }
                }
            }

          else // kpath->abpos <= jpath->abpos

            { if (jpath->abpos <= kpath->aepos && jpath->bbpos <= kpath->bepos)
                { dist = Entwine(kpath,jpath,tbuf,&awhen);
                  if (dist == 0)
                    { if (kpath->abpos == jpath->abpos)
                        { if (kpath->aepos >jpath->aepos)
                            { *jpath = *kpath;
                              if (hasB)
                                *jmath = *kmath;;
                            }
                        }
                      else if (jpath->aepos > kpath->aepos)
                        { if (hasB)
                            { if (cm)
                                { dist = Entwine(jmath,kmath,tbuf,&bwhen);
                                  if (dist != 0)
                                    continue;
                                  Fusion(kpath,awhen,jpath,tbuf);
                                  *jpath = *kpath;
                                  Fusion(jmath,bwhen,kmath,tbuf);
#ifdef TEST_CONTAIN
                                  printf("  Really 4");
#endif
                                }
                              else
                                { dist = Entwine(kmath,jmath,tbuf,&bwhen);
                                  if (dist != 0)
                                    continue;
                                  Fusion(kpath,awhen,jpath,tbuf);
                                  *jpath = *kpath;
                                  Fusion(kmath,bwhen,jmath,tbuf);
                                  *jmath = *kmath;
#ifdef TEST_CONTAIN
                                  printf("  Really 5");
#endif
                                }
                            }
                          else
                            { Fusion(kpath,awhen,jpath,tbuf);
                              *jpath = *kpath;
#ifdef TEST_CONTAIN
                              printf("  Really 6");
#endif
                            }
                          k = j;
                        }
                      else
                        { *jpath = *kpath;
                          if (hasB)
                            *jmath = *kmath;
                        }
                      kpath->abpos = -1;
#ifdef TEST_CONTAIN
                      printf("  Fuse! B %d %d\n",j,k);
#endif
                      break;
                    }
                }
            }
        }
    }

  no = 0;
  for (j = 0; j < novls; j++)
    if (amatch[j].path.abpos >= 0)
      { if (hasB)
          bmatch[no] = bmatch[j];
        amatch[no++] = amatch[j];
      }
  novls = no;

#ifdef TEST_CONTAIN
  for (j = 0; j < novls; j++)
    printf("  %3d: [%5d,%5d] x [%5d,%5d]\n",j,amatch[j].path.abpos,amatch[j].path.aepos,
                                              amatch[j].path.bbpos,amatch[j].path.bepos);
#endif

  return (novls);
}

int special_log(int cover)
{ static double Spow[51];
  static int    first = 1;

  int    l, r, m;
  double x;

  if (first)
    { first = 0;
      for (m = 0; m <= 40; m++)
        Spow[m] = pow(10.,m/10.);
    }
 
  if (cover <= 1)
    return (cover);
  else if (cover >= 10000)
    return (40);

  // smallest k s.t. Spow[k] > cover (or 41 if does not exist)

  x = cover;
  l = 0;
  r = 41;
  while (l < r)
    { m = ((l+r) >> 1);
      if (Spow[m] <= x)
        l = m+1;
      else
        r = m;
    }

  return (l-1);
}

static int AMATCH_SORT(const void *l, const void *r)
{ Overlap *x = (Overlap *) l;
  Overlap *y = (Overlap *) r;
  if (y->path.abpos != x->path.abpos)
    return (y->path.abpos - x->path.abpos);
  else if (y > x)
    return (1);
  else if (y < x)
    return (-1);
  else
    return (0);
}

static int BN_MATCH_SORT(const void *l, const void *r)
{ Overlap *x = (Overlap *) l;
  Overlap *y = (Overlap *) r;
  if (y->path.bbpos != x->path.bbpos)
    return (y->path.bbpos - x->path.bbpos);
  else if (y > x)
    return (1);
  else if (y < x)
    return (-1);
  else
    return (0);
}

static int BC_MATCH_SORT(const void *l, const void *r)
{ Overlap *x = (Overlap *) l;
  Overlap *y = (Overlap *) r;
  if (y->path.bepos != x->path.bepos)
    return (x->path.bepos - y->path.bepos);
  else if (y > x)
    return (1);
  else if (y < x)
    return (-1);
  else
    return (0);
}

typedef struct
  { int   score;
    int   link;
    int   mark;
  } Links;

typedef struct
  { int beg;
    int end;
    int top;
  } Zones;

static int LINK_SORT(const void *l, const void *r)
{ Links *x = *((Links **) l);
  Links *y = *((Links **) r);

  return (y->score - x->score);
}

static void *report_thread(void *arg)
{ Report_Arg  *data   = (Report_Arg *) arg;
  Work_Data   *work  = data->work;
  FILE        *afile = data->afile;
  FILE        *bfile = data->bfile;
  Candidate   *cbase = data->cbase;
  Jump        *jbase = data->jbase;

  DAZZ_READ   *aread  = MR_ablock->reads;
  DAZZ_READ   *bread  = MR_bblock->reads;
  int          astart = MR_ablock->tfirst;
  char        *aseq   = (char *) (MR_ablock->bases);
  char        *bseq   = (char *) (MR_bblock->bases);

  Path        _path;
  Alignment   _align, *align = &_align;
  Path        *apath  = &_path;
  Path        *bpath;

  int64  ahits, bhits;
  int    ar, br, cm, c, d;
  int    j, n, m, t;
  int    apos, bpos, alast;
  int    alen, blen;
  int    hascomp;
  char  *acomp;

  int16 *cnt;
  uint8 *log;
  int    atck;

  int   small, tbytes, hithr;

  int      Omax;
  int      novl, lovl;
  Overlap *amatch, *bmatch;

  int      COmax;
  int      nparts;
  Links   *linker, **perm;
  Zones   *part;

  Trace_Buffer _tbuf, *tbuf = &_tbuf;

  ahits = 0;
  bhits = 0;
  align->path = apath;
  acomp = New_Read_Buffer(MR_ablock);

  if (SPACING <= TRACE_XOVR)
    { small  = 1;
      tbytes = sizeof(uint8);
    }
  else
    { small  = 0;
      tbytes = sizeof(uint16);
    }
  hithr = HITMIN*Kmer;

  if (MR_doA)
    { fwrite(&ahits,sizeof(int64),1,afile);
      fwrite(&SPACING,sizeof(int),1,afile);
    }
  if (MR_doB)
    { fwrite(&bhits,sizeof(int64),1,bfile);
      fwrite(&SPACING,sizeof(int),1,bfile);
    }

  Omax = COmax = MATCH_CHUNK;
  amatch = Malloc(sizeof(Overlap)*Omax,"Allocating match vector");
  if (MR_doB)
    bmatch = Malloc(sizeof(Overlap)*Omax,"Allocating match vector");
  else
    bmatch = amatch;
  linker = Malloc(sizeof(Links)*COmax,"Allocating link vector");
  perm   = Malloc(sizeof(Links *)*COmax,"Allocating permutation vector");
  part   = Malloc(sizeof(Zones)*COmax,"Allocating zones vector");

  tbuf->max   = 2*TRACE_CHUNK;
  tbuf->trace = Malloc(sizeof(short)*tbuf->max,"Allocating trace vector");

  if (amatch == NULL || bmatch == NULL || linker == NULL ||
      perm == NULL || part == NULL || tbuf->trace == NULL)
    exit (1);

  cnt = data->cover;
  log = (uint8 *) data->cover;
  for (ar = data->abeg; ar < data->aend; ar++)
    { hascomp = 0;
      alen    = aread[ar].rlen;
      atck    = (alen-1)/SPACING + 1; 
      novl    = 0;
      lovl    = 0;
      tbuf->top = 0;

#if defined(TEST_ALIGN) || defined(TEST_CONTAIN) || defined(TEST_SELECT)
      printf("\nRead %d:\n",ar+astart);
#endif
      for (c = aread[ar].coff; c >= 0; c = d)
        { br = cbase[c].bread;
          cm = cbase[c].comp;

#ifdef TEST_ALIGN
          printf("    %5d = %5d%c [%5d,%5d] vs [%8d,%8d]  |%4d|\n",
                 cbase[c].score,br,(cm?'c':'n'),cbase[c].afirst,cbase[c].alast,
                 cbase[c].bfirst,cbase[c].blast,cbase[c].length);
          fflush(stdout);
#endif

          align->bseq = bseq + bread[br].boff;
          align->blen = blen = bread[br].rlen;
          align->alen = alen;
          if (cm)
            { if ( ! hascomp)
                { complement(acomp,aseq+aread[ar].boff,alen);
                  hascomp = 1;
                }
              align->flags = ACOMP_FLAG;
              align->aseq  = acomp;
            }
          else
            { align->flags = 0;
              align->aseq  = aseq + aread[ar].boff;
            }

          apos  = cbase[c].alast;
          bpos  = cbase[c].blast;
          alast = alen + 1;
          m = cbase[c].length;
          for (j = cbase[c].chain; j >= 0; j = jbase[j].next)
            { t = 5;
              if (m < 5)
                t = m;
              for (n = 0; n < t; n++)
                { apos -= jbase[j].adisp[n];
                  bpos -= jbase[j].bdisp[n];
                  if (apos < alast)
                    { int dg, ad;

                      if (cm)
                        { int ac, bc;

                          ac = alen - apos;
                          bc = blen - bpos;
                          dg = ac-bc;
                          ad = ac+bc;
                        }
                      else
                        { dg = apos-bpos;
                          ad = apos+bpos;
                        }
                      bpath = Local_Alignment(align,work,MR_spec,dg,dg,ad,-1,-1);
                      if (apath->aepos - apath->abpos >= hithr)
                        { alast  = apath->abpos;
#ifdef TEST_ALIGN
                          if (cm)
                            printf("      Overlap %5dc",br);
                          else
                            printf("      Overlap %5dn",br);
                          printf(" [%5d,%5d] x [%5d,%5d] = %4d",
                                 apath->abpos,apath->aepos,apath->bbpos,apath->bepos,apath->diffs);
                          printf(" >> %d <<\n",alast);
#endif

                          if (novl >= Omax)
                            { Omax = 1.2*novl + MATCH_CHUNK;
                              amatch = Realloc(amatch,sizeof(Overlap)*Omax,
                                               "Reallocating match vector");
                              if (amatch == NULL)
                                exit (1);
                              if (MR_doB)
                                { bmatch = Realloc(bmatch,sizeof(Overlap)*Omax,
                                                          "Reallocating match vector");
                                  if (bmatch == NULL)
                                    exit (1);
                                }
                            }

                          amatch[novl].aread = ar + astart;
                          amatch[novl].bread = br;
                          amatch[novl].flags = (cm ? COMP_FLAG : 0);
                          amatch[novl].path  = *apath;
                       
                          if (tbuf->top + apath->tlen > tbuf->max)
                            { tbuf->max = 1.2*(tbuf->top+apath->tlen) + TRACE_CHUNK;
                              tbuf->trace = Realloc(tbuf->trace,sizeof(short)*tbuf->max,
                                                    "Reallocating trace vector");
                              if (tbuf->trace == NULL)
                                exit (1);
                            }
                          amatch[novl].path.trace = (void *) (tbuf->top);
                          memmove(tbuf->trace+tbuf->top,apath->trace,sizeof(short)*apath->tlen);
                          tbuf->top += apath->tlen;

                          if (MR_doB)
                            { bmatch[novl].aread = br;
                              bmatch[novl].bread = ar + astart;
                              bmatch[novl].flags = (cm ? COMP_FLAG : 0);
                              bmatch[novl].path  = *bpath;

                              if (tbuf->top + bpath->tlen > tbuf->max)
                                { tbuf->max = 1.2*(tbuf->top+bpath->tlen) + TRACE_CHUNK;
                                  tbuf->trace = Realloc(tbuf->trace,sizeof(short)*tbuf->max,
                                                        "Reallocating trace vector");
                                  if (tbuf->trace == NULL)
                                    exit (1);
                                }
                              bmatch[novl].path.trace = (void *) (tbuf->top);
                              memmove(tbuf->trace+tbuf->top,bpath->trace,sizeof(short)*bpath->tlen);
                              tbuf->top += bpath->tlen;
                            }

                          novl += 1;
                        }
                    }
                }
              m -= 5;
            }

          d = cbase[c].next;
          if (d < 0 || cbase[d].bread != br || cbase[d].comp != cm)
            {
#ifdef TEST_CONTAIN
              if (novl > 1)
                printf("\n%5d vs %5d:\n",ar+astart,br);
#endif

              //  Remove redundant alignment and fuse entwined alignments

              if (novl-lovl > 1)
                { if (MR_doB)
                    novl = lovl + Handle_Redundancies(amatch+lovl,novl-lovl,bmatch+lovl,tbuf,cm);
                  else
                    novl = lovl + Handle_Redundancies(amatch+lovl,novl-lovl,NULL,tbuf,cm);
                }

              if (novl - lovl > 1)
                { qsort(amatch+lovl,novl-lovl,sizeof(Overlap),AMATCH_SORT);
                  if (MR_doB)
                    { if (cm)
                        qsort(bmatch+lovl,novl-lovl,sizeof(Overlap),BC_MATCH_SORT);
                      else
                        qsort(bmatch+lovl,novl-lovl,sizeof(Overlap),BN_MATCH_SORT);
                    }
                }
  
              lovl = novl;
            }
        }

      if (novl >= COmax)
        { COmax  = 1.2*novl + MATCH_CHUNK;
          linker = Realloc(linker,sizeof(Links)*COmax,"Reallocating link vector");
          perm   = Realloc(perm,sizeof(Links *)*COmax,"Reallocating permutation vector");
          part   = Realloc(part,sizeof(Zones)*COmax,"Reallocating zones vector");
          if (linker == NULL || perm == NULL || part == NULL)
            exit (1);
        }

#ifdef TEST_SELECT
      for (c = 0; c < novl; c++)
        printf("  %3d: vs %5d [%5d,%5d] x [%5d,%5d] %.1f\n",c,amatch[c].bread,
               amatch[c].path.abpos,amatch[c].path.aepos,amatch[c].path.bbpos,amatch[c].path.bepos,
               (100.*amatch[c].path.diffs)/(amatch[c].path.aepos - amatch[c].path.abpos));
      printf("\n");
      fflush(stdout);
#endif

      lovl = 0;
      linker[0].link = -1;
      linker[0].score = (amatch[0].path.aepos - amatch[0].path.abpos) 
                      - DIFF_SCORE * amatch[0].path.diffs;
      linker[0].mark = 1;
      perm[0] = linker;
      br = amatch[0].bread;
      for (c = 1; c < novl; c++) 
        { Path *cpath;
          int   cor, dor;

          cpath = &(amatch[c].path);

          linker[c].link  = -1;
          linker[c].score = (cpath->aepos - cpath->abpos) - DIFF_SCORE * cpath->diffs;
          linker[c].mark  = 1;
          perm[c] = linker+c;

          if (amatch[c].bread != br) 
            { br = amatch[c].bread;
              lovl = c;
              continue;
            }

          cor = COMP(amatch[c].flags);
          for (d = c-1; d >= lovl; d--)
            if ((dor = COMP(amatch[d].flags)) == cor)
              { Path  *dpath = &(amatch[d].path);
                int    scr, scr2;
                int    gap, gap2;
                double rat;

                if (dor)
                  { if (dpath->bepos < cpath->bepos)
                      continue;
                  }
                else
                  { if (dpath->bbpos < cpath->bbpos)
                      continue;
                  }

                if (dpath->abpos <= cpath->aepos - CHAIN_OVL ||
                    dpath->bbpos <= cpath->bepos - CHAIN_OVL)
                  continue;

                rat = ( (dpath->abpos - cpath->aepos + CHAIN_OFF) /
                        (dpath->bbpos - cpath->bepos + CHAIN_OFF) );

                if (1. > rat*CHAIN_PLAY || rat > CHAIN_PLAY)
                  continue;

                scr  = linker[d].score + (cpath->aepos - cpath->abpos)
                                       - DIFF_SCORE * cpath->diffs;
                scr2 = linker[c].score;
                if (scr < scr2 - TIE_SCORE)
                  continue;

                if (scr <= scr2 + TIE_SCORE)
                  { gap = dpath->abpos - cpath->aepos;
                    if (linker[d].link >= 0)
                      gap2 = amatch[linker[c].link].path.aepos - dpath->abpos;
                    else
                      gap2 = 0;
                    if (gap > gap2 + TIE_GAP)
                      continue;

                    if (gap >= gap2 - TIE_GAP)
                      { if (scr < scr2)
                          continue;
                        if (scr == scr2)
                          { if (gap >= gap2)
                              continue;
                          }
                      }
                  }

                linker[c].link  = d;
                linker[c].score = scr;
                linker[d].mark  = 0;
              }
        }

      qsort(perm,novl,sizeof(Links *),LINK_SORT);

      nparts = 0;
      for (c = 0; c < novl && perm[c]->score >= 0; c++)
        if (perm[c]->mark == 1)
          { int p, b, e, q, n;
            int best;

            b = e = perm[c] - linker;
            for (p = linker[b].link; p >= 0 && linker[p].mark >= 0; p = linker[p].link)
              e = p;

            for (p = 0; p < nparts; p++)
              if (amatch[b].path.abpos < part[p].end-100 && amatch[e].path.aepos > part[p].beg+100)
                break;
            if (p >= nparts)
              { part[p].beg = amatch[b].path.abpos;
                part[p].end = amatch[e].path.aepos;
                part[p].top = linker[b].score;
                best = 1;
                nparts += 1;
              }
            else
              { if (linker[b].score < BEST_TIE * part[p].top)
                  continue;
                best = (linker[b].score == part[p].top);
              }

            q = -1;
            for (p = b; 1; p = n)
              { linker[p].mark = -1;
                if (MR_doA)
                  { amatch[p].path.trace = tbuf->trace + (uint64) (amatch[p].path.trace);
                    if (small)
                      Compress_TraceTo8(amatch+p);
                    if (p == b)
                      { amatch[p].flags |= START_FLAG;
                        if (best)
                          amatch[p].flags |= BEST_FLAG;
                      }
                    else
                      amatch[p].flags |= NEXT_FLAG;
                    Write_Overlap(afile,amatch+p,tbytes);
                    ahits += 1;
                  }
                n = linker[p].link;
                if (MR_doB)
                  { if (COMP(bmatch[p].flags))
                      { linker[p].link = q;
                        q = p;
                      } 
                    else
                      { bmatch[p].path.trace = tbuf->trace + (uint64) (bmatch[p].path.trace);
                        if (small)
                          Compress_TraceTo8(bmatch+p);
                        if (p == b)
                          { bmatch[p].flags |= START_FLAG;
                            if (best)
                              bmatch[p].flags |= BEST_FLAG;
                          }
                        else
                          bmatch[p].flags |= NEXT_FLAG;
                        Write_Overlap(bfile,bmatch+p,tbytes);
                        bhits += 1;
                      }
                  }
#ifdef TEST_SELECT
                if (p == b)
                  printf("  %5d : %5d ",amatch[p].bread,linker[p].score);
                else
                  printf(" -> ");
                printf("[%5d,%5d] %5d %.1f",
                       amatch[p].path.abpos,amatch[p].path.aepos,
                       (amatch[p].path.aepos-amatch[p].path.abpos),
                       (100.*amatch[p].path.diffs)/(amatch[p].path.aepos-amatch[p].path.abpos));
                fflush(stdout);
#endif
                if (p == e)
                  break;
              }
#ifdef TEST_SELECT
            printf("\n");
            fflush(stdout);
#endif
            if (MR_doB && COMP(bmatch[b].flags))
              { e = b;
                b = q;
                for (p = b; 1; p = linker[p].link)
                  { bmatch[p].path.trace = tbuf->trace + (uint64) (bmatch[p].path.trace);
                    if (small)
                      Compress_TraceTo8(bmatch+p);
                    if (p == b)
                      { bmatch[p].flags |= START_FLAG;
                        if (best)
                          bmatch[p].flags |= BEST_FLAG;
                      }
                    else
                      bmatch[p].flags |= NEXT_FLAG;
                    Write_Overlap(bfile,bmatch+p,tbytes);
                    bhits += 1;
                    if (p == e)
                      break;
                  }
              }
          }

#ifdef TEST_SELECT
      for (d = 0; d < nparts; d++)
        printf("  %2d: [%5d,%5d] %d\n",d,part[d].beg,part[d].end,part[d].top);
      if (PROFILE)
        { int i, c;

          printf(" Profile: %ld\n",cnt-data->cover);
          c = 0;
          for (i = 0; i <= atck; i++)
            if (cnt[i] != 0)
              { c += cnt[i];
                printf(" %5d: %3d\n",i*SPACING,c);
              }
        }
      fflush(stdout);
#endif

      if (PROFILE)
        { int i, c;

          c = 0;
          for (i = 0; i <= atck; i++)
            { c += cnt[i];
              log[i] = special_log(c);
            }
          cnt = cnt + atck + 1;
          log = log + atck + 1;
        }
    }

  free(tbuf->trace);
  free(part);
  free(perm);
  free(linker);
  if (MR_doB)
    free(bmatch);
  free(amatch);
  free(acomp-1);
 
  if (MR_doA)
    { rewind(afile);
      fwrite(&ahits,sizeof(int64),1,afile);
      data->ncheck = ahits;
    }
  else
    data->ncheck = bhits;

  if (MR_doB)
    { rewind(bfile);
      fwrite(&bhits,sizeof(int64),1,bfile);
    }

  return (NULL);
}


/*******************************************************************************************
 *
 *  THE ALGORITHM
 *    Match_Filter is called with each block of the reference (b-block) against a target
 *    block (a-block) with start non-zero for the first block of the reference.  Significant
 *    chains are cached and updates as all blocks of the reference are processed, and then
 *    Reporter is called to convert the significant chains into alignments and decide
 *    upon the final hits reported for each target read.
 *
 ********************************************************************************************/

static Report_Arg *parmr;
static int         tfilt;
static int         firstime = 1;

void Match_Filter(DAZZ_DB *ablock, DAZZ_DB *bblock,
                  void *vasort, int alen, void *vbsort, int blen, int comp, int start)

{ THREAD     threads[NTHREADS];
  Merge_Arg  parmm[NTHREADS];
  Lex_Arg    parmx[NTHREADS];
  int        pairsort[16];

  SeedPair *khit, *hhit;
  SeedPair *work1, *work2;
  int64     nhits;
  int64     nfilt;

  KmerPos  *asort, *bsort;
  int64     atot, btot;

  asort = (KmerPos *) vasort;
  bsort = (KmerPos *) vbsort;

  atot = ablock->totlen;
  btot = bblock->totlen;

  { int64 powr;
    int   i, nbyte;

    for (i = 0; i < NTHREADS; i++)
      parmx[i].sptr = (int64 *) alloca(NTHREADS*BPOWR*sizeof(int64));

    for (i = 0; i < 16; i++)
      pairsort[i] = 0;

    powr = 1;
    for (nbyte = 0; powr < ablock->maxlen; nbyte += 1)
      powr <<= 8;
    for (i = 4; i < 4+nbyte; i++)
      pairsort[i] = 1;

    powr = 1;
    for (nbyte = 0; powr < bblock->nreads; nbyte += 1)
      powr <<= 8;
    for (i = 8; i < 8+nbyte; i++)
      pairsort[i] = 1;

    powr = 1;
    for (nbyte = 0; powr < ablock->nreads; nbyte += 1)
      powr <<= 8;
    for (i = 12; i < 12+nbyte; i++)
      pairsort[i] = 1;
  }

  if (firstime)
    { int i;

      firstime = 0;
      parmr = (Report_Arg *) Malloc(sizeof(Report_Arg)*NTHREADS,"Allocating chain thread records");
      for (i = 0; i < NTHREADS; i++)
        { parmr[i].cmax  = 0;
          parmr[i].jmax  = 0;
          parmr[i].comax = 0;
          parmr[i].cover = NULL;
        }
      tfilt = 0;
    }

  nfilt = nhits = 0;

  if (alen == 0 || blen == 0)
    goto epilogue;

  { int    i, j, p;
    uint64 c;
    int    limit;

    MG_alist = asort;
    MG_blist = bsort;
    MG_comp  = comp;

    parmm[0].abeg = parmm[0].bbeg = 0;
    for (i = 1; i < NTHREADS; i++)
      { p = (int) ((((int64) alen) * i) >> NSHIFT);
        if (p > 0)
          { c = asort[p-1].code;
            while (asort[p].code == c)
              p += 1;
          }
        parmm[i].abeg = parmm[i-1].aend = p;
        parmm[i].bbeg = parmm[i-1].bend = find_tuple(asort[p].code,bsort,blen);
      }
    parmm[NTHREADS-1].aend = alen;
    parmm[NTHREADS-1].bend = blen;

    for (i = 0; i < NTHREADS; i++)
      for (j = 0; j < MAXGRAM; j++)
        parmm[i].hitgram[j] = 0;

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,count_thread,parmm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

    if (VERBOSE)
      printf("\n");
    if (MEM_LIMIT > 0)
      { int64 histo[MAXGRAM];
        int64 tom, avail;

        for (j = 0; j < MAXGRAM; j++)
          histo[j] = parmm[0].hitgram[j];
        for (i = 1; i < NTHREADS; i++)
          for (j = 0; j < MAXGRAM; j++)
            histo[j] += parmm[i].hitgram[j];

        avail = (int64) (MEM_LIMIT - (sizeof_DB(ablock) + sizeof_DB(bblock))) / sizeof(Double);
        if (avail > alen + 2*blen)
          avail = (avail - alen) / 2;
        else
          avail = avail - (alen + blen);
        avail *= .98;

        tom = 0;
        for (j = 0; j < MAXGRAM; j++)
          { tom += j*histo[j];
            if (tom > avail)
              break;
          }
        limit = j;

        if (limit <= 1)
          { fprintf(stderr,"\nError: Insufficient ");
            if (MEM_LIMIT == MEM_PHYSICAL)
              fprintf(stderr," physical memory (%.1fGb), reduce block size\n",
                             (1.*MEM_LIMIT)/0x40000000ll);
            else
              { fprintf(stderr," memory allocation (%.1fGb),",(1.*MEM_LIMIT)/0x40000000ll);
                fprintf(stderr," reduce block size or increase allocation\n");
              }
            fflush(stderr);
            exit (1);
          }
        if (limit < 10)
          { fprintf(stderr,"\nWarning: Sensitivity hampered by low ");
            if (MEM_LIMIT == MEM_PHYSICAL)
              fprintf(stderr," physical memory (%.1fGb), reduce block size\n",
                             (1.*MEM_LIMIT)/0x40000000ll);
            else
              { fprintf(stderr," memory allocation (%.1fGb),",(1.*MEM_LIMIT)/0x40000000ll);
                fprintf(stderr," reduce block size or increase allocation\n");
              }
            fflush(stderr);
          }
        if (VERBOSE)
          { printf("   Capping mutual k-mer matches over %d (effectively -t%d)\n",
                   limit,(int) sqrt(1.*limit));
            fflush(stdout);
          }

        for (i = 0; i < NTHREADS; i++)
          { parmm[i].nhits = 0;
            for (j = 1; j < limit; j++)
              parmm[i].nhits += j * parmm[i].hitgram[j];
            parmm[i].limit = limit;
          }
      }
    else
      for (i = 0; i < NTHREADS; i++)
        parmm[i].limit = INT32_MAX;

    nhits = parmm[0].nhits;
    for (i = 1; i < NTHREADS; i++)
      parmm[i].nhits = nhits += parmm[i].nhits;

    if (VERBOSE)
      { printf("   Hit count = ");
        Print_Number(nhits,0,stdout);
        if (nhits >= blen)
          printf("\n   Highwater of %.2fGb space\n",
                       (1. * (alen + 2*nhits)) / 67108864);
        else
          printf("\n   Highwater of %.2fGb space\n",
                       (1. * (alen + blen + nhits)) / 67108864);
        fflush(stdout);
      }

    if (nhits >= blen)
      bsort = (KmerPos *) Realloc(bsort,sizeof(SeedPair)*(nhits+1),
                                   "Reallocating damapper sort vectors");
    hhit = work1 = (SeedPair *) bsort;
    khit = work2 = (SeedPair *) Malloc(sizeof(SeedPair)*(nhits+1),
                                        "Allocating damapper hit vectors");
    if (hhit == NULL || khit == NULL || bsort == NULL)
      exit (1);

    MG_blist = bsort;
    MG_hits  = khit;

    for (i = NTHREADS-1; i > 0; i--)
      parmm[i].nhits = parmm[i-1].nhits;
    parmm[0].nhits = 0;

    for (i = 0; i < NTHREADS; i++)
      { parmm[i].kptr = parmx[i].tptr;
        for (p = 0; p < BPOWR; p++)
          parmm[i].kptr[p] = 0;
      }

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,merge_thread,parmm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

#ifdef TEST_PAIRS
    printf("\nSETUP SORT:\n");
    for (i = 0; i < HOW_MANY && i < nhits; i++)
      { SeedPair *c = khit+i;
        printf(" %5d / %5d / %5d /%5d\n",c->aread,c->bread,c->apos,c->apos-c->diag);
      }
#endif
  }

  { int   i;
    int64 x;

    x = 0;
    for (i = 0; i < NTHREADS-1; i++)
      { parmx[i].beg = x;
        parmx[i].end = x = parmm[i+1].nhits;
      }
    parmx[NTHREADS-1].beg = x;
    parmx[NTHREADS-1].end = nhits;

    khit = (SeedPair *) lex_sort(pairsort,(Double *) khit,(Double *) hhit,parmx);

    khit[nhits].aread = 0x7fffffff;
    khit[nhits].bread = 0x7fffffff;
    khit[nhits].diag  = 0x7fffffff;
    khit[nhits].apos  = 0;

#ifdef TEST_CSORT
    printf("\nCROSS SORT %lld:\n",nhits);
    for (i = 0; i < HOW_MANY && i <= nhits; i++)
      { SeedPair *c = khit+i;
        printf(" %5d / %5d / %5d /%5d\n",c->aread,c->bread,c->apos,c->apos-c->diag);
      }
#endif
  }

  { int    i, a;
    int64  p;

    MR_ablock = ablock;
    MR_bblock = bblock;
    MR_hits   = khit;
    MR_start  = start;

    parmr[0].hbeg = 0;
    parmr[0].abeg = 0;
    for (i = 1; i < NTHREADS; i++)
      { a = (i * ablock->nreads) >> NSHIFT;
        if (a > 0)
          p = find_aread(a,khit,nhits);
        else
          p = 0;
        parmr[i].hbeg = parmr[i-1].hend = p;
        parmr[i].abeg = parmr[i-1].aend = a;
      }
    parmr[NTHREADS-1].hend = nhits;
    parmr[NTHREADS-1].aend = ablock->nreads;

#ifdef NOTHREAD

    for (i = 0; i < NTHREADS; i++)
      chain_thread(parmr+i);

#else

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,chain_thread,parmr+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

#endif

    if (VERBOSE)
      { for (i = 0; i < NTHREADS; i++)
          nfilt += parmr[i].nfilt;
        tfilt += nfilt;
      }
  }

  free(work2);
  free(work1);

epilogue:
  if (VERBOSE)
    { int width;

      if (nhits <= 0)
        width = 1;
      else
        width = ((int) log10((double) nhits)) + 1;
      width += (width-1)/3;

      printf("\n     ");
      Print_Number(nhits,width,stdout);
      printf(" %d-mers (%e of matrix)\n     ",Kmer,(1.*nhits/atot)/btot);
      if (nfilt < 0)
        { Print_Number(-nfilt,width,stdout);
          printf(" candidates removed\n     ");;
        }
      else
        { Print_Number(nfilt,width,stdout);
          printf(" candidates added\n     ");;
        }
      Print_Number(tfilt,width,stdout);
      printf(" candidates (%e of matrix)\n",(1.*tfilt/atot)/btot);
      fflush(stdout);
    }
}

static char *NameBuffer(char *aname, char *bname)
{ static char *cat = NULL;
  static int   max = -1;
  int len;

  len = strlen(aname) + strlen(bname) + 100;
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      if ((cat = (char *) realloc(cat,max+1)) == NULL)
        { fprintf(stderr,"%s: Out of memory (Making path name)\n",Prog_Name);
          exit (1);
        }
    }
  return (cat);
}

void Reporter(char *aname, DAZZ_DB *ablock, char *bname, DAZZ_DB *bblock,
              Align_Spec *aspec, int mflag)
{
#ifndef NOTHREAD
  THREAD threads[NTHREADS];
#endif
  int   i;
  int   ncheck;
  char *fname;

  ncheck = 0;
  fname  = NameBuffer(aname,bname);

  MR_ablock = ablock;
  MR_bblock = bblock;
  MR_spec   = aspec;
  MR_doA    = ((mflag & FLAG_DOA) != 0);
  MR_doB    = ((mflag & FLAG_DOB) != 0);

  ablock->maxlen += SPACING;

  for (i = 0; i < NTHREADS; i++)
    { parmr[i].work  = New_Work_Data();
      if (mflag & 0x1)
        { sprintf(fname,"/tmp/%s.%s.M%d.las",aname,bname,i+1);
          parmr[i].afile = Fopen(fname,"w");
          if (parmr[i].afile == NULL)
            exit (1);
        }
      if (mflag & 0x2)
        { sprintf(fname,"/tmp/%s.%s.R%d.las",bname,aname,i+1);
          parmr[i].bfile = Fopen(fname,"w");
          if (parmr[i].bfile == NULL)
            exit (1);
        }
    }

#ifdef NOTHREAD

  for (i = 0; i < NTHREADS; i++)
    report_thread(parmr+i);

#else

  for (i = 0; i < NTHREADS; i++)
    pthread_create(threads+i,NULL,report_thread,parmr+i);

  for (i = 0; i < NTHREADS; i++)
    pthread_join(threads[i],NULL);

#endif

  for (i = 0; i < NTHREADS; i++)
    { Free_Work_Data(parmr[i].work);
      if (mflag & 0x1)
        fclose(parmr[i].afile);
      if (mflag & 0x2)
        fclose(parmr[i].bfile);
      if (VERBOSE)
        ncheck += parmr[i].ncheck;
    }

  if (VERBOSE)
    { printf("      ");
      Print_Number(ncheck,0,stdout);
      printf(" mapped segments\n");
    }

  if ( ! PROFILE) return;
   
  { int   i, a, size;
    int64 beg, cnt;
    FILE *afile, *dfile;

    afile = Fopen(Catenate(".",aname,".prof",".anno"),"w");
    dfile = Fopen(Catenate(".",aname,".prof",".data"),"w");

    size = sizeof(int64);
    fwrite(&(ablock->nreads),sizeof(int),1,afile);
    fwrite(&size,sizeof(int),1,afile);

    cnt = 0;
    for (i = 0; i < NTHREADS; i++)
      { beg = cnt;
        for (a = parmr[i].abeg; a < parmr[i].aend; a++)
          { fwrite(&cnt,sizeof(int64),1,afile);
            cnt += (ablock->reads[a].rlen-1)/SPACING + 2;
          }
        fwrite(parmr[i].cover,1,cnt-beg,dfile);
      }
    fwrite(&cnt,sizeof(int64),1,afile);
  }
}
