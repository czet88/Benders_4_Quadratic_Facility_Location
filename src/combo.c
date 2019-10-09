
/* ======================================================================
      	     combo.c,    S.Martello, D.Pisinger, P.Toth     feb 1997
   ====================================================================== */

/* Revised version dec 2002, updated rudidiv() */

/* This is the COMBO algorithm described in 
 *
 *   S.Martello, D.Pisinger, P.Toth
 *   "Dynamic Programming and Strong Bounds for the 0-1 Knapsack Problem",
 *   submitted Management Science (1997)
 *
 * (c) copyright S.Martello, D.Pisinger, P.Toth.
 * The code may only be used for academic or non-commercial purposes.
 *
 * Further details on the project can also be found in
 *
 *   S.Martello, D.Pisinger, P.Toth
 *   "Dynamic programming and tight bounds for the 0-1 knapsack problem",
 *   Report 97/11, DIKU, University of Copenhagen
 *   Universitetsparken 1
 *   DK-2100 Copenhagen
 *
 * The code has been tested on a hp9000/735, and conforms with the
 * ANSI-C standard apart from some of the timing routines (which may
 * be removed).
 *  
 * Errors and questions are refered to:
 *   David Pisinger, associate professor
 *   DIKU, University of Copenhagen,
 *   Universitetsparken 1,
 *   DK-2100 Copenhagen.
 *   e-mail: pisinger@diku.dk
 *   fax: +45 35 32 14 01
 */

/* ======================================================================
                                  definitions
   ====================================================================== */

#define MINRUDI      1000    /* parameter M1 from paper: 1000 normally */
#define MINSET       2000    /* parameter M2 from paper: 2000 normally */
#define MINHEUR     10000    /* parameter M3 from paper: 10000 normally */
#define MAXSTATES 1500000   

#undef HASCHANCE             /* should strong upper bounds be used? */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "combo.h" /* mari */


/* ======================================================================
				   macros
   ====================================================================== */

#define srand(x)     srand48(x)
#define randm(x)    (lrand48() % (long) (x))

#define SYNC            5      /* when to switch to linear scan in bins */
#define SORTSTACK     200      /* stack for saving discarded intervals */
#define MINMED       1000      /* find exact median if larger size */
#define MAXV  (8*sizeof(long)) /* number of bits in a long integer */

#define TRUE  1                /* boolean value */
#define FALSE 0

#define LEFT  1                /* expansion of core in given direction */
#define RIGHT 2

#define PARTITION 1            /* should sort routine partition or sort */
#define SORTALL   2

#define DET(a1, a2, b1, b2)    ((a1) * (prod)(b2) - (a2) * (prod)(b1))
#define SWAP(a, b)   { register item q; q = *(a); *(a) = *(b); *(b) = q; }
#define NO(a,p)                ((int) ((p) - (a)->fitem + 1))
#define DIFF(a,b)              ((int) (((b)+1) - (a)))
/*#define TIME(t)                ((double) t / 1000)*/
#define MIN(a,b)               ((a) < (b) ? (a) : (b))
#define MAX(a,b)               ((a) > (b) ? (a) : (b))


/* ======================================================================
				 type declarations
   ====================================================================== */

//typedef int           boolean; /* logical variable         */
//typedef int           ntype;   /* number of states/items   */
//typedef long          itype;   /* item profits and weights */  
//typedef long          stype;   /* sum of profit or weight  */  
typedef unsigned long btype;   /* binary solution vector   */
typedef double        prod;    /* product of state, item   */

typedef int (*funcptr) (const void *, const void *);

/* item record */
/*typedef struct {*/
 /* itype   p;              /* profit                  */
 /* itype   w;              /* weight                  */
 /* boolean x;              /* solution variable       */
/*} item;*/   /*mari:comentat, ja definit a combo.h */

/* interval record */
typedef struct {
  item  *f;               /* first item in interval  */
  item  *l;               /* last item in interval   */
} interval;

/* state */
typedef struct {
  stype psum;             /* profit sum of state     */
  stype wsum;             /* weight sum of state     */
  btype vect;             /* corresponding (partial) solution vector */
} state;

/* set of partial vectors */
typedef struct pset {
  ntype size;             /* set size                */
  state *fset;            /* first element in set    */
  state *lset;            /* last element in set     */
  state *set1;            /* first element in array  */
  state *setm;            /* last element in array   */

  btype    vno;           /* current vector number   */
  item     *vitem[MAXV];  /* current last MAXV items */
  item     *ovitem[MAXV]; /* optimal set of items    */
  btype    ovect;         /* optimal solution vector */
} partset;


typedef struct { /* all info for solving separated problem */
  item     *fitem;        /* first item in problem */
  item     *litem;        /* last item in problem */
  item     *s;            /* core is [s,t] */
  item     *t;
  item     *b;            /* break item */
  item     *fpart;        /* first element in sorted core */
  item     *lpart;        /* last element in sorted core */
  stype    wfpart;        /* weight sum up to sorted core */
  item     *fsort;
  item     *lsort;
  stype    wfsort;
  stype    c;             /* capacity of problem */
  stype    z;             /* incumbent solution */
  stype    zwsum;         /* weight sum of incumbent solution */
  stype    lb;            /* lower bound */

  /* solutions may be represented in one of two ways: either a complete */
  /* array of items (fullsol = TRUE), or as last changes in dynamic     */
  /* programming enumeration (fullsol = FALSE) See description of partset */

  boolean  fullsol;       /* which representation of solution */
  item     *fsol;         /* first item in opt solution (fullsol=FALSE) */
  item     *lsol;         /* last item in opt solution (fullsol=FALSE) */
  item     *ffull;        /* start of item array (fullsol=TRUE) */
  item     *lfull;        /* end of item array (fullsol=TRUE) */
  partset  d;             /* set of states, including solution */

  stype    dantzig;       /* dantzig upper bound */
  stype    ub;            /* global upper bound */
  stype    psumb;         /* profit sum up to break item */
  stype    wsumb;         /* weight sum up to break item */

  stype    ps, ws, pt, wt;

  interval *intv1, *intv2;
  interval *intv1b, *intv2b;

  boolean relx;
  boolean master;
  int coresize;
} allinfo;


/* ======================================================================
                                debug variables
   ====================================================================== */

long simpreduced;
long iterates;
long maxstates;
long coresize;
long optsur;
long relaxations;
long relfeasible;
/* long reltime; */
long pitested;
long pireduced;
long dynheur;


/* ======================================================================
			      forward declarations
   ====================================================================== */

stype combo(item *f, item *l, stype c, stype lb, stype ub,
            boolean def, boolean relx);


/* ======================================================================
			          debug routines
   ====================================================================== */

static void cb_error(char *str, ...)
{
  va_list args;

  va_start(args, str);
  vprintf(str, args); printf("\n");
  va_end(args);
  printf("PROGRAM TERMINATED (combo)!!!\n\n");
  exit(-1);
}


/* ======================================================================
                                  timing routines
   ====================================================================== */

/* These timing routines are written for HP-UNIX, but should be portable
 * with minor changes. The timing function "times" should be portable, but
 * since it returns the number of clock-ticks used by the process, one must
 * convert them to seconds. The function "sysconf" returns the number of
 * clock-ticks per second when given the argument _SC_CLK_TCK, but apparently
 * there is no convention on what _SC_CLK_TCK should be (HP-UNIX is using
 * the value 2, while digital UNIX has value 3). The constant is defined
 * in <unistd.h>, but only if _POSIX_SOURCE is defined.
 */

/*#define _POSIX_SOURCE   */      /* to read <unistd.h> on digital UNIX */
/* #define _INCLUDE_POSIX_SOURCE *//* to read <unistd.h> on HP-UX */
/*#include <unistd.h>           /* define the constant _SC_CLK_TCK */
/*#include <sys/times.h>        /* timing routines */

/*
void give_time(long *time)
{ /* return the number of milliseconds used */
 /* struct tms timeend;
  double t1;
  times(&timeend);
  t1 = (double) (timeend.tms_utime) / sysconf(_SC_CLK_TCK);
  *time = t1 * 1000;
}*/


/* ======================================================================
				  palloc
   ====================================================================== */

static void pfree(void *p)
{ /* semi-own free routine which makes additional test */
  if (p == NULL) cb_error("freeing null");
  free(p);
}


static void *palloc(long size, long no)
{ /* semi-own alloc routine which makes additional test */
  long *p;
  size *= no;
  if (size == 0) size = 1;
  if (size != (size_t) size) cb_error("Alloc too big %ld", size);
  //printf("warning palloc \n");
  //p = malloc(size);
  //printf("warning palloc %ld \n", size);
  p = malloc(size);
  if(p == NULL) printf ("\nError: Memoria insuficiente\n");
  //printf("warning palloc2 \n");
  if (p == NULL) cb_error("no memory size %ld", size);
  return p;
}


/* ======================================================================
				push/pop
   ====================================================================== */

static void push(allinfo *a, int side, item *f, item *l)
{
  interval *pos;
  switch (side) {
    case LEFT : pos = a->intv1; (a->intv1)++; break;
    case RIGHT: pos = a->intv2; (a->intv2)--; break;
  }
  if (a->intv1 == a->intv2) cb_error("interval stack full");
  pos->f = f; pos->l = l;
}


static void pop(allinfo *a, int side, item **f, item **l)
{
  interval *pos;
  switch (side) {
    case LEFT : if (a->intv1 == a->intv1b) cb_error("pop left");
                (a->intv1)--; pos = a->intv1;
                *f = pos->f; *l = pos->l;
                break;
    case RIGHT: if (a->intv2 == a->intv2b) cb_error("pop right");
                (a->intv2)++; pos = a->intv2;
                *f = pos->f; *l = pos->l;
                break;
  }
}


/* ======================================================================
				improvesol
   ====================================================================== */

static void improvesol(allinfo *a, state *v)
{
//  long t;
  char *txt;

  if (v->wsum  > a->c) cb_error("wrong improvesol");
  if (v->psum <= a->z) cb_error("not improved solution");

  a->z     = v->psum;
  a->zwsum = v->wsum;
  a->fsol  = a->s;
  a->lsol  = a->t;
  a->d.ovect = v->vect;
  a->fullsol = FALSE;
  memcpy(a->d.ovitem, a->d.vitem, sizeof(item *) * MAXV);
}


/* ======================================================================
				definesolution
   ====================================================================== */

static void definesolution(allinfo *a)
{
  register item *i, *h, *m;
  register stype psum, wsum, ws;
  item *f, *l;
  btype j, k;
//  long t;

  /* endtime(&t); */

  /* full solutions are easy to update */
  if (a->fullsol) {
    for (i = a->fitem, m = a->litem+1, h = a->ffull; i != m; i++,h++) *i = *h;
    return;
  }

  /* partial solutions are more difficult as we only have MAXV last changes */

  /* initialize break solution */
  for (i = a->fitem, m = a->b; i != m; i++) i->x = 1;
  for (i = a->b, m = a->litem+1; i != m; i++) i->x = 0;

  /* backtrack vector */
  psum = a->z - a->psumb;
  wsum = a->zwsum - a->wsumb;
  f    = a->fsol;
  l    = a->lsol;

  /* backtrack */
  for (j = 0; j < MAXV; j++) {
    i = a->d.ovitem[j]; if (i == NULL) continue;
    k = a->d.ovect & ((btype) 1 << j);
    if (i->x == 1) {
      if (i > f) f = i;
      if (k) { psum += i->p; wsum += i->w; i->x = 0; }
    } else {
      if (i < l) l = i;
      if (k) { psum -= i->p; wsum -= i->w; i->x = 1; }
    }
  }
  if ((psum == 0) && (wsum == 0)) return;

  f++; l--; /* new core */
  psum = a->z; wsum = a->zwsum; iterates++;
  if (f > l) cb_error("wrong backtrack");
  for (i=a->fitem, m = f; i != m; i++) if (i->x) { psum-=i->p; wsum-=i->w; }
  for (i=l+1, m=a->litem+1; i != m; i++) if (i->x) { psum-=i->p; wsum-=i->w; }
  for (i = f, m = l+1, ws = 0; i != m; i++) ws += i->w; 
  if (ws == wsum) { 
    for (i = f, m = l+1; i != m; i++) i->x = 1;
  } else {
    combo(f, l, wsum, psum-1, psum, TRUE, TRUE);
  }
}


/* ======================================================================
			        rudidiv
   ====================================================================== */

static void rudidiv(allinfo *a)
{
  register item *i, *m, *b;
  register itype x, y, r;
  register prod pb, wb, q;
  register stype ws;

  b = a->b; pb = b->p; wb = b->w;
  q = DET(a->z+1-a->psumb, a->c-a->wsumb, pb, wb);
  x = a->fitem->w; ws = 0;
  for (i = a->fitem, m = a->litem+1; i != m; i++) {
    if ((i < b) && (DET(-i->p, -i->w, pb, wb) < q)) { ws += i->w; continue; }
    if ((i > b) && (DET(i->p, i->w, pb, wb) < q)) { continue; }
    y = x; x = i->w;
    while (y != 0) { r = x % y; x = y; y = r; }
    if (x == 1) return;
  }
  a->c = ws + x * ((a->c-ws) / x);
}


/* ======================================================================
				partsort
   ====================================================================== */

static void partsort(allinfo *a, item *f, item *l, stype ws, stype c, int what)
{
  register itype mp, mw;
  register item *i, *j, *m;
  register stype wi;
  int d;

  d = l - f + 1;
  if (d > 1) {
    m = f + d / 2;
    if (DET(f->p, f->w, m->p, m->w) < 0) SWAP(f, m);
    if (d > 2) {
      if (DET(m->p, m->w, l->p, l->w) < 0) {
        SWAP(m, l);
        if (DET(f->p, f->w, m->p, m->w) < 0) SWAP(f, m);
      }
    }
  }

  if (d > 3) {
    mp = m->p; mw = m->w; i = f; j = l; wi = ws;
    for (;;) {
      do { wi += i->w; i++; } while (DET(i->p, i->w, mp, mw) > 0);
      do {             j--; } while (DET(j->p, j->w, mp, mw) < 0);
      if (i > j) break;
      SWAP(i, j);
    }

    if (wi <= c) {
      if (what ==   SORTALL) partsort(a, f, i-1, ws, c, what);
      if (what == PARTITION) push(a, LEFT, f, i-1);
      partsort(a, i, l, wi, c, what);
    } else {
      if (what ==   SORTALL) partsort(a, i, l, wi, c, what);
      if (what == PARTITION) push(a, RIGHT, i,  l);
      partsort(a, f, i-1, ws, c, what);
    }
  }

  if ((d <= 3) || (what == SORTALL)) {
    a->fpart  = f;
    a->lpart  = l;
    a->wfpart = ws;
  }
}


/* ======================================================================
				minweights
   ====================================================================== */

static item *minweights(item *f, item *l, stype c)
{
  register itype mw;
  register item *i, *j, *m;
  register stype ws;
  int d;

  for (;;) {
    d = l - f + 1;
    if (d > 1) {
      m = f + d / 2;
      if (f->w > m->w) SWAP(f, m);
      if (d > 2) {
        if (m->w > l->w) { SWAP(m, l); if (f->w > m->w) SWAP(f, m); }
      }
    }
    if (d <= 3) break;
    mw = m->w; i = f; j = l; ws = 0;
    for (;;) {
      do { ws += i->w; i++; } while (i->w < mw);
      do {             j--; } while (j->w > mw);
      if (i > j) break;
      SWAP(i, j);
    }
    if (ws <= c) { f = i; c -= ws; } else l = i-1;
  }
  while (f->w <= c) { c -= f->w; f++; }
  return f;
}


/* ======================================================================
				maxprofits
   ====================================================================== */

static item *maxprofits(item *f, item *l, stype z)
{
  register itype mp;
  register item *i, *j, *m;
  register stype ps;
  int d;

  for (;;) {
    d = l - f + 1;
    if (d > 1) {
      m = f + d / 2;
      if (f->p < m->p) SWAP(f, m);
      if (d > 2) {
        if (m->p < l->p) { SWAP(m, l); if (f->p < m->p) SWAP(f, m); }
      }
    }
    if (d <= 3) break;
    mp = m->p; i = f; j = l; ps = 0;
    for (;;) {
      do { ps += i->p; i++; } while (i->p > mp);
      do {             j--; } while (j->p < mp);
      if (i > j) break;
      SWAP(i, j);
    }
    if (ps <= z) { f = i; z -= ps; } else l = i-1;
  }
  while (f->p <= z) { z -= f->p; f++; }
  return f;
}


/* ======================================================================
				sursort
   ====================================================================== */

static void sursort(item *f, item *l, itype sur, stype c,
             stype *p1, stype *w1, item **b)
{
  register itype s;
  register prod mp, mw;
  register item *i, *j, *m;
  register stype ws, ps;
  static item nn;
  item *l1;
  stype psum;
  int d;

  psum = 0; s = sur; l1 = l + 1;
  for (;;) {
    d = l - f + 1;
    if (d > 1) {
      m = f + d / 2;
      if (DET(f->p, f->w+s, m->p, m->w+s) < 0) SWAP(f, m);
      if (d > 2) {
        if (DET(m->p, m->w+s, l->p, l->w+s) < 0) {
          SWAP(m, l); if (DET(f->p, f->w+s, m->p, m->w+s) < 0) SWAP(f, m);
        }
      }
    }
    if (d <= 3) break;
    mp = m->p; mw = m->w+s; i = f; j = l; ws = ps = 0;
    for (;;) {
      do { ws+=i->w+s; ps+=i->p; i++; } while (DET(i->p,i->w+s,mp,mw) > 0);
      do {                       j--; } while (DET(j->p,j->w+s,mp,mw) < 0);
      if (i > j) break;
      SWAP(i, j);
    }
    if (ws <= c) { f = i; c -= ws; psum += ps; } else l = i-1;
  }
  for ( ; f != l1; f++) {
    if (f->w+s > c) { *p1 = psum; *w1 = c; *b  = f; return; }
    c -= f->w+s; psum += f->p;
  }
  nn.p = 0; nn.w = 1; *p1 = psum; *w1 = c; *b  = &nn;
}


/* ======================================================================
                                  haschance
   ====================================================================== */

static boolean haschance(allinfo *a, item *i, int side)
{
  register itype p, w;
//  register state *j, *m;
//  register stype pp, ww;

  if (a->d.size == 0) return FALSE;

#ifdef HASCHANCE
  if (side == RIGHT) {
    if (a->d.fset->wsum <= a->c - i->w) return TRUE;
    p = a->ps; w = a->ws; pitested++;
    pp = i->p - a->z - 1; ww = i->w - a->c;
    for (j = a->d.fset, m = a->d.lset + 1; j != m; j++) {
      if (DET(j->psum + pp, j->wsum + ww, p, w) >= 0) return TRUE;
    }
  } else {
    if (a->d.lset->wsum > a->c + i->w) return TRUE;
    p = a->pt; w = a->wt; pitested++;
    pp = -i->p - a->z - 1; ww = -i->w - a->c;
    for (j = a->d.lset, m = a->d.fset - 1; j != m; j--) {
      if (DET(j->psum + pp, j->wsum + ww, p, w) >= 0) return TRUE;
    }
  }
  pireduced++;
  return FALSE;
#else
  p = a->b->p; w = a->b->w;
  if (side == LEFT) {
    return (DET(a->psumb - i->p - a->z-1, a->wsumb - i->w - a->c, p, w) >= 0);
  } else {
    return (DET(a->psumb + i->p - a->z-1, a->wsumb + i->w - a->c, p, w) >= 0);
  }
#endif
}

/* ======================================================================
				  moveset
   ====================================================================== */

static void moveset(allinfo *a)
{
  register state *i, *j, *m;
  partset *d;

  /* move array to end if necessary */
  d = &a->d;
  if (d->lset != d->setm-1) {
    for (i = d->setm, j = d->lset, m = d->fset-1; j != m; j--) {
      i--; *i = *j;
    }
    d->fset = i; d->lset = d->setm-1;
  }
}


/* ======================================================================
				  multiply
   ====================================================================== */

static void multiply(allinfo *a, item *h, int side)
{
  register state *i, *j, *k, *m;
  register itype p, w;
  register btype mask0, mask1;
  state *r1, *rm;
  partset *d;

  d = &a->d; if (d->size == 0) return;
  if (side == RIGHT) { p = h->p; w = h->w; } else { p = -h->p; w = -h->w; }

  /* keep track on solution vector */
  d->vno++; if (d->vno == MAXV) d->vno = 0;
  mask1 = ((btype) 1 << d->vno); mask0 = ~mask1;
  d->vitem[d->vno] = h;

  /* initialize limits */
  r1 = d->fset; rm = d->lset; k = d->set1; m = rm + 1;
  k->psum = -1;
  k->wsum = r1->wsum + h->w + 1;
  m->wsum = rm->wsum + h->w + 1;

  for (i = r1, j = r1; (i != m) || (j != m); ) {
    if (i->wsum <= j->wsum + w) {
      if (i->psum > k->psum) {
        if (i->wsum > k->wsum) { k++; if ((k == i) || (k == j)) break; } 
        k->psum = i->psum; k->wsum = i->wsum;
        k->vect = i->vect & mask0;
      }
      i++;
    } else {
      if (j->psum + p > k->psum) {
        if (j->wsum + w > k->wsum) { k++; if ((k == i) || (k == j)) break; }
        k->psum = j->psum + p; k->wsum = j->wsum + w;
        k->vect = j->vect | mask1;
      }
      j++;
    }
  }
  if ((k == i) || (k == j)) cb_error("multiply, no space");

  d->fset = d->set1;
  d->lset = k;
  d->size = DIFF(d->fset,d->lset);

  if (d->size > maxstates) maxstates = d->size;
  a->coresize++;
  if (a->master) coresize++;
}


/* ======================================================================
			           surbin
   ====================================================================== */

static void surbin(item *f, item *l, itype s1, itype s2, stype c, 
                   stype dantzig, ntype card, stype *sur, stype *u)
{
  item *b;
  stype csur, r, psum, s, d, suropt;
  double ua, ub, gr, e, uopt;

  /* iterate surr. multiplier */
  uopt = dantzig; suropt = 0;
  for ( ; s1 <= s2; ) {
    s = (s2 + s1) / 2;
    csur = c + s * (long) card;
    if (csur < 0) csur = 0;
    sursort(f, l, s, csur, &psum, &r, &b);

    /* derive bound and gradient */
    e = 1; d = b-f;
    ua = psum + r * (prod) b->p / (b->w+s);
    ub = psum + (r + (card-d)*e) * (prod) b->p / (b->w+s+e);
    gr = (ub - ua) / e;

    if (ua < uopt) { suropt = s; uopt = ua; }
    if (gr > 0) s2 = s-1; else s1 = s+1;
  }
  *sur = suropt; *u = uopt;
}


/* ======================================================================
				solvesur
   ====================================================================== */

static void solvesur(allinfo *a, item *f, item *l, stype minsur, stype maxsur, 
              ntype card, stype *ub)
{
  register item *i, *k, *m;
  register stype ps, csur;
  register ntype no;
  stype sur, u, z1;
  boolean feasible;
 
  /* find optimal surrogate multiplier, and update ub */
  surbin(f, l, minsur, maxsur, a->c, a->dantzig, card, &sur, &u);
  optsur = sur;

  /* if bound <= current solution return */
  if ((u <= a->z) || (sur == 0)) {
    if (u > *ub) *ub = u;
    return;
  }
 
  /* add sur to weights and remove items with negative weight */
  csur = a->c + sur * (long) card; ps = 0;
  for (i = f, k = f, m = l+1; i != m; i++) {
    i->w += sur; i->x = 0;
    if (i->w <= 0) { ps += i->p; csur -= i->w; i->x = 1; SWAP(i, k); k++; }
  }

  /* solve problem to optimality */
  z1 = ps + combo(k, l, csur, a->z-ps, 0, TRUE, TRUE);

  /* subtract weight and check cardinality */
  for (i = f, m = l+1, no = 0; i != m; i++) { i->w -= sur; if (i->x) no++; }

  /* if feasible cardinality and improved, save solution in extra array */
  feasible = (no == card); if (feasible) relfeasible++;
  if ((z1 > a->z) && (feasible)) {
    for (i = f, k = a->ffull, m = l+1; i != m; i++, k++) *k = *i;
    a->z = z1; a->fullsol = TRUE;
  }

  /* output: maintain global upper bound */
  if (z1 > *ub) *ub = z1;
}


/* ======================================================================
				surrelax
   ====================================================================== */

static void surrelax(allinfo *a)
{
  register item *i, *j, *m;
  item *f, *l, *b;
  ntype n, card1, card2, b1;
  stype u, minsur, maxsur, wsum;
  itype minw, maxp, maxw;
//  long t1, t2;

  /* copy table */
/*  give_time(&t1);*/
  relaxations++;
  n = DIFF(a->fitem, a->litem);
  f = palloc(n, sizeof(item));
  l = f + n - 1;
  minw = a->fitem->w; maxp = maxw = wsum = 0;
  for (j = f, i = a->fitem, m = l+1; j != m; i++, j++) {
    *j = *i; wsum += i->w; 
    if (i->w < minw) minw = i->w;
    if (i->w > maxw) maxw = i->w;
    if (i->p > maxp) maxp = i->p;
  }

  /* find cardinality */
  b = a->b; b1 = DIFF(a->fitem, b-1);
  card1 = minweights(f, l, a->c) - f;      /* sum_{j=1}^{n} x_{j} \leq card1 */
  card2 = maxprofits(f, l, a->z) - f+1;    /* sum_{j=1}^{n} x_{j} \geq card2 */

  /* delimiters on sur.multipliers */
  maxsur = maxw;  /* should ideally be: maxp*maxp, but may cause overflow */
  minsur = -maxw; /* should ideally be: -maxp*maxw, but may cause overflow */

  /* choose strategy */
  u = 0;
  for (;;) {
    if (card2 == b1+1) { 
      solvesur(a, f, l, minsur, 0, b1+1, &u); /* min card constr */
      if (u < a->z) u = a->z; /* since bound for IMPROVED solution */
      break;
    }
    if (card1 == b1  ) { 
      solvesur(a, f, l, 0, maxsur, b1, &u); /* max card constr */
      break;
    }
    if (card1 == b1+1) { /* dichothomy: card <= b1 or card >= b1+1 */ 
      solvesur(a, f, l, minsur, 0, b1+1, &u); 
      solvesur(a, f, l, 0, maxsur, b1, &u); 
      break;
    }
    if (card2 == b1  ) { /* dichothomy: card <= b1 or card >= b1+1 */ 
      solvesur(a, f, l, 0, maxsur, b1, &u); 
      solvesur(a, f, l, minsur, 0, b1+1, &u); 
      break;
    }
    u = a->dantzig; break;
  }
  if (u < a->ub) a->ub = u;
  pfree(f);
/*  give_time(&t2); reltime = t2 - t1;*/
}


/* =========================================================================
				   simpreduce
   ========================================================================= */

static void simpreduce(int side, item **f, item **l, allinfo *a)
{
  register item *i, *j, *k;
  register prod pb, wb;
  register prod q;
  register int red;

  if (a->d.size == 0) { *f = *l+1; return; }
  if (*l < *f) return;

  pb = a->b->p; wb = a->b->w;
  q = DET(a->z+1-a->psumb, a->c-a->wsumb, pb, wb);
  i = *f; j = *l; red = 0;
  if (side == LEFT) {
    k = a->fsort - 1;
    while (i != j+1) {
      if (DET(-j->p, -j->w, pb, wb) < q) {
	SWAP(i, j); i++;       /* not feasible */
	red++;
      } else {
        SWAP(j, k); j--; k--;  /* feasible */
      }
    }
    *l = a->fsort - 1; *f = k + 1;
  } else {
    k = a->lsort + 1;
    while (i != j+1) {
      if (DET(i->p, i->w, pb, wb) < q) {
        SWAP(i, j); j--;       /* not feasible */
	red++;
      } else {
        SWAP(i, k); i++; k++;  /* feasible */
      }
    }
    *f = a->lsort + 1; *l = k - 1;
  }
  if (a->master) simpreduced += red;
}


/* ======================================================================
                                  findvect
   ====================================================================== */

static state *findvect(stype ws, state *f, state *l)
{
  /* find vector i, so that i->wsum <= ws < (i+1)->wsum */
  state *m;

  /* a set should always have at least one vector */
  /* if (f > l) cb_error("findvect: empty set"); */
  if (f->wsum >  ws) return NULL;
  if (l->wsum <= ws) return l;
  while (l - f > SYNC) {
    m = f + (l - f) / 2;
    if (m->wsum > ws) l = m-1; else f = m;
  }
  while (l->wsum > ws) l--;
  /* if (l->wsum     >  ws) cb_error("findvect: too big l"); */
  /* if ((l+1)->wsum <= ws) cb_error("findvect: too small l"); */
  return l;
}


/* ======================================================================
				  expandcore
   ====================================================================== */

static void expandcore(allinfo *a, boolean *atstart, boolean *atend)
{
  item *f, *l;

  /* expand core */
  *atstart = FALSE;
  if (a->s < a->fsort) {
    if (a->intv1 == a->intv1b) {
      *atstart = TRUE;
    } else {
      pop(a, LEFT, &f, &l); a->ps = f->p; a->ws = f->w;
      simpreduce(LEFT, &f, &l, a);
      if (f != l+1) {
	partsort(a, f, l, 0, 0, SORTALL); a->fsort = f;
	a->ps = a->s->p; a->ws = a->s->w;
      }
    }
  } else { a->ps = a->s->p; a->ws = a->s->w; }

  /* expand core */
  *atend = FALSE;
  if (a->t > a->lsort) {
    if (a->intv2 == a->intv2b) {
      *atend = TRUE;
    } else {
      pop(a, RIGHT, &f, &l); a->pt = l->p; a->wt = l->w;
      simpreduce(RIGHT, &f, &l, a);
      if (f != l+1) {
	partsort(a, f, l, 0, 0, SORTALL); a->lsort = l;
	a->pt = a->t->p; a->wt = a->t->w;
      }
    }
  } else { a->pt = a->t->p; a->wt = a->t->w; }
}


/* ======================================================================
				  reduceset
   ====================================================================== */

static void reduceset(allinfo *a)
{
  register state *i, *m, *k;
  register stype c, z;
  register prod p, w;
  state *v, *r1, *rm;
  boolean atstart, atend;

  if (a->d.size == 0) return;

  /* find break point and improve solution */
  r1 = a->d.fset; rm = a->d.lset;
  v = findvect(a->c, r1, rm);
  if (v == NULL) v = r1 - 1; else if (v->psum > a->z) improvesol(a, v);

  /* expand core, and choose ps, ws, pt, wt */
  expandcore(a, &atstart, &atend);

  /* now do the reduction */
  /* NB! This is the most efficient implementation, no product q is needed */
  c = a->c; z = a->z + 1; k = a->d.setm;
  if (!atstart) {
    p = a->ps; w = a->ws; 
    for (i = rm, m = v; i != m; i--) {
      if (DET(i->psum-z, i->wsum-c, p, w) >= 0) { k--; *k = *i; } 
    }
  }
  if (!atend) {
    p = a->pt; w = a->wt; 
    for (i = v, m = r1 - 1; i != m; i--) {
      if (DET(i->psum-z, i->wsum-c, p, w) >= 0) { k--; *k = *i; } 
    }
  }

  /* save limit */
  a->d.fset = k;
  a->d.lset = a->d.setm - 1; /* reserve one record for multiplication */
  a->d.size = DIFF(a->d.fset, a->d.lset);
}


/* ======================================================================
				  initfirst
   ====================================================================== */

static void initfirst(allinfo *a, stype pb, stype wb)
{
  partset *d;
  state *k;
  btype i;

  /* create table */
  d = &(a->d);
  d->set1 = palloc(sizeof(state), MAXSTATES);
  d->setm = d->set1 + MAXSTATES - 1;
  d->size = 1;
  d->fset = d->set1;
  d->lset = d->set1;

  /* init first state */
  k = d->fset;
  k->psum = pb;
  k->wsum = wb;
  k->vect = 0;

  /* init solution vector */
  for (i = 0; i < MAXV; i++) d->vitem[i] = NULL;
  d->vno = MAXV-1;

  /* init full solution */
  a->fullsol = FALSE;
  a->ffull = palloc(DIFF(a->fitem,a->litem), sizeof(item));
  a->lfull = a->ffull + DIFF(a->fitem,a->litem);
}


/* ======================================================================
                                 swapout
   ====================================================================== */

static void swapout(allinfo *a, item *i, int side)
{
  interval *pos, *k;

  if (side == LEFT) {
    for (pos = a->intv1b; pos != a->intv1; pos++) {
      if ((pos->f <= i) && (i <= pos->l)) {
        SWAP(i, pos->l); pos->l--; pos++; break;
      }
    }
    for (; pos != a->intv1; pos++) {
      SWAP(pos->f-1, pos->l); pos->f--; pos->l--;
    }
    /* remove empty intervals */
    for (pos = k = a->intv1b; pos != a->intv1; pos++) {
      if (pos->f <= pos->l) { *k = *pos; k++; }
    }
    a->intv1 = k;
  } else {
    for (pos = a->intv2b; pos != a->intv2; pos--) {
      if ((pos->f <= i) && (i <= pos->l)) {
        SWAP(i, pos->f); pos->f++; pos--; break;
      }
    }
    for (; pos != a->intv2; pos--) {
      SWAP(pos->f, pos->l+1); pos->f++; pos->l++;
    }
    /* remove empty intervals */
    for (pos = k = a->intv2b; pos != a->intv2; pos--) {
      if (pos->f <= pos->l) { *k = *pos; k--; }
    }
    a->intv2 = k;
  }
}


/* ======================================================================
                                 findcore
   ====================================================================== */

void findcore(allinfo *a)
{
  register item *i, *m;
  register itype p, r;
  item *j, *s, *t, *b;

  /* all items apart from b must be in intervals */
  s = t = b = a->b;
  if (a->fpart <= b-1) push(a, LEFT,  a->fpart, b-1);
  if (b+1 <= a->lpart) push(a, RIGHT, b+1, a->lpart);

  /* initial core is b-1, b, b+1 (if these exist) */
  if (b-1 >= a->fitem) { swapout(a, b-1, LEFT); s--; }
  if (b+1 <= a->litem) { swapout(a, b+1, RIGHT); t++; }

  /* forward greedy solution */
  if (b-1 >= a->fitem) {
    p = 0; r = a->c - a->wsumb + (b-1)->w;
    for (i = t+1, m = a->litem+1, j = NULL; i != m; i++) {
      if ((i->w <= r) && (i->p > p)) { p = i->p; j = i; }
    }
    if (j != NULL) { swapout(a, j, RIGHT); t++; }
  }

  /* second forward greedy solution */
  if (TRUE) {
    p = 0; r = a->c - a->wsumb;
    for (i = t+1, m = a->litem+1, j = NULL; i != m; i++) {
      if ((i->w <= r) && (i->p > p)) { p = i->p; j = i; }
    }
    if (j != NULL) { swapout(a, j, RIGHT); t++; }
  }

  /* backward greedy solution */
  if (TRUE) {
    j = NULL; r = a->wsumb - a->c + b->w;
    for (i = a->fitem, m = s; i != m; i++) if (i->w >= r) { p = i->p+1; break; }
    for (; i != m; i++) if ((i->w >= r) && (i->p < p)) { p = i->p; j = i; }
    if (j != NULL) { swapout(a, j, LEFT); s--; }
  }

  /* second backward solution */
  if (b+1 <= a->litem) {
    j = NULL; r = a->wsumb - a->c + b->w + (b+1)->w;
    for (i = a->fitem, m = s; i != m; i++) if (i->w >= r) { p = i->p+1; break; }
    for (; i != m; i++) if ((i->w >= r) && (i->p < p)) { p = i->p; j = i; }
    if (j != NULL) { swapout(a, j, LEFT); s--; }
  }

  /* add first and last item to ensure some variation in weights */
  if (a->fitem < s) { s--; swapout(a, a->fitem, LEFT); }
  if (a->litem > t) { t++; swapout(a, a->litem, RIGHT); }


  /* enumerate core: reductions are not allowed! */
  initfirst(a, a->psumb, a->wsumb); moveset(a);
  for (i = b, j = b-1; (i <= t) || (j >= s); ) {
    if (i <= t) { multiply(a, i,  RIGHT); moveset(a); i++; }
    if (j >= s) { multiply(a, j,  LEFT);  moveset(a); j--; }
  }
  a->s = s-1; a->fsort = s;
  a->t = t+1; a->lsort = t;
}


/* ======================================================================
				  heuristic
   ====================================================================== */

static void heuristic(allinfo *a)
{
  register item *i, *j, *m;
  register stype c, z, ub;
  register state *v, *r1, *rm;
  item *red, d;

  if (a->d.size == 0) return;

  /* define limits */
  dynheur++;
  r1 = a->d.fset; rm = a->d.lset;
  c = a->c; z = a->z; ub = a->ub;

  /* forward solution with dyn prog */
  if (a->intv2 != a->intv2b) {
    red = (a->intv2+1)->f;
    for (i = red, m = a->litem+1, j = NULL; i != m; i++) {
      v = findvect(c - i->w, r1, rm); if (v == NULL) continue;
      if (v->psum+i->p > z) { j = i; z = v->psum+i->p; if (z == ub) break; }
    }
    if (j != NULL) {
      swapout(a, j, RIGHT); d = *red;
      for (i = red, m = a->t; i != m; i--) *i = *(i-1);
      *(a->t) = d; a->lsort++;
      multiply(a, a->t, RIGHT); a->t++;
      reduceset(a);
    }
  }

  /* backward solution with dyn prog */
  if (a->intv1 != a->intv1b) {
    red = (a->intv1-1)->l;
    for (i = a->fitem, m = red+1, j = NULL; i != m; i++) {
      v = findvect(c + i->w, r1, rm); if (v == NULL) continue;
      if (v->psum-i->p > z) { j = i; z = v->psum-i->p; if (z == ub) break; }
    }
    if (j != NULL) {
      swapout(a, j, LEFT); d = *red;
      for (i = red, m = a->s; i != m; i++) *i = *(i+1);
      *(a->s) = d; a->fsort--;
      multiply(a, a->s, LEFT); a->s--;
      reduceset(a);
    }
  }
}


/* ======================================================================
				findbreak
   ====================================================================== */

static void findbreak(allinfo *a)
{
  register item *i;
  register stype psum, r;
  stype wsum;

  /* find break item */
  psum = 0; r = a->c;
  for (i = a->fitem; i->w <= r; i++) { psum += i->p; r -= i->w; }
  wsum = a->c - r;

  a->b       = i;
  a->psumb   = psum;
  a->wsumb   = wsum;
  a->dantzig = psum + (r * (prod) i->p) / i->w;
}


/* ======================================================================
				combo
   ====================================================================== */

extern stype combo(item *f, item *l, stype c, stype lb, stype ub,
            boolean def, boolean relx)
/* f,l : first, last item                                               */
/* c   : capacity of knapsack                                           */
/* lb  : lower bound. Solution vector is only updated if better z found */
/* ub  : upper bound. When upper bound is reached terminate immediately */
/* def : should solution vector be defined or just find objective value */
/* relx: relaxed problem is solved (no more relaxations will be made)   */
/* returns the objective value of the problem                           */
{
  allinfo a;
  interval *inttab;
  boolean heur, rudi;

  if ((ub != 0) && (lb == ub)) return lb;



  heur     = FALSE;
  rudi     = FALSE;
  //printf("warning combo13 \n");
  inttab   = palloc(sizeof(interval), SORTSTACK);
  // printf("warning combo14 \n");
  a.intv1b = &inttab[0];
  a.intv2b = &inttab[SORTSTACK - 1];
  a.intv1  = a.intv1b;
  a.intv2  = a.intv2b;
  a.fitem  = f;
  a.litem  = l;
  a.c      = c;
  a.z      = lb;
  a.lb     = lb;
  a.relx   = relx;
  a.master = (def && !relx);
  a.coresize = 0;
 
  partsort(&a, a.fitem, a.litem, 0, a.c, PARTITION);
  findbreak(&a);
  a.ub = (ub == 0 ? a.dantzig : ub);
 
  /* find and enumerate core */
  findcore(&a);
  reduceset(&a);

  while ((a.d.size > 0) && (a.z < a.ub)) {
    if (a.t <= a.lsort) {
      if (haschance(&a, a.t, RIGHT)) multiply(&a, a.t,  RIGHT);
       a.t++;
    }
    reduceset(&a);
    if (a.s >= a.fsort) {
      if (haschance(&a, a.s, LEFT)) multiply(&a, a.s,  LEFT);
      a.s--;
    }
    reduceset(&a);

    /* find better lower bound when needed */
    if ((!heur) && (a.d.size > MINHEUR)) { heuristic(&a); heur = TRUE; }

    /* find tight bound when needed */
    if ((!relx) && (a.d.size > MINSET)) { surrelax(&a); relx = TRUE; }

    /* use rudimentary divisibility to decrease c */
    if ((!rudi) && (a.d.size > MINRUDI)) { rudidiv(&a); rudi = TRUE; }
  }
  pfree(a.d.set1);
  pfree(inttab);

  if ((def) && (a.z > a.lb)) definesolution(&a);
  pfree(a.ffull);

  return a.z;
}


void pinta_mochila(item *mochila,stype C,int n)
{
	int i;

	printf("Maximizar  %dx0",mochila[0].p);
	for(i=1;i<n;i++){
		printf("+%dx%d",mochila[i].p,i);
	}
	printf("\n sujeto a  %dx0",mochila[0].w);
	for(i=1;i<n;i++){
		printf("+%dx%d",mochila[i].w,i);
	}
	printf("<= %d\n",(int)C);


	return;
}
