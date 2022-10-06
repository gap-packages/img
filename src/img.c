/****************************************************************************
 *
 * img.c                                                    Laurent Bartholdi
 *
 * Copyright (c) 2009-2012, Laurent Bartholdi
 *
 ****************************************************************************
 *
 * Call cpoly to compute the roots of a univariate polynomial
 * Call solver to normalize barycenter of points on S^2
 *
 ****************************************************************************/

#undef DEBUG_COMPLEX_ROOTS

#include "img.h"
#include <complex.h>

#ifdef MALLOC_HACK
#include <malloc.h>
#endif

/****************************************************************************
 * capture code that exits uncleanly rather than returning error message
 ****************************************************************************/
#ifdef CAPTURE_EXITS
static jmp_buf e_t_go_home;

static void baby_please_dont_go (void) {
  longjmp(e_t_go_home, 1);
}

/* in code:
...
#include <setjmp.h>
...
   atexit (baby_please_dont_go);
   if (setjmp(e_t_go_home))
     return __result;

   __result = fail;

   __result = call_bad_function();
   exit(0);
*/

#endif

/****************************************************************************
 * complex_roots of polynomial (as increasing-degree list of pairs (real,imag)
 ****************************************************************************/
#define cpoly cpoly_Cdouble
typedef double xreal;
typedef _Complex double xcomplex;
static const struct { Cdouble ZERO, INFIN; int MIN_EXP, MAX_EXP; }
  xdata = { 0.0, DBL_MAX, DBL_MIN_EXP, DBL_MAX_EXP };
static xreal xnorm(xcomplex z) { return __real__(z)*__real(z)+__imag__(z)*__imag__(z);}
static xreal xabs(xcomplex z) { return sqrt(xnorm(z)); }
static xreal xroot(xreal x, int n) { return pow(x,1.0/n); }
static int xlogb(xcomplex z) { return ilogb(xnorm(z)) / 2; }
#define xbits(z) DBL_MANT_DIG
#define xeta(z) DBL_EPSILON
typedef enum { false = 0, true = 1 } bool;
static void xscalbln (xcomplex *z, int e) {
  __real__(*z) = scalbln(__real__(*z), e);
  __imag__(*z) = scalbln(__imag__(*z), e);
}
#include "cpoly.C"

static Obj COMPLEX_ROOTS (Obj self, Obj coeffs)
{
  Obj result;
  Int i, numroots, degree = LEN_PLIST(coeffs)-1;
  xcomplex op[degree+1], zero[degree];
  bool real = true;
  
  if (degree < 1)
    return Fail;

  for (i = 0; i <= degree; i++) {
    __real__(op)[degree-i] = VAL_MACFLOAT(ELM_PLIST(ELM_PLIST(coeffs,i+1),1));
    __imag__(op)[degree-i] = VAL_MACFLOAT(ELM_PLIST(ELM_PLIST(coeffs,i+1),2));
    if (isnan(__real__(op)[degree-i]) || isnan(__imag__(op)[degree-i]))
      return Fail;
    if (__imag__(op)[degree-i] != 0.0)
      real = false;
  }

#ifdef DEBUG_COMPLEX_ROOTS
  fprintf(stderr,"coeffs");
  for (i = 0; i <= degree; i++)
    fprintf(stderr," %g+I*%g",(double)opr[i],(double)opi[i]);
   /* __asm__ __volatile__ ("int3"); */
  fprintf(stderr,"\n");
#endif

  numroots = cpoly (degree, op, zero);

  if (numroots == -1)
    return Fail;

#ifdef DEBUG_COMPLEX_ROOTS
  fprintf(stderr,"roots");
  for (i = 0; i < numroots; i++)
    fprintf(stderr," %g+I*%g",__real__(zero)[i],__imag__(zero)[i]);
  fprintf(stderr,"\n");
#endif

  result = ALLOC_PLIST(numroots);
  for (i = 1; i <= numroots; i++) {
    Obj t = ALLOC_PLIST(2);
    if (real && fabs(__imag__(zero)[i-1]) < 2*degree*DBL_EPSILON*fabs(__real__(zero)[i-1]))
      __imag__(zero)[i-1] = 0.0;
    set_elm_plist(t,1, NEW_MACFLOAT(__real__(zero)[i-1]));
    set_elm_plist(t,2, NEW_MACFLOAT(__imag__(zero)[i-1]));
    set_elm_plist(result,i, t);
  }
  return result;
}

/****************************************************************************
 * real_roots of polynomial (in increasing degree)
 ****************************************************************************/
static Obj REAL_ROOTS (Obj self, Obj coeffs)
{
  Obj result;
  Int i, numroots;
  int degree = LEN_PLIST(coeffs)-1;
  Cdouble opr[degree+1], zeror[degree], zeroi[degree];

  if (degree < 1)
    return Fail;

  for (i = 0; i <= degree; i++) {
    opr[degree-i] = VAL_MACFLOAT(ELM_PLIST(coeffs,i+1));
    if (isnan(opr[degree-i]))
      return Fail;
  }

  rpoly (opr, &degree, zeror, zeroi);
  numroots = degree;

  if (numroots < 0)
    return Fail;

  result = ALLOC_PLIST(numroots);
  for (i = 1; i <= numroots; i++) {
    if (zeroi[i-1] == 0.0)
      set_elm_plist(result,i, NEW_MACFLOAT(zeror[i-1]));
    else {
      Obj t = ALLOC_PLIST(2);
      set_elm_plist(t,1, NEW_MACFLOAT(zeror[i-1]));
      set_elm_plist(t,2, NEW_MACFLOAT(zeroi[i-1]));
      set_elm_plist(result,i, t);
    }
  }
  return result;
}

/****************************************************************************
 * NFFUNCTION reduces a list according to an IMG relation
 ****************************************************************************/
#define PUSH_LETTER(__v) {				\
    Int v = __v;					\
    Obj w = INTOBJ_INT(v);				\
    if (resulti && v == -INT_INTOBJ(ELM_PLIST(result,resulti)))	\
      resulti--;					\
    else {						\
      resulti++;					\
      if (resulti > allocn) {				\
        allocn *= 2;					\
	GROW_PLIST(result,allocn);			\
      }							\
      SET_ELM_PLIST(result,resulti,w);			\
    }							\
  }

#define MATCH_POS(p,__v) {					\
    Int v = __v;						\
    if (v > 0)							\
      p = INT_INTOBJ(ELM_PLIST(posind,v));			\
    else							\
      p = INT_INTOBJ(ELM_PLIST(negind,-v));			\
  }

Int Intabs(Int v) { return v > 0 ? v : -v; }

#define EXP_LETTER(v) (INT_INTOBJ(ELM_PLIST(exp,Intabs(v))))

#define FIND_MATCH(match,matchlen) {					\
  matchlen = 1;								\
  while (resulti > matchlen && ELM_PLIST(result,resulti-matchlen) == ELM_PLIST(rel,match+n-1)) { \
    matchlen++;								\
    if (!--match)							\
      match = n;							\
  }									\
}

static Obj NFFUNCTION(Obj self, Obj rel, Obj exp, Obj dir, Obj word)
{
  /* word is an integer lists. dir is true/false.
     rel is a list of lists: square of positive relator+square of negative
     relator; positions in 1st of letter i; position in 1st of letter -i
     exp is a list of exponents of the corresponding generators
   * if dir=true, replace all (>=1/2)-cyclic occurrences of rel in word by the shorter half
   * if dir=false, replace all occurrences of the last generator in word by the corresponding bit of rel
   */

  Obj posind = ELM_PLIST(rel,2), negind = ELM_PLIST(rel,3);
  rel = ELM_PLIST(rel,1);
  Int n = LEN_PLIST(posind), allocn = n, resulti = 0, match = 0, matchlen = 0, j, vlast = 0, power = 0, maxexp = 0;
  Obj result = ALLOC_PLIST(allocn);

  if (dir == False) { /* get rid of last generator */
    for (Int i = 1; i <= LEN_PLIST(word); i++) {
      Obj wi = ELM_PLIST(word,i);
      Int vi = INT_INTOBJ(wi);

      if (vi == n) {
	match = INT_INTOBJ(ELM_PLIST(negind,n));
	for (j = 1; j < n; j++)
	  PUSH_LETTER(INT_INTOBJ(ELM_PLIST(rel,j+match)));
      } else if (vi == -n) {
	match = INT_INTOBJ(ELM_PLIST(posind,n));
	for (j = 1; j < n; j++)
	  PUSH_LETTER(INT_INTOBJ(ELM_PLIST(rel,j+match)));
      } else
	PUSH_LETTER(vi);
    }
    SET_LEN_PLIST(result,resulti);
    return result;
  }

  for (Int i = 1; i <= LEN_PLIST(word); i++) {
    /* we produced result[1..resulti] as the compressed version of word[1..i].
       additionally, matchlen is maximal such that
       rel[match..match+matchlen-1] = result[resulti-matchlen+1..resulti]
       power is such that the last power letters in result are vlast^power.
    */
    Int vi = INT_INTOBJ(ELM_PLIST(word,i));
    if (EXP_LETTER(vi)==2)
      vi = Intabs(vi);
    Int idle;

    if (vi == -vlast) { /* pop letter */
      resulti--;
      matchlen--;
      power--;
    } else {
      PUSH_LETTER(vi);
      if (vi == vlast)
	power++;
      else
	power = 0; /* force recompute */
      if (match && vi == INT_INTOBJ(ELM_PLIST(rel,match+matchlen)))
	matchlen++;
      else
	matchlen = 0; /* force recompute */
    }

    do {
      idle = 1;
      if (maxexp && power >= (maxexp+1+(vlast > 0))/2) { /* apply power relation */
	Int delta = 2*power - maxexp;
	resulti -= delta;
	power -= delta;
	vlast = -vlast;
	for (Int j=0; j<power; j++)
	  SET_ELM_PLIST(result,resulti-j,INTOBJ_INT(vlast));
	matchlen = 0; /* force recompute */
	idle = 0;
      }

      if (matchlen >= (n+1+(match < 2*n))/2) { /* more than half a relation */
	resulti -= matchlen;
	for (j = n-1; j >= matchlen; j--) {
	  Int letter = INT_INTOBJ(ELM_PLIST(rel,j+match));
	  if (EXP_LETTER(letter)!=2) letter = -letter;
	  PUSH_LETTER(letter);
	}
	matchlen = n-matchlen;
	match = 4*n+1 - (match+n-1);
	power = 0; /* force recompute */
	idle = 0;
      }

      if (power == 0 && resulti) { /* recompute power */
	vlast = INT_INTOBJ(ELM_PLIST(result,resulti));
	maxexp = EXP_LETTER(vlast);
	power = 1;
	while (resulti > power && INT_INTOBJ(ELM_PLIST(result,resulti-power)) == vlast)
	  power++;
      }

      if (matchlen == 0 && resulti) { /* recompute match */
	Int last = INT_INTOBJ(ELM_PLIST(result,resulti)), match0, matchlen0;
	MATCH_POS(match,last);
	FIND_MATCH(match,matchlen);
	if (EXP_LETTER(last)==2) {
	  MATCH_POS(match0,-last);
	  FIND_MATCH(match0,matchlen0);
	  if (matchlen0 > matchlen) {
	    matchlen = matchlen0;
	    match = match0;
	  }
	}
	idle = 0;
      }
    } while (!idle);
  }

  SET_LEN_PLIST(result,resulti);
  return result;
}

/****************************************************************************
 * FIND_BARYCENTER finds a mobius transformation that centers points
 ****************************************************************************/
typedef struct {
  int n;
  Double (*points)[3];
} bparams;

#ifdef MALLOC_HACK
void *old_free_hook, *old_malloc_hook;

static void *
my_malloc_hook (size_t size, const void *caller)
{
  fprintf(stderr,"allocating %d\n", size);

  return *NewBag(T_DATOBJ, size + sizeof(Int));
}

static void
my_free_hook (void *ptr, const void *caller)
{
  fprintf(stderr,"freeing pointer %p\n", ptr);
}
#endif

static void barycenter (const double x[], double y[], int m, int n, const bparams *param)
{
  /* x is a "shifting" parameter; it is a vector in R^3, and
     describes the M\"obius transformation with north-south dynamics.
     more precisely, let t=|x|. in R^3, the transformation sends
     P to (2(1-t)P+(2-t+(v*P))v)/(1+(1-t)^2+(2-t)(v*P)).
     In particular, for t=0 it sends everything to v, and for t=1 it fixes P.

     The M\"obius transformation is 
  */
  int i, j;
  long double v[3]; /* a little extra precision */

  for (i = 0; i < 3; i++) v[i] = x[i];
  long double t = sqrtl(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  long double sum[3] = { 0.0, 0.0, 0.0 };

  for (j = 0; j < param->n; j++) {
    long double x[3], z = 0.0;

    for (i = 0; i < 3; i++) /* z = v*P */
      z += param->points[j][i] * v[i];

    long double d = 1.0 + (1.0-t)*(1.0-t) + (2.0 - t)*z;
    for (i = 0; i < 3; i++)
      x[i] = (2.0*(1.0-t)*param->points[j][i] + (2.0-t+z)*v[i]) / d;

    for (i = 0; i < 3; i++) sum[i] += x[i];
  }

  for (i = 0; i < 3; i++) y[i] = sum[i] / param->n;
}

/* given a set of points on S^2 \subset R^3, there is, up to rotations
   of the sphere, a unique M\"obius transformation that centers these
   points, i.e. such that their barycentre is (0,0,0). This follows
   from GIT, as Burt Totaro told me:

   Dear Laurent,

   Geometric invariant theory (GIT) gives a complete answer to your
   question. More concretely, the answer follows from the Kempf-Ness
   theorem in GIT, as I think Frances Kirwan first observed.

   Namely, given a sequence of N points p_1,...,p_N on the 2-sphere,
   there is a Mobius transformation that moves these points
   to have center of mass at the origin of R^3 if and only if either
   (1) fewer than N/2 of the points are equal to any given point
   in the 2-sphere; or
   (2) N is even, N/2 of the points are equal to one point
   in the 2-sphere, and the other N/2 points are equal
   to a different point in the 2-sphere.
   (In GIT terminology, condition (1) describes which
   N-tuples of points in S^2 = CP^1 are "stable",
   and (2) describes which N-tuples are "polystable"
   but not stable.) In particular, if p_1,...,p_N are all distinct
   and N is at least 2, then they can be centered
   by some Mobius transformation.
   
   This result was the beginning of many developments in GIT,
   such as Donaldson's notion of "balanced" metrics.
   Here is a good survey (where Theorem 4.13 is the statement
   above).

   R. P. Thomas. Notes on GIT and symplectic reduction for bundles
   and varieties. arXiv:math/0512411

   Burt Totaro

   This is also proven in <Cite Ref="MR2121737">.
*/

#include <levmar.h>

static Obj FIND_BARYCENTER (Obj self, Obj gap_points, Obj gap_iter)
{
#ifdef MALLOC_HACK
  old_malloc_hook = __malloc_hook;
  old_free_hook = __free_hook;
  __malloc_hook = my_malloc_hook;
  __free_hook = my_free_hook;
#endif

  UInt n = LEN_PLIST(gap_points);

  Double x[3];
  for (int j = 0; j < 3; j++)
    x[j] = 0.0;

  Double __points[n][3];
  bparams param = { n, __points };
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < 3; j++)
      param.points[i][j] = VAL_MACFLOAT(ELM_PLIST(ELM_PLIST(gap_points,i+1),j+1));

  int iter, max_iter = INT_INTOBJ(gap_iter);
  double info[LM_INFO_SZ];
  double opts[LM_OPTS_SZ] = { LM_INIT_MU, 1.e-33, 1.e-33, 0., LM_DIFF_DELTA };

  for (int i = 0; i < n; i++)
    for (int j = 0; j < 3; j++)
      x[j] += param.points[i][j];

  iter = dlevmar_dif ((void (*)(double*, double*, int, int, void*)) barycenter,
		      (double*) x, NULL, 3, 3, max_iter, opts, info, NULL, NULL, (void *) &param);

  Obj result = ALLOC_PLIST(3);
  Obj list = ALLOC_PLIST(3); set_elm_plist(result, 1, list);
  for (int i = 0; i < 3; i++)
    set_elm_plist(list, i+1, NEW_MACFLOAT(x[i]));
  list = ALLOC_PLIST(LM_INFO_SZ); set_elm_plist(result, 2, list);
  for (int i = 0; i < LM_INFO_SZ; i++)
    set_elm_plist(list, i+1, NEW_MACFLOAT(info[i]));
  set_elm_plist(result, 3, INTOBJ_INT(iter));

#ifdef MALLOC_HACK
  __malloc_hook = old_malloc_hook;
  __free_hook = old_free_hook;
#endif
  return result;
}

/****************************************************************************
 * other
 ****************************************************************************/

/****************************************************************************
 * interface to GAP
 ****************************************************************************/
static StructGVarFunc GVarFuncs [] = {
  { "COMPLEX_ROOTS_FR", 1, "coeffs", COMPLEX_ROOTS, "img.c:COMPLEX_ROOTS" },
  { "REAL_ROOTS_FR", 1, "coeffs", REAL_ROOTS, "img.c:REAL_ROOTS" },
  { "NFFUNCTION_FR", 4, "rel, exp, dir, word", NFFUNCTION, "img.c:NFFUNCTION" },
  { "FIND_BARYCENTER", 2, "points, iter", FIND_BARYCENTER, "img.c:FIND_BARYCENTER" },
  { 0 }
};

static Int InitKernel ( StructInitInfo * module )
{
  InitHdlrFuncsFromTable( GVarFuncs );
  InitP1Kernel();
  return 0;
}

/* 'InitLibrary' sets up gvars, rnams, functions */
static Int InitLibrary ( StructInitInfo * module )
{
  InitGVarFuncsFromTable( GVarFuncs );
  InitP1Library();
  return 0;
}

static StructInitInfo module = {
    .type = MODULE_DYNAMIC,
    .name = "img.c",
    .initKernel = InitKernel,
    .initLibrary = InitLibrary,
};

StructInitInfo * Init__Dynamic ( void )
{
 return &module;
}
/* img.c . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here */
