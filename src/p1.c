/****************************************************************************
 *
 * p1.c                                                     Laurent Bartholdi
 *
 * Copyright (c) 2009-2012, Laurent Bartholdi
 *
 ****************************************************************************
 *
 * handle points on P1:
 * compute points to and from complex number; point antipode, barycentre;
 * rational maps, including Möbius transformations; convert to and from
 * rational map. Compose, invert, evaluate, compute preimages;
 * compute intersections of algebraic curves.
 ****************************************************************************/

#undef DEBUG_DELAUNAY
#undef DEBUG_COMPLEX_ROOTS

#include <math.h>
#include <complex.h>
#include "img.h"
#define IS_MACFLOAT(obj) (TNUM_OBJ(obj) == T_MACFLOAT)

typedef long double ldouble;
typedef _Complex long double ldcomplex;
typedef ldcomplex p1point;

Obj TYPE_P1POINT, TYPE_P1MAP, IsP1Point, IsP1Map,
  NewFloat, IsPMComplex, P1infinity, P1zero;

static Obj NEW_DATOBJ (size_t size, Obj type)
{
  Obj obj = NewBag(T_DATOBJ,sizeof(Obj)+size);
  SET_TYPE_DATOBJ(obj, type);
  return obj;
}
  
static void guarantee(Obj filter, const char *name, Obj obj)
{
  if (!IS_DATOBJ(obj) || DoFilter(filter, obj) != True) {
    ErrorMayQuit("FR: object must be a %s", (Int)name, 0);
  }
}

static int cisfinite(ldcomplex c)
{
#ifdef isfinite
  return isfinite(creall(c));
#else
  return __finitel(creall(c));
#endif
}

static ldouble cnorm (ldcomplex c)
{
  if (cisfinite(c))
    return c*~c;
  else
    return FLT_MAX;
}

static ldcomplex p1map_eval (int deg, ldcomplex *numer, ldcomplex *denom, p1point p);

/****************************************************************
 * points
 ****************************************************************/
static p1point GET_P1POINT(Obj obj) {
  guarantee(IsP1Point, "P1 point", obj);
  p1point p;
  memcpy (&p, ADDR_OBJ(obj)+1, sizeof p);
  return p;
}

static Obj NEW_P1POINT (p1point p)
{
  Obj obj = NEW_DATOBJ(sizeof p, TYPE_P1POINT);
  memcpy (ADDR_OBJ(obj)+1, &p, sizeof p);
  return obj;
}

static Obj NEW_P1POINT2 (ldcomplex n, ldcomplex d)
{
  if (d == 0.)
    return P1infinity;
  else
    return NEW_P1POINT(n/d);
}

static Obj NEW_COMPLEX (ldcomplex c) {
  Obj r = NEW_FLOAT(creal(c));
  Obj i = NEW_FLOAT(cimag(c));
  return CALL_3ARGS(NewFloat,IsPMComplex,r,i);
}

static Obj P1POINT2STRING(Obj self, Obj objprec, Obj obj)
{
  p1point q = GET_P1POINT(obj);
  Obj str = NEW_STRING(100);
  int len;
  int prec = INT_INTOBJ(objprec);

  if (cisfinite(q))
    len = sprintf(CSTR_STRING(str),"%.*Lg%+.*Lgi",prec,creall(q),prec,cimagl(q));
  else
    len = sprintf(CSTR_STRING(str),"P1infinity");

  SET_LEN_STRING(str, len);
  SHRINK_STRING(str);
  return str;
}

static Obj STRINGS2P1POINT(Obj self, Obj gapre, Obj gapim)
{
  while (!IsStringConv(gapre) || !IsStringConv(gapim)) {
    ErrorQuit("STRINGS2POINT: Expected 2 strings", 0L, 0L);
  }
  ldouble re, im;
  if (sscanf(CSTR_STRING(gapre),"%Lg", &re) != 1)
    return Fail;
  if (sscanf(CSTR_STRING(gapim),"%Lg", &im) != 1)
    return Fail;
  return NEW_P1POINT(re+im*1.0i);
}

static void p1point_c2 (ldcomplex *c, p1point p)
{
  if (!cisfinite(p))
    c[0] = 1.0, c[1] = 0.0;
  else if (cnorm(p) <= 1.0)
    c[0] = p, c[1] = 1.0;
  else
    c[0] = 1.0, c[1] = 1.0/p;
}

static Obj P1POINT2C2(Obj self, Obj obj)
{
  p1point q = GET_P1POINT(obj);
  Obj one = NEW_COMPLEX(1.);

  obj = ALLOC_PLIST(2);
  if (!cisfinite(q)) {
    set_elm_plist(obj,1, one);
    set_elm_plist(obj,2, NEW_COMPLEX(0.));
  } else if (cnorm(q) <= 1.0) {
    set_elm_plist(obj,1, NEW_COMPLEX(q));
    set_elm_plist(obj,2, one);
  } else {
    set_elm_plist(obj,1, one);
    set_elm_plist(obj,2, NEW_COMPLEX(1.0/q));
  }
  return obj;
}

static ldcomplex VAL_COMPLEX(Obj obj)
{
  return VAL_FLOAT(ELM_PLIST(obj,1)) + 1.0i*VAL_FLOAT(ELM_PLIST(obj,2));
}

static Obj C22P1POINT(Obj self, Obj num, Obj den)
{
  ldcomplex n = VAL_COMPLEX(num), d = VAL_COMPLEX(den);
  if (d == 0.)
    return NEW_P1POINT(1.0/0.0);
  else
    return NEW_P1POINT2(n,d);
}

static Obj EQ_P1POINT(Obj self, Obj p, Obj q)
{
  p1point pp = GET_P1POINT(p), pq = GET_P1POINT(q);
  if (!cisfinite(pp) || !cisfinite(pq))
    return cisfinite(pp) == cisfinite(pq) ? True : False;
  return pp == pq ? True : False;
}

static Obj LT_P1POINT(Obj self, Obj p, Obj q)
{
  p1point pp = GET_P1POINT(p), pq = GET_P1POINT(q);
  if (!cisfinite(pp)) return False;
  if (!cisfinite(pq)) return True;
  if (creall(pp) < creall(pq)) return True;
  if (creall(pp) > creall(pq)) return False;
  return cimagl(pp) < cimagl(pq) ? True : False;
}

p1point p1point_sphere (ldouble u[3])
{
  ldouble n = sqrtl(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  ldouble v[3] = { u[0]/n, u[1]/n, u[2]/n };

  if (v[2] > 0.0)
    return (v[0] + v[1]*1.0i) / (1.0 + v[2]);
  else if (v[0] == 0.0 && v[1] == 0.0)
    return 1.0/0.0;
  else
    return (1.0-v[2]) / (v[0] - v[1]*1.0i);  
}

static Obj P1Sphere(Obj self, Obj obj)
{
  if (!IS_PLIST(obj) || LEN_PLIST(obj)!=3 || !IS_MACFLOAT(ELM_PLIST(obj,1))
	 || !IS_MACFLOAT(ELM_PLIST(obj,2)) || !IS_MACFLOAT(ELM_PLIST(obj,3))) {
    ErrorMayQuit("FR: object must be a floatean list of length 3", 0, 0);
  }
  ldouble v[3];
  int i;
  for (i = 0; i < 3; i++)
    v[i] = VAL_FLOAT(ELM_PLIST(obj,i+1));
  return NEW_P1POINT(p1point_sphere(v));
}

void sphere_p1point (p1point p, ldouble s[3])
{
  if (cisfinite(p)) {
    ldouble n = cnorm(p);
    s[0] = 2.0*creall(p)/(1.0+n);
    s[1] = 2.0*cimagl(p)/(1.0+n);
    s[2] = (1.0-n)/(1.0+n);
  } else
    s[0] = s[1] = 0.0, s[2] = -1.0;
}

static Obj SphereP1(Obj self, Obj obj)
{
  ldouble s[3];
  int i;
  sphere_p1point (GET_P1POINT(obj), s);
  obj = ALLOC_PLIST(3);
  for (i = 0; i < 3; i++)
    set_elm_plist(obj,i+1, NEW_FLOAT((Double)s[i]));
  return obj;
}

static Obj SphereP1Y(Obj self, Obj obj)
{
  p1point p = GET_P1POINT(obj);
  if (cisfinite(p))
    return NEW_FLOAT((Double)2.0*cimagl(p)/(1.0+cnorm(p)));
  else
    return NEW_FLOAT(0.0);
}

static Obj P1Antipode(Obj self, Obj obj)
{
  p1point p = GET_P1POINT(obj);
  if (p == 0.)
    return P1infinity;
  else if (cisfinite(p))
    p = -~(1.0/p);
  else
    return P1zero;
  return NEW_P1POINT(p);
}

void clean_complex (ldcomplex *v, ldouble prec)
{
  if (cisfinite(*v)) {
    if (fabsl(cimagl(*v)) < prec*fabsl(creall(*v))) {
      ldouble z = creall(*v);
      if (fabsl(z-1.0) < prec)
	z = 1.0;
      if (fabsl(z+1.0) < prec)
	z = -1.0;
      *v = z;
    }
    if (fabsl(creall(*v)) < prec*fabsl(cimagl(*v))) {
      ldouble z = cimagl(*v);
      if (fabsl(z-1.0) < prec)
	z = 1.0;
      if (fabsl(z+1.0) < prec)
	z = -1.0;
      *v = z*1.0i;
    }
  }
}

static Obj CLEANEDP1POINT(Obj self, Obj objp, Obj objprec)
{
  ldcomplex p = GET_P1POINT(objp);
  ldouble prec = VAL_FLOAT(objprec);
  clean_complex(&p, prec);
  ldouble n = cnorm(p);
  if (n > 0.5/(prec*prec))
    return P1infinity;
  if (n < 2.0*prec*prec)
    p = 0.0;
  return NEW_P1POINT(p);
}

static Obj P1BARYCENTRE(Obj self, Obj list)
{
  int n = LEN_PLIST(list), i, j;
  ldouble s[3] = { 0.0, 0.0, 0.0 }, t[3];
  for (i = 0; i < n; i++) {
    sphere_p1point (GET_P1POINT(ELM_PLIST(list,i+1)), t);
    for (j = 0; j < 3; j++) s[j] += t[j];
  }
  for (j = 0; j < 3; j++) s[j] /= n;
  return NEW_P1POINT(p1point_sphere(s));
}

static Obj P1Midpoint(Obj self, Obj objp, Obj objq)
{
  p1point p = GET_P1POINT(objp), q = GET_P1POINT(objq);
  if (!cisfinite(p) && !cisfinite(q))
    return P1infinity;
  else if (!cisfinite(p)) {
    if (q == 0.0) return Fail;
    return NEW_P1POINT(q*(1.0+sqrtl(1.0+1.0/cnorm(q))));
  } else if (!cisfinite(q)) {
    if (p == 0.0) return Fail;
    return NEW_P1POINT(p*(1.0+sqrtl(1.0+1.0/cnorm(p))));
  } else {
    ldcomplex d = 1.0 + q*~p;
    if (d == 0.0) return Fail;
    ldouble a = sqrtl((1.0+cnorm(p))/(1.0+cnorm(q))*cnorm(d));
    return NEW_P1POINT2(a*q+d*p,a+d);
  }
}

static Obj P1Distance(Obj self, Obj objp, Obj objq)
{
  p1point p = GET_P1POINT(objp);
  p1point q = GET_P1POINT(objq);
  ldouble v;

  if (p == q || (!cisfinite(p) && !cisfinite(q)))
    return NEW_FLOAT(0.0);
  else if (!cisfinite(p)) {
    v = 1.0/cabsl(q);
  } else if (!cisfinite(q))
    v = 1.0/cabsl(p);
  else {
    ldcomplex d = 1.0+(~p)*q;
    if (d == 0.)
      v = 1.0/0.0;
    else
      v = cabsl((p-q)/d);
  }
  return NEW_FLOAT((Double)2.0*atan(v));
}

static Obj XRatio(Obj self, Obj p1, Obj p2, Obj p3, Obj p4)
{
  ldcomplex x[4][2];
  p1point_c2 (x[0], GET_P1POINT(p1));
  p1point_c2 (x[1], GET_P1POINT(p2));
  p1point_c2 (x[2], GET_P1POINT(p3));
  p1point_c2 (x[3], GET_P1POINT(p4));

  return NEW_COMPLEX ((x[0][0]*x[2][1]-x[2][0]*x[0][1])
		      / (x[1][0]*x[2][1]-x[2][0]*x[1][1])
		      * (x[1][0]*x[3][1]-x[3][0]*x[1][1])
		      / (x[0][0]*x[3][1]-x[3][0]*x[0][1]));
}

static Obj P1XRatio(Obj self, Obj p1, Obj p2, Obj p3, Obj p4)
{
  ldcomplex x[4][2];
  p1point_c2 (x[0], GET_P1POINT(p1));
  p1point_c2 (x[1], GET_P1POINT(p2));
  p1point_c2 (x[2], GET_P1POINT(p3));
  p1point_c2 (x[3], GET_P1POINT(p4));

  return NEW_P1POINT ((x[0][0]*x[2][1]-x[2][0]*x[0][1])
		      / (x[1][0]*x[2][1]-x[2][0]*x[1][1])
		      * (x[1][0]*x[3][1]-x[3][0]*x[1][1])
		      / (x[0][0]*x[3][1]-x[3][0]*x[0][1]));
}

static Obj P1Circumcentre(Obj self, Obj obja, Obj objb, Obj objc)
{
/* circumcentre := function(a,b,c)
 * local p;
 * p := Norm(a)*(b-c)+Norm(b)*(c-a)+Norm(c)*(a-b);
 * q := a*~b*(1+Norm(c))+b*~c*(1+Norm(a))+c*~a*(1+Norm(b));
 * centres := solve(p+(q-~q)*z+~p*z^2);
 * d := P1Distance(centres[1],a);
 * if d<pi/2 then
 *     return [centres[1],d];
 * else
 *     return [centres[2],pi-d];
 * fi;
 */
  ldcomplex p, q, v[3][2];
  int i;
  p1point_c2 (v[0], GET_P1POINT(obja));
  p1point_c2 (v[1], GET_P1POINT(objb));
  p1point_c2 (v[2], GET_P1POINT(objc));

  p = q = 0.0;
  for (i = 0; i < 3; i++) {
    ldcomplex *a = v[i], *b = v[(i+1)%3], *c = v[(i+2)%3];
    p += cnorm(a[0])*b[1]*c[1]*~(b[0]*c[1]-c[0]*b[1]);
    q += a[0]*~b[0]*~a[1]*b[1]*(cnorm(c[0])+cnorm(c[1]));
  }
  q = (q - ~q) / 2.0;

  ldcomplex centre;

  if (p == 0.0)
    centre = 0.0;
  else
    centre = (-q + csqrtl(q*q-cnorm(p))) / p;

  ldouble d = cabsl(centre*v[0][1] - v[0][0]) / cabsl(~centre*v[0][0]+v[0][1]);
  if (d > 1.0)
    d = 1.0/d, centre = -1.0/~centre;

  Obj result = ALLOC_PLIST(2);
  set_elm_plist(result,1, NEW_P1POINT(centre));
  set_elm_plist(result,2, NEW_FLOAT((Double)2.0*atan(d)));

  return result;
}

/****************************************************************
 * rational maps
 ****************************************************************/
static int p1map_degree(Obj obj) {
  guarantee (IsP1Map, "P1 map", obj);
  return (SIZE_OBJ(obj)-sizeof(Obj))/2/sizeof(ldcomplex)-1;
}

static ldcomplex *p1map_numer(Obj obj) {
  guarantee (IsP1Map, "P1 map", obj);
  return (ldcomplex *) (ADDR_OBJ(obj)+1);
}

static ldcomplex *p1map_denom(Obj obj) {
  guarantee (IsP1Map, "P1 map", obj);
  return (ldcomplex *) ((char *) (ADDR_OBJ(obj)+1) + (SIZE_OBJ(obj)-sizeof(Obj))/2);
}

static Obj NEW_P1MAP (int degree, ldcomplex *oldnumer, ldcomplex *olddenom)
{
  Obj obj = NEW_DATOBJ((degree+1)*2*sizeof(ldcomplex), TYPE_P1MAP);
  ldcomplex *numer = p1map_numer(obj), *denom = p1map_denom(obj);
  for (int i = 0; i <= degree; i++)
    numer[i] = oldnumer[i], denom[i] = olddenom[i];
  return obj;
}

static Obj NEW_SHRUNK_P1MAP (int degree, ldcomplex *numer, ldcomplex *denom)
{
  while (degree >= 0 && numer[degree] == 0. && denom[degree] == 0.) degree--;
  while (degree >= 0 && numer[0] == 0. && denom[0] == 0.) degree--, numer++, denom++;
  return NEW_P1MAP(degree, numer, denom);
}

static Obj MAT2P1MAP(Obj self, Obj obj)
{
  int deg = LEN_PLIST(ELM_PLIST(obj,1))-1, i, j;
  ldcomplex coeff[2][deg+1];
  for (i = 0; i < 2; i++)
    for (j = 0; j <= deg; j++)
      coeff[i][j] = VAL_COMPLEX(ELM_PLIST(ELM_PLIST(obj,i+1),j+1));
  return NEW_P1MAP(deg, coeff[0], coeff[1]);
}

static Obj P1MAP2MAT(Obj self, Obj map)
{
  Obj mat = ALLOC_PLIST(2);
  int deg = p1map_degree(map);
  Obj objnumer = ALLOC_PLIST(deg+1), objdenom = ALLOC_PLIST(deg+1);
  ldcomplex *numer = p1map_numer(map), *denom = p1map_denom(map);
  int i;
  for (i = 0; i <= deg; i++) {
    set_elm_plist(objnumer,i+1, NEW_COMPLEX(numer[i]));
    set_elm_plist(objdenom,i+1, NEW_COMPLEX(denom[i]));
  }
  set_elm_plist(mat,1, objnumer);
  set_elm_plist(mat,2, objdenom);
  return mat;
}

static Obj P1MAPDEGREE(Obj self, Obj obj)
{
  return INTOBJ_INT(p1map_degree(obj));
}

static Obj P1MAP3(Obj self, Obj objp, Obj objq, Obj objr)
{ /* Möbius transformation 0->p, 1->q, infty->r */
  ldcomplex p[2], q[2], r[2];
  p1point_c2 (p, GET_P1POINT(objp));
  p1point_c2 (q, GET_P1POINT(objq));
  p1point_c2 (r, GET_P1POINT(objr));
  ldcomplex pq = q[0]*p[1]-p[0]*q[1], qr = r[0]*q[1]-q[0]*r[1];
  ldcomplex numer[2] = { p[0]*qr, r[0]*pq }, denom[2] = { p[1]*qr, r[1]*pq };
  return NEW_P1MAP(1, numer, denom);
}

static Obj P1PATH(Obj self, Obj objp, Obj objq)
{ /* Möbius transformation 0->p, 1->q, infty->P1Antipode(p) */
  ldcomplex p[2], q[2], r[2];
  p1point_c2 (p, GET_P1POINT(objp));
  p1point_c2 (q, GET_P1POINT(objq));
  r[0] = -~p[1]; r[1] = ~p[0];
  ldcomplex pq = q[0]*p[1]-p[0]*q[1], qr = r[0]*q[1]-q[0]*r[1];
  ldcomplex numer[2] = { p[0]*qr, r[0]*pq }, denom[2] = { p[1]*qr, r[1]*pq };
  return NEW_P1MAP(1, numer, denom);
}

static Obj P1MAP2(Obj self, Obj objp, Obj objq)
{ /* Möbius transformation 0->p, infty->q */
  ldcomplex p[2], q[2];
  p1point_c2 (p, GET_P1POINT(objp));
  p1point_c2 (q, GET_P1POINT(objq));
  ldcomplex numer[2] = { p[0], q[0] }, denom[2] = { p[1], q[1] };
  return NEW_P1MAP(1, numer, denom);
}

static Obj P1MAP1(Obj self, Obj objp)
{ /* Möbius transformation infty->p */
  ldcomplex p[2];
  p1point_c2 (p, GET_P1POINT(objp));
  ldcomplex numer[2] = { -~p[1], p[0] }, denom[2] = { ~p[0], p[1] };
  return NEW_P1MAP(1, numer, denom);
}

static Obj P1ISPOLYNOMIAL(Obj self, Obj map)
{
  int deg = p1map_degree(map);
  ldcomplex *denom = p1map_denom(map);

  for (int i = 1; i <= deg; i++)
    if (denom[i] != 0.) return False;
  return True;
}

static void copy_poly (int deg, ldcomplex *coeff, int dega, ldcomplex *coeffa)
{ /* coeff = coeffa */
  int i;
  for (i = 0; i <= dega; i++)
    coeff[i] = coeffa[i];
  for (i = dega+1; i <= deg; i++)
    coeff[i] = 0.0;
}

static void der_poly (int deg, ldcomplex *coeff, int dega, ldcomplex *coeffa)
{ /* coeff = coeffa' */
  int i;
  for (i = 1; i <= dega; i++)
    coeff[i-1] = i*coeffa[i];
  for (i = dega; i <= deg; i++)
    coeff[i] = 0.0;
}

static void mul_poly (int deg, ldcomplex *coeff, int dega, ldcomplex *coeffa, int degb, ldcomplex *coeffb)
{ /* coeff = coeffa * coeffb */
  int i, j;
  for (i = 0; i <= deg; i++) {
    ldcomplex sum = 0.0;
    for (j = 0; j <= i; j++)
      if (j <= dega && i-j <= degb) sum += coeffa[j]*coeffb[i-j];
    coeff[i] = sum;
  }
}

static void mul_zminusb (int deg, ldcomplex *coeff, p1point b)
{
  coeff[deg+1] = 0.0;

  if (cisfinite(b)) {
    int i;
    for (i = deg+1; i >= 1; i--)
      coeff[i] = coeff[i-1] - b*coeff[i];
    coeff[0] *= -b;
  }
}

static void xplusequalay_poly (int deg, ldcomplex *coeff, ldcomplex k, int dega, ldcomplex *coeffa)
{ /* coeff += k*coeffa */
  int i;
  for (i = 0; i <= dega; i++)
    coeff[i] += k*coeffa[i];
}

static ldcomplex eval_poly (int deg, ldcomplex *coeff, ldcomplex x)
{ /* evaluate coeff at x */
  ldcomplex v = coeff[deg];
  int i;
  for (i = deg-1; i >= 0; i--)
    v = v*x + coeff[i];
  return v;
}

static ldcomplex eval_ylop (int deg, ldcomplex *coeff, ldcomplex x)
{ /* evaluate x^deg*coeff at 1/x */
  ldcomplex v = coeff[0];
  int i;
  for (i = 1; i <= deg; i++)
    v = v*x + coeff[i];
  return v;
}

static ldcomplex eval_d_poly (int deg, ldcomplex *coeff, ldcomplex x)
{ /* evaluate coeff' at x */
  ldcomplex v = deg*coeff[deg];
  int i;
  for (i = deg-1; i >= 1; i--)
    v = v*x + i*coeff[i];
  return v;
}

static ldcomplex eval_d_ylop (int deg, ldcomplex *coeff, ldcomplex x)
{ /* evaluate x^(deg-1)*coeff' at 1/x */
  ldcomplex v = coeff[1];
  int i;
  for (i = 2; i <= deg; i++)
    v = v*x + i*coeff[i];
  return v;
}

static Obj CLEANUPP1MAP(Obj self, Obj map, Obj objprec)
{
  int deg = p1map_degree(map), i, j;
  ldouble prec = VAL_FLOAT(objprec);
  ldcomplex coeff[2][deg+1];

  copy_poly (deg, coeff[0], deg, p1map_numer(map));
  copy_poly (deg, coeff[1], deg, p1map_denom(map));

  ldouble m;
  for (i = 0; i < 2; i++) {
    ldouble norm[deg+1];
    m = 0.0;
    for (j = 0; j <= deg; j++) {
      norm[j] = cnorm(coeff[i][j]);
      if (norm[j] > m) m = norm[j];
    }
    for (j = 0; j <= deg; j++)
      if (norm[j] < prec*m)
	coeff[i][j] = 0.0;
  }

  for (i = 0; i < deg && coeff[1][i] == 0.0; i++);
  ldcomplex c = coeff[1][i];
  for (i = 0; i < 2; i++)
    for (j = 0; j <= deg; j++) {
      coeff[i][j] /= c;
      clean_complex(&coeff[i][j], prec);
    }

  return NEW_P1MAP(deg, coeff[0], coeff[1]);
}

void compose_rat (int deg, ldcomplex *numer, ldcomplex *denom, int dega, ldcomplex *numera, ldcomplex *denoma, int degb, ldcomplex *numerb, ldcomplex *denomb)
{ /* compute numer/denom = numera/denoma @ numerb/denomb */
  int i, j;
  ldcomplex powb[dega+1][deg+1], temp[deg+1]; /* powb[i] will be numerb^i denomb^(dega-i) */
  copy_poly(deg, powb[0], degb, denomb);
  copy_poly(deg, powb[1], degb, numerb);
  for (i = 2; i <= dega; i++) {
    for (j = i; j > 0; j--)
      mul_poly (deg, powb[j], deg, powb[j-1], degb, numerb);
    mul_poly (deg, temp, deg, powb[0], degb, denomb);
    copy_poly (deg, powb[0], deg, temp);
  }
  copy_poly (deg, numer, -1, NULL);
  copy_poly (deg, denom, -1, NULL);
  for (i = 0; i <= dega; i++) {
    xplusequalay_poly (deg, numer, numera[i], deg, powb[i]);
    xplusequalay_poly (deg, denom, denoma[i], deg, powb[i]);
  }
}

static Obj COMPOSEP1MAP(Obj self, Obj a, Obj b)
{
  int dega = p1map_degree(a), degb = p1map_degree(b), deg = dega*degb;
  ldcomplex numer[deg+1], denom[deg+1];

  compose_rat (deg, numer, denom, dega, p1map_numer(a), p1map_denom(a), degb, p1map_numer(b), p1map_denom(b));
  return NEW_P1MAP(deg, numer, denom);
}

void invert_rat (ldcomplex *numer, ldcomplex *denom, ldcomplex *numera, ldcomplex *denoma)
{
  numer[0] = -numera[0];
  numer[1] = denoma[0];
  denom[0] = numera[1];
  denom[1] = -denoma[1];
}

static Obj INVERTP1MAP(Obj self, Obj map)
{
  if (p1map_degree(map) != 1)
    return Fail;
  ldcomplex numer[2], denom[2];
  invert_rat (numer, denom, p1map_numer(map), p1map_denom(map));
  return NEW_P1MAP(1, numer, denom);
}

static ldcomplex p1map_eval (int deg, ldcomplex *numer, ldcomplex *denom, p1point p)
{
  if (!cisfinite(p))
    return numer[deg] / denom[deg];
  if (cnorm(p) <= 1.0)
    return eval_poly (deg, numer, p) / eval_poly(deg, denom, p);
  p = 1.0 / p;
  return eval_ylop (deg, numer, p) / eval_ylop(deg, denom, p);
}

static Obj P1IMAGE(Obj self, Obj map, Obj objp)
{
  return NEW_P1POINT(p1map_eval(p1map_degree(map), p1map_numer(map), p1map_denom(map),
				GET_P1POINT(objp)));
}

#define cpoly cpoly_ldouble
typedef ldouble xreal;
typedef ldcomplex xcomplex;
static const struct { Cdouble ZERO, INFIN; int MIN_EXP, MAX_EXP; }
  xdata = { 0.0, LDBL_MAX, LDBL_MIN_EXP, LDBL_MAX_EXP };
static xreal xnorm(xcomplex z) { return __real__(z)*__real(z)+__imag__(z)*__imag__(z);}
static xreal xabs(xcomplex z) { return sqrtl(xnorm(z)); }
static xreal xroot(xreal x, int n) { return powl(x,1.0l/n); }
static int xlogb(xcomplex z) { return ilogbl(xnorm(z)) / 2; }
#define xbits(z) DBL_MANT_DIG
#define xeta(z) DBL_EPSILON
typedef enum { false = 0, true = 1 } bool;
static void xscalbln (xcomplex *z, int e) {
  __real__(*z) = scalblnl(__real__(*z), e);
  __imag__(*z) = scalblnl(__imag__(*z), e);
}
#include "cpoly.C"
static int roots_poly (int degree, ldcomplex *coeff, ldcomplex *zero)
{
  xcomplex op[degree+1];
  int i;
  for (i = 0; i <= degree; i++) /* high-degree coefficient first for cpoly */
    op[i] = coeff[degree-i];
  return cpoly_ldouble (degree, op, zero);
}

static int roots_rpoly (int degree, ldouble *coeff, ldcomplex *zero)
{
  long double opr[degree+1], zeror[degree], zeroi[degree];
  int i;
  while (degree > 0 && coeff[degree] == 0.0)
    degree--;
  for (i = 0; i <= degree; i++)
    opr[i] = coeff[degree-i];
  rpoly (opr, &degree, zeror, zeroi);
  for (i = 0; i < degree; i++)
    zero[i] = zeror[i] + 1.0i*zeroi[i];
  return degree;
}

static Obj P1PREIMAGE(Obj self, Obj map, Obj objp)
{
  while (p1map_degree(map) != 1)
    ErrorQuit("P1PREIMAGE: map must have degree 1, not %d", p1map_degree(map), 0);

  ldcomplex numer[2], denom[2];
  invert_rat (numer, denom, p1map_numer(map), p1map_denom(map));
  return NEW_P1POINT(p1map_eval(p1map_degree(map), p1map_numer(map), p1map_denom(map),
				GET_P1POINT(objp)));
}

static Obj P1PREIMAGES(Obj self, Obj map, Obj objp)
{
  ldcomplex p[2];
  p1point_c2 (p, GET_P1POINT(objp));
  int deg = p1map_degree(map), i;
  ldcomplex poly[deg+1], zero[deg];
  copy_poly (deg, poly, -1, NULL);
  xplusequalay_poly (deg, poly, p[1], deg, p1map_numer(map));
  xplusequalay_poly (deg, poly, -p[0], deg, p1map_denom(map));

  int numroots = roots_poly (deg, poly, zero);
  if (numroots < 0)
    return Fail;

  Obj obj = ALLOC_PLIST(deg);
  for (i = 0; i < numroots; i++) {
    // clean?
    set_elm_plist(obj,i+1, NEW_P1POINT(zero[i]));
  }
  for (i = numroots; i < deg; i++)
    set_elm_plist(obj,i+1, P1infinity);
  return obj;
}

static Obj P1CRITICAL(Obj self, Obj map)
{
  int deg = p1map_degree(map), i;
  ldcomplex poly[2*deg], temp[2*deg], der[deg], zero[2*deg-2];
  der_poly (deg-1, der, deg, p1map_numer(map));
  mul_poly (2*deg-1, poly, deg-1, der, deg, p1map_denom(map));
  der_poly (deg-1, der, deg, p1map_denom(map));
  mul_poly (2*deg-1, temp, deg-1, der, deg, p1map_numer(map));
  xplusequalay_poly (2*deg-2, poly, -1.0, 2*deg-2, temp);
  int numroots = roots_poly (2*deg-2, poly, zero);
  if (numroots < 0)
    return Fail;
  Obj obj = ALLOC_PLIST(2*deg-2);
  for (i = 0; i < numroots; i++) {
    // clean?
    set_elm_plist(obj,i+1, NEW_P1POINT(zero[i]));
  }
  for (i = numroots; i < 2*deg-2; i++)
    set_elm_plist(obj,i+1, P1infinity);
  return obj;
}

static Obj P1CONJUGATE(Obj self, Obj map)
{
  int deg = p1map_degree(map), i;
  ldcomplex num[deg+1], den[deg+1],
    *n = p1map_numer(map), *d = p1map_denom(map);

  for (i = 0; i <= deg; i++)
    num[i] = ~n[i];
  for (i = 0; i <= deg; i++)
    den[i] = ~d[i];

  return NEW_P1MAP(deg,num,den);
}

static Obj P1PRIMITIVE(Obj self, Obj map)
{
  int deg = p1map_degree(map), i, j;
  ldcomplex num[deg+2], den[deg+2], *n, *d;

  n = p1map_numer(map);
  d = p1map_denom(map);
  for (i = 0; i <= deg && d[i] == 0.; i++);
  for (j = deg; j >= i && d[j] == 0.; j--);
  if (i != j)
    return Fail;
  if (j == 0) { /* polynomial */
    num[0] = 0.;
    for (i = 0; i <= deg; i++)
      num[i+1] = n[i]/(i+1);
    for (i = 0; i <= deg+1; i++)
      den[i] = 0.;
    den[0] = d[0];
    return NEW_P1MAP(deg+1,num,den);
  } else { /* laurent polynomial */
   for (i = 0; i <= deg; i++) den[i] = 0.;
   den[j-1] = d[j];
   for (i = 0; i <= deg; i++) {
     if (i == j-1) {
       if (n[i] == 0.) num[i] = 0.; else return Fail;
     } else
       num[i] = n[i]/(i+1-j);
   }
   return NEW_SHRUNK_P1MAP(deg,num,den);
 }
}

static Obj P1DERIVATIVE(Obj self, Obj map)
{
  int deg = p1map_degree(map), i, j;
  ldcomplex *n, *d;

  if (deg == 0) {
    ldcomplex num[1], den[1];
    num[0] = 0.;
    den[0] = 1.;
    return NEW_P1MAP(0,num,den);
  }

  n = p1map_numer(map);
  d = p1map_denom(map);
  for (i = 0; i <= deg && d[i] == 0.; i++);
  for (j = deg; j >= i && d[j] == 0.; j--);
  if (j == 0) { /* polynomial */
    ldcomplex num[deg];
    for (i = 1; i <= deg; i++)
      num[i-1] = n[i]*i;
    return NEW_P1MAP(deg-1,num,d);
  } else if (i == j) { /* laurent polynomial */
    ldcomplex num[deg+2], den[deg+2];
    for (i = 0; i <= deg+1; i++) den[i] = 0.;
    den[j+1] = d[j];
    for (i = 0; i <= deg; i++)
      num[i] = n[i]*(i-j);
    num[deg+1] = 0.;
    return NEW_SHRUNK_P1MAP(deg+1,num,den);
  } else { /* rational */
    ldcomplex num[2*deg+1], den[2*deg+1], temp[2*deg+1], der[deg];
    der_poly (deg-1, der, deg, n);
    mul_poly (2*deg, num, deg-1, der, deg, d);
    der_poly (deg-1, der, deg, d);
    mul_poly (2*deg, temp, deg-1, der, deg, n);
    xplusequalay_poly (2*deg, num, -1.0, 2*deg, temp);
    mul_poly (2*deg, den, deg, d, deg, d);
    return NEW_SHRUNK_P1MAP(2*deg,num,den);
  }
}

static Obj P1MAPSCALING(Obj self, Obj map, Obj objp)
{
  p1point p = GET_P1POINT(objp);
  int deg = p1map_degree(map);
  ldcomplex *numer = p1map_numer(map), *denom = p1map_denom(map), diff;

  if (cisfinite(p)) {
    ldcomplex d = eval_poly(deg, denom, p);
    ldcomplex n = eval_poly(deg, numer, p);
    diff = eval_d_poly(deg, numer, p)*d - eval_d_poly(deg, denom, p)*n;
    diff *= (1.0 + cnorm(p));
    diff /= (cnorm(d) + cnorm(n));
  } else {
    diff = numer[deg]*denom[deg-1] - numer[deg-1]*denom[deg];
    diff /= (cnorm(numer[deg]) + cnorm(denom[deg]));
  }
  return NEW_FLOAT(cabsl(diff));
}

static Obj P1NUMERATOR(Obj self, Obj map)
{
  int deg = p1map_degree(map), i;
  ldcomplex den[deg+1];

  for (i = 0; i <= deg; i++)
    den[i] = 0.;
  den[0] = 1.;
  return NEW_SHRUNK_P1MAP(deg,p1map_numer(map),den);
}

static Obj P1DENOMINATOR(Obj self, Obj map)
{
  int deg = p1map_degree(map), i;
  ldcomplex den[deg+1];

  for (i = 0; i <= deg; i++)
    den[i] = 0.;
  den[0] = 1.;
  return NEW_SHRUNK_P1MAP(deg,p1map_denom(map),den);
}

/* can be optimized very much, O(n) instead of O(n^2), for polynomials... */
static Obj P1_SUM(Obj self, Obj a, Obj b)
{
  int dega = p1map_degree(a), degb = p1map_degree(b);
  ldcomplex num[dega+degb+1], den[dega+degb+1], temp[dega+degb+1];

  mul_poly (dega+degb, num, dega, p1map_numer(a), degb, p1map_denom(b));
  mul_poly (dega+degb, temp, dega, p1map_denom(a), degb, p1map_numer(b));
  xplusequalay_poly (dega+degb, num, 1., dega+degb, temp);
  mul_poly (dega+degb, den, dega, p1map_denom(a), degb, p1map_denom(b));
  return NEW_SHRUNK_P1MAP(dega+degb,num,den); 
}

static Obj P1_DIFF(Obj self, Obj a, Obj b)                                                
{
  int dega = p1map_degree(a), degb = p1map_degree(b);
  ldcomplex num[dega+degb+1], den[dega+degb+1], temp[dega+degb+1];

  mul_poly (dega+degb, num, dega, p1map_numer(a), degb, p1map_denom(b));
  mul_poly (dega+degb, temp, dega, p1map_denom(a), degb, p1map_numer(b));
  xplusequalay_poly (dega+degb, num, -1., dega+degb, temp);
  mul_poly (dega+degb, den, dega, p1map_denom(a), degb, p1map_denom(b));
  return NEW_SHRUNK_P1MAP(dega+degb,num,den);  
} 

static Obj P1_PROD(Obj self, Obj a, Obj b)                                                
{
  int dega = p1map_degree(a), degb = p1map_degree(b);
  ldcomplex num[dega+degb+1], den[dega+degb+1];

  mul_poly (dega+degb, num, dega, p1map_numer(a), degb, p1map_numer(b));
  mul_poly (dega+degb, den, dega, p1map_denom(a), degb, p1map_denom(b));
  return NEW_SHRUNK_P1MAP(dega+degb,num,den);  
} 

static Obj P1_QUO(Obj self, Obj a, Obj b)
{
  int dega = p1map_degree(a), degb = p1map_degree(b);
  ldcomplex num[dega+degb+1], den[dega+degb+1];

  mul_poly (dega+degb, num, dega, p1map_numer(a), degb, p1map_denom(b));
  mul_poly (dega+degb, den, dega, p1map_denom(a), degb, p1map_numer(b));
  return NEW_SHRUNK_P1MAP(dega+degb,num,den);
}

static Obj P1_AINV(Obj self, Obj map)
{
  int deg = p1map_degree(map), i;
  ldcomplex num[deg+1], *n;

  n = p1map_numer(map);
  for (i = 0; i <= deg; i++)
    num[i] = -n[i];
  return NEW_P1MAP(deg,num,p1map_denom(map));
}

static Obj P1_INV(Obj self, Obj map)
{
  int deg = p1map_degree(map);
  return NEW_P1MAP(deg,p1map_denom(map),p1map_numer(map));
}


static Obj P1MAPBYZEROSPOLES(Obj self, Obj gapzeros, Obj gappoles, Obj gapsrc, Obj gapdst)
{ /* construct a rational map with specified zeros and poles, and sending
     src to dst */
  p1point src = GET_P1POINT(gapsrc), dst = GET_P1POINT(gapdst);
  while (!IS_PLIST(gapzeros) || !IS_PLIST(gappoles) || LEN_PLIST(gapzeros) != LEN_PLIST(gappoles)) {
    ErrorQuit("P1MapByZerosPoles: first 2 arguments must be lists of same length", 0,0);
  }
  int deg = LEN_PLIST(gappoles), i;
  ldcomplex numer[deg+1], denom[deg+1];

  numer[0] = denom[0] = 1.0;
  for (i = 1; i <= deg; i++) {
    mul_zminusb (i-1, numer, GET_P1POINT(ELM_LIST(gapzeros,i)));
    mul_zminusb (i-1, denom, GET_P1POINT(ELM_LIST(gappoles,i)));
  }
  
  ldcomplex v;
  if (cisfinite(dst)) {
    v = dst / p1map_eval (deg, numer, denom, src);
  } else
    v = src;

  for (i = 0; i <= deg; i++)
    numer[i] *= v;

  return NEW_P1MAP(deg, numer, denom);
}

static Obj P1INTERSECT(Obj self, Obj gamma, Obj ratmap, Obj delta)
{ /* computes the (t,u) in [t0,1]x[0,1] such that gamma(t) = ratmap(delta(u)).
   * returns a list of [t,u,Im(gamma^-1*ratmap*delta)'(u),gamma(t),delta(u)]
   * gamma, delta are Möbius transformations, and ratmap is a rational map.
   */
  const ldouble eps = 1.0e-8;

  if (p1map_degree(gamma) != 1 || p1map_degree(delta) != 1)
    return Fail;

  int deg = p1map_degree(ratmap), i;
  ldcomplex numer[deg+1], denom[deg+1]; /* numer/denom = gamma^-1*ratmap*delta */
  {
    ldcomplex numer2[deg+1], denom2[deg+1];
    invert_rat (numer, denom, p1map_numer(gamma), p1map_denom(gamma));
    compose_rat (deg, numer2, denom2, 1, numer, denom, deg, p1map_numer(ratmap), p1map_denom(ratmap));
    compose_rat (deg, numer, denom, deg, numer2, denom2, 1, p1map_numer(delta), p1map_denom(delta));
  }

  ldcomplex poly[2*deg+1], zero[2*deg], conjdenom[deg+1]; /* poly = numer*~denom */
  for (i = 0; i <= deg; i++)
    conjdenom[i] = ~denom[i];
  mul_poly (2*deg, poly, deg, numer, deg, conjdenom);

  ldouble rpoly[2*deg+1]; /* rpoly = imag(poly) */
  for (i = 0; i <= 2*deg; i++)
    rpoly[i] = cimagl(poly[i]);

  int numroots = roots_rpoly (2*deg, rpoly, zero);
  if (numroots < 0)
    return Fail;

  Obj res = ALLOC_PLIST(0);
  for (i = 0; i < numroots; i++)
    if (cimagl(zero[i]) == 0.0 && creall(zero[i]) >= -eps && creall(zero[i]) <= 1.0+eps) {
      ldcomplex t = p1map_eval (deg, numer, denom, zero[i]); /* t = gamma^-1*ratmap*delta(u) */
      
      if (cisfinite(t) && cimagl(t) >= -1.0 && cimagl(t) <= 1.0 && creall(t) >= -eps && creall(t) <= 1.0+eps) { /* in fact, cimag(t) is microscopic; just avoid infinity */
	Obj tu = ALLOC_PLIST(5);
	set_elm_plist(tu,1, NEW_FLOAT(creal(t))); /* t */
	set_elm_plist(tu,2, NEW_FLOAT(creal(zero[i]))); /* u */
	ldcomplex d = eval_d_poly(2*deg, poly, zero[i]);
	ldouble y = cimagl(d) / cabsl(d); /* direction of approach */
	set_elm_plist(tu,3, INTOBJ_INT(y < -eps ? -1 : (y > eps ? 1 : 0)));
	set_elm_plist(tu,4, NEW_P1POINT(p1map_eval(1, p1map_numer(gamma), p1map_denom(gamma), t)));
	set_elm_plist(tu,5, NEW_P1POINT(p1map_eval(1, p1map_numer(delta), p1map_denom(delta), zero[i])));
	AddPlist(res,tu);
      }
    }
  return res;
}

static Obj P1ROTATION(Obj self, Obj points, Obj oldpoints)
{ /* find a Möbius transformation that matches points and oldpoints as well as possible.
   * points and oldpoints are normalized in that the last element is P1infinity.
   */
  int len = LEN_PLIST(points), i;

  ldcomplex proj[len];
  for (i = 0; i < len; i++) {
    ldcomplex p = GET_P1POINT(ELM_PLIST(points,i+1));
    if (cisfinite(p))
      proj[i] = 2.0*p / (1.0 + cnorm(p));
    else
      proj[i] = 0.0;
  }

  ldcomplex theta = 0.0;

  ldcomplex oldproj[len];
  for (i = 0; i < len; i++) {
    ldcomplex p = GET_P1POINT(ELM_PLIST(oldpoints,i+1));
    if (cisfinite(p))
      oldproj[i] = 2.0*p / (1.0 + cnorm(p));
    else
      oldproj[i] = 0.0;
  }
  ldouble n = 0.0;
  for (i = 0; i < len; i++)
    theta += ~proj[i]*oldproj[i], n += cnorm(proj[i]);
  theta /= n;
  
  if (n == 0.0 || cnorm(theta) < 0.7) /* no good rotation */
    theta = 0.0;

  if (theta == 0.0) { /* now just force the point of largest projection to be
			 on the positive real axis */
    ldouble q = 0.1;
    theta = 1.0;
    for (i = 0; i < len; i++) {
      ldouble n = cnorm(proj[i]);
      if (n > q)
	q = n, theta = ~proj[i];
    }
  }    
  theta /= cabsl(theta); /* make it of norm 1 */

  ldcomplex numer[2] = { 0., theta }, denom[2] = { 1., 0. };
  return NEW_P1MAP(1, numer, denom);
}

/* data to be passed to img */
static StructGVarFunc GVarFuncs[] = {
  { "P1POINT2C2", 1, "p1point", P1POINT2C2, "p1.c:P1POINT2C2" },
  { "C22P1POINT", 2, "num, den", C22P1POINT, "p1.c:C22P1POINT" },
  { "P1POINT2STRING", 2, "digits, p1point", P1POINT2STRING, "p1.c:P1POINT2STRING" },
  { "STRINGS2P1POINT", 2, "re, im", STRINGS2P1POINT, "p1.c:STRINGS2P1POINT" },
  { "EQ_P1POINT", 2, "p1point, p1point", EQ_P1POINT, "p1.c:EQ_P1POINT" },
  { "LT_P1POINT", 2, "p1point, p1point", LT_P1POINT, "p1.c:LT_P1POINT" },
  { "P1SPHERE", 1, "list", P1Sphere, "p1.c:P1Sphere" },
  { "SPHEREP1", 1, "p1point", SphereP1, "p1.c:SphereP1" },
  { "SPHEREP1Y", 1, "p1point", SphereP1Y, "p1.c:SphereP1Y" },
  { "CLEANEDP1POINT", 2, "p1point, prec", CLEANEDP1POINT, "p1.c:CLEANEDP1POINT" },
  { "P1ANTIPODE", 1, "p1point", P1Antipode, "p1.c:P1Antipode" },
  { "P1BARYCENTRE", 1, "list", P1BARYCENTRE, "p1.c:P1BARYCENTRE" },
  { "P1MIDPOINT", 2, "p1point, p1point", P1Midpoint, "p1.c:P1Midpoint" },
  { "P1DISTANCE", 2, "p1point, p1point", P1Distance, "p1.c:P1Distance" },
  { "XRATIO", 4, "p1point, p1point, p1point, p1point", XRatio, "p1.c:XRatio" },
  { "P1XRATIO", 4, "p1point, p1point, p1point, p1point", P1XRatio, "p1.c:P1XRatio" },
  { "P1CIRCUMCENTRE", 3, "p1point, p1point, p1point", P1Circumcentre, "p1.c:P1Circumcentre" },

  { "MAT2P1MAP", 1, "matrix", MAT2P1MAP, "p1.c:MAT2P1MAP" },
  { "P1MAP2MAT", 1, "p1map", P1MAP2MAT, "p1.c:P1MAP2MAT" },
  { "DEGREEOFP1MAP", 1, "p1map", P1MAPDEGREE, "p1.c:P1MAPDEGREE" },
  { "P1PATH", 2, "p1point, p1point", P1PATH, "p1.c:P1PATH" },
  { "P1MAP1", 1, "p1point", P1MAP1, "p1.c:P1MAP1" },
  { "P1MAP2", 2, "p1point, p1point", P1MAP2, "p1.c:P1MAP2" },
  { "P1MAP3", 3, "p1point, p1point, p1point", P1MAP3, "p1.c:P1MAP3" },
  { "P1ISPOLYNOMIAL", 1, "p1map", P1ISPOLYNOMIAL, "p1.c:P1ISPOLYNOMIAL" },
  { "CLEANEDP1MAP", 2, "p1map, float", CLEANUPP1MAP, "p1.c:CLEANUPP1MAP" },
  { "COMPOSEP1MAP", 2, "p1map, p1map", COMPOSEP1MAP, "p1.c:COMPOSEP1MAP" },
  { "INVERSEP1MAP", 1, "p1map", INVERTP1MAP, "p1.c:INVERTP1MAP" },
  { "P1IMAGE", 2, "p1map, p1point", P1IMAGE, "p1.c:P1IMAGE" },
  { "P1PREIMAGE", 2, "p1map, p1point", P1PREIMAGE, "p1.c:P1PREIMAGE" },
  { "P1PREIMAGES", 2, "p1map, p1point", P1PREIMAGES, "p1.c:P1PREIMAGES" },
  { "P1MAPCRITICALPOINTS", 1, "p1map", P1CRITICAL, "p1.c:P1CRITICAL" },
  { "P1MAPBYZEROSPOLES", 4, "zeros, poles, src, dst", P1MAPBYZEROSPOLES, "p1.c:P1MAPBYZEROSPOLES" },
  { "P1MAPCONJUGATE", 1, "p1map", P1CONJUGATE, "p1.c:P1CONJUGATE" },
  { "P1MAPPRIMITIVE", 1, "p1map", P1PRIMITIVE, "p1.c:P1PRIMITIVE" },
  { "P1MAPDERIVATIVE", 1, "p1map", P1DERIVATIVE, "p1.c:P1DERIVATIVE" },
  { "P1MAPSCALING", 2, "p1map, p1point", P1MAPSCALING, "p1.c:P1MAPSCALING" },
  { "P1MAPNUMERATOR", 1, "p1map", P1NUMERATOR, "p1.c:P1NUMERATOR" },
  { "P1MAPDENOMINATOR", 1, "p1map", P1DENOMINATOR, "p1.c:P1DENOMINATOR" },
  { "P1MAP_SUM", 2, "p1map,p1map", P1_SUM, "p1.c:P1_SUM" },
  { "P1MAP_DIFF", 2, "p1map,p1map", P1_DIFF, "p1.c:P1_DIFF" },
  { "P1MAP_PROD", 2, "p1map,p1map", P1_PROD, "p1.c:P1_PROD" },
  { "P1MAP_QUO", 2, "p1map,p1map", P1_QUO, "p1.c:P1_QUO" },
  { "P1MAP_INV", 1, "p1map", P1_INV, "p1.c:P1_INV" },
  { "P1MAP_AINV", 1, "p1map", P1_AINV, "p1.c:P1_AINV" },
  { "P1INTERSECT_IEEE754", 3, "p1map, p1map, p1map", P1INTERSECT, "p1.c:P1INTERSECT" },
  { "P1ROTATION_IEEE754", 2, "p1points, p1points/float", P1ROTATION, "p1.c:P1ROTATION" },
  { 0 } };

void InitP1Kernel(void)
{
  InitHdlrFuncsFromTable (GVarFuncs);

  ImportGVarFromLibrary ("TYPE_IEEE754P1POINT", &TYPE_P1POINT);
  ImportGVarFromLibrary ("TYPE_IEEE754P1MAP", &TYPE_P1MAP);  
  ImportGVarFromLibrary ("IsIEEE754P1Point", &IsP1Point);
  ImportGVarFromLibrary ("IsIEEE754P1Map", &IsP1Map);
  ImportGVarFromLibrary ("NewFloat", &NewFloat);  
  ImportGVarFromLibrary ("IsPMComplex", &IsPMComplex);  
  ImportGVarFromLibrary ("P1infinity", &P1infinity);
  ImportGVarFromLibrary ("P1zero", &P1zero);
}

void InitP1Library(void)
{
  InitGVarFuncsFromTable (GVarFuncs);
}

/* p1.c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here */
