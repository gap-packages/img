#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <complex.h>
#include <stdio.h>

_Complex double gsl_get (const gsl_vector *vreal, int n)
{
  return gsl_vector_get (vreal, 2*n) + 1.0i * gsl_vector_get (vreal, 2*n+1);
}

void gsl_set (gsl_vector *vreal, int n, _Complex double x)
{
  gsl_vector_set (vreal, 2*n, creal(x));
  gsl_vector_set (vreal, 2*n+1, cimag(x));
}

void cprint (FILE *f, _Complex double z)
{
  fprintf(f,"%.10lf", creal(z));
  if (cimag(z)>=0.0)
    fprintf(f,"+%.10lf*I", cimag(z));
  else
    fprintf(f,"-%.10lf*I", -cimag(z));
}

_Complex double f(_Complex double a, _Complex double z)
{
  return 1.0 - (a+1.0)/z + a/(z*z);
}

_Complex double df(_Complex double a, _Complex double z)
{
  return ((a+1.0) - 2.0*a/z) / (z*z);
}

int parabolic (const gsl_vector *xreal, void *param, gsl_vector *freal)
{
  /* x = (a,c1,...,cn)
   * compute f = (prod f'_a(c_i)-multiplier,f_a(c_1)-c_2,...,f_a(c_n)-c_1
  */
  int n = xreal->size/2 - 1;
  double _Complex multiplier = *(_Complex double *)param;
  _Complex double a = gsl_get (xreal, 0), z[n];
  for (int i = 0; i < n; i++)
    z[i] = gsl_get (xreal, i+1);

  _Complex double m = 1.0;
  for (int i = 0; i < n; i++)
    m *= df(a,z[i]);

  gsl_set (freal, 0, m - multiplier);
  for (int i = 0; i < n; i++)
    gsl_set (freal, i+1, f(a,z[i]) - z[(i+1) % n]);

  return GSL_SUCCESS;
}

int main (int argc, char *argv[])
{
  if (argc != 6) {
    fprintf(stderr, "Use: %s a.real a.imag n mult.real mult.imag\n", argv[0]);
    return -1;
  }

  _Complex double a = strtod(argv[1],NULL)+1.0i*strtod(argv[2],NULL);
  int n = strtol(argv[3],NULL,10);
  _Complex double m0 = strtod(argv[4],NULL)+1.0i*strtod(argv[5],NULL);
  _Complex double z[100] = { 2.0*a/(a+1.0) };
  for (int i = 1; i < n; i++)
    z[i] = f(a,z[i-1]);

  const double precision = 1.e-10;
  double error;
  _Complex double multiplier = 0.0;
  double speed = 0.5; /* speed=1 means "instant jump" */

  do {
    _Complex double m1 = (1.0-speed)*multiplier + speed*m0;

    gsl_vector *x = gsl_vector_alloc (2*(1+n));
    gsl_set (x, 0, a);
    for (int i = 0; i < n; i++) gsl_set (x, i+1, z[i]);

    gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (gsl_multiroot_fsolver_hybrids, 2*(1+n));
    gsl_multiroot_function f = {&parabolic, 2*(1+n), &m1};
    gsl_multiroot_fsolver_set (s, &f, x);

    for (int iter = 0; ; iter++) {
      int status = gsl_multiroot_fsolver_iterate (s);

      if (iter == 100 || status) {
	speed *= 0.5; /* 50% slower */
	break;
      }
      if (gsl_multiroot_test_residual (s->f, precision) != GSL_CONTINUE) {
	multiplier = m1; /* go further */
	speed = (0.1+speed)/1.1; /* try 10% faster */
	break;
      }
    }

    a = gsl_get (s->x, 0);
    for (int i = 0; i < n; i++)
      z[i] = gsl_get (s->x, i+1);

    error = cabs(gsl_get (s->f, 0));

    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
  } while (cabs(multiplier - m0) > precision);

  cprint(stdout,a);
#if 1
  printf("\n");
#else
  printf(" error=%g\n",error);
#endif

  return 0;
}
