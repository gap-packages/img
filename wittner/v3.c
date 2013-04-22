#include <complex.h>
extern "C" {
#include <ppm.h>
}

int MAXITER = 1000;

double cnorm (_Complex double z)
{
  return creal(z)*creal(z) + cimag(z)*cimag(z);
}

int color (_Complex double a)
{
  _Complex double c = 2*a/(a+1), z[MAXITER];
  z[0] = c;
  double dz = 1.0;
  int i, j;

  for (i = 1; i < MAXITER; i++) {
    _Complex double zz = 1.0/(z[i-1]*z[i-1]);
    z[i] = 1.0 - (a+1.0)/z[i-1] + a*zz;
    dz *= cnorm(((a+1.0) - 2.0*a/z[i-1])*zz);

    if (cnorm(z[i]) < 1.e-10)
      return (MAXITER-i) % 3;
    if (cnorm(z[i]) > 1.e10)
      return (MAXITER+1-i) % 3;
    if (cnorm(z[i]-1.0) < 1.e-10)
      return (MAXITER+2-i) % 3;
    /* if (cnorm(z[i]-z[i/2]) < dz)
       return -2; */
    if (cnorm(z[i]-z[i/2]) < 1.e-10)
      return -1;
  }
  return -2; /* unknown color */
}

main(int argc, char *argv[]) {
  ppm_init (&argc, argv);

  if (argc != 6) {
    printf("Use: LL UR DX DY MAXITER\n");
    return -1;
  }

  _Complex double LL, UR;
  { double a, b; sscanf(argv[1], "%lf+%lfi", &a, &b); LL = a+1.0i*b; }
  { double a, b; sscanf(argv[2], "%lf+%lfi", &a, &b); UR = a+1.0i*b; }
  double DX = strtod (argv[3],NULL);
  double DY = strtod (argv[4],NULL);
  MAXITER = strtol (argv[5],NULL,10);

  int CSIZE = (creal(UR)-creal(LL))/DX + 1;
  int RSIZE = (cimag(UR)-cimag(LL))/DY + 1;

  pixel *array[RSIZE];
  for (int i = 0; i < RSIZE; array[i++] = ppm_allocrow(CSIZE));
  
  for (int i = 0; i < RSIZE; i++)
    for (int j = 0; j < CSIZE; j++) {
      _Complex double z = LL + j*DX + (RSIZE-i-1)*1.0i*DY;
      int c = color(z/(2.0-z));
      switch (c) {
      case -1:
	PPM_ASSIGN(array[i][j],2*PPM_MAXMAXVAL/3,2*PPM_MAXMAXVAL/3,2*PPM_MAXMAXVAL/3);
	break;
      case -2:
	PPM_ASSIGN(array[i][j],2*PPM_MAXMAXVAL/3,2*PPM_MAXMAXVAL/3,2*PPM_MAXMAXVAL/3);
	break;
      case 0:
	PPM_ASSIGN(array[i][j],PPM_MAXMAXVAL,PPM_MAXMAXVAL/2,PPM_MAXMAXVAL/2);
	break;
      case 1:
	PPM_ASSIGN(array[i][j],PPM_MAXMAXVAL/2,PPM_MAXMAXVAL,PPM_MAXMAXVAL/2);
	break;
      case 2:
	PPM_ASSIGN(array[i][j],PPM_MAXMAXVAL/2,PPM_MAXMAXVAL/2,PPM_MAXMAXVAL);
	break;
      }
    }

  {
    FILE *f = popen("convert - avs:-", "w");
    ppm_writeppm(f,array,CSIZE,RSIZE,PPM_MAXMAXVAL,0);
    pclose(f);
  }
  return 0;
}
