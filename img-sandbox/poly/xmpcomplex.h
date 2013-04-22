#include <limits.h>
#include <mpfr.h>
#include <mpc.h>
#include <gmpxx.h>

#ifdef MPFR_REALS
typedef mpf_class xreal;
#else
typedef double xreal;
#endif

struct xcomplex {
  mpc_t z;

  static const mp_rnd_t default_rnd;
  static const int default_prec = 128;
  static const long int MAX_EXP;
  static const long int MIN_EXP;
  static const double INFIN;
  static const xcomplex ZERO;

  // constructor
  xcomplex(){ mpc_init2(z,default_prec); mpc_set_d_d(z,0.0,0.0,default_rnd); }
#ifdef MPFR_REALS
  xcomplex(const double a){ mpc_init2(z,default_prec); mpc_set_d(z,a,default_rnd); }
  xcomplex(const xreal r){ mpc_init2(z,default_prec); mpc_set_f(z,r.get_mpf_t(),default_rnd); };
  xcomplex(const xreal a,const xreal b){ mpc_init2(z,default_prec); mpc_set_f_f(z,a.get_mpf_t(),b.get_mpf_t(),default_rnd); };
#else
  xcomplex(const xreal r){ mpc_init2(z,default_prec); mpc_set_d(z,r,default_rnd); }
  xcomplex(const xreal a,const xreal b){ mpc_init2(z,default_prec); mpc_set_d_d(z,a,b,default_rnd); }
#endif
  xcomplex(const mpc_t newz){ mpc_init2(z,mpc_get_prec(newz)); mpc_set(z,newz,default_rnd); };
  ~xcomplex() { mpc_clear(z); }

  // operations
  xcomplex operator + () const { return(xcomplex(z)); };
  xcomplex operator - () const { xcomplex newz; mpc_neg(newz.z,z,default_rnd); return(newz); };

  void operator += (const xcomplex a){ mpc_add(z,z,a.z,default_rnd); };
  void operator -= (const xcomplex a){ mpc_sub(z,z,a.z,default_rnd); };
  void operator *= (const xcomplex a){ mpc_mul(z,z,a.z,default_rnd); };
  void operator /= (const xcomplex a){ mpc_div(z,z,a.z,default_rnd); };
  
  void operator = (const mpc_t newz){ mpc_init2(z,mpc_get_prec(newz)); mpc_set(z,newz,default_rnd) ;  };
  void operator = (const xcomplex newz){ mpc_set(z,newz.z,default_rnd); };

  unsigned int const operator == (const xcomplex newz) const { return(mpc_cmp(z,newz.z) == 0); };
  unsigned int operator != (const xcomplex newz) const { return(mpc_cmp(z,newz.z) != 0); };

  void cscalbln(long int);
};

const mp_rnd_t xcomplex::default_rnd = mpfr_get_default_rounding_mode();
const long int xcomplex::MAX_EXP = mpfr_get_emax();
const long int xcomplex::MIN_EXP = mpfr_get_emin();
const double xcomplex::INFIN = pow(2,xcomplex::MAX_EXP);
const xcomplex xcomplex::ZERO = xcomplex(0.0);

std::ostream& operator<<(std::ostream& os,const xcomplex& newz){
  return(os << mpc_get_str(10,mpc_get_prec(newz.z),newz.z,xcomplex::default_rnd));
  //  return(os << "(" << real(newz) << " " << imag(newz) << ")");
}

xcomplex operator + (const xcomplex a, const xcomplex b) {
  xcomplex newz; mpc_add(newz.z,a.z,b.z,xcomplex::default_rnd); return(newz);
}
xcomplex operator - (const xcomplex a, const xcomplex b) {
  xcomplex newz; mpc_sub(newz.z,a.z,b.z,xcomplex::default_rnd); return(newz);
}
xcomplex operator * (const xcomplex a, const xcomplex b) {
  xcomplex newz; mpc_mul(newz.z,a.z,b.z,xcomplex::default_rnd); return(newz);
}
xcomplex operator / (const xcomplex a, const xcomplex b) {
  xcomplex newz; mpc_div(newz.z,a.z,b.z,xcomplex::default_rnd); return(newz);
}

xreal abs(const xcomplex newz){ 
  xreal tmp;
  mpfr_t tmpfr;
  mpfr_init2(tmpfr,xcomplex::default_prec);
  mpc_abs(tmpfr,newz.z,xcomplex::default_rnd);
#ifdef MPFR_REALS
  mpfr_get_f(tmp.get_mpf_t(),tmpfr,xcomplex::default_rnd);
#else
  tmp = mpfr_get_d(tmpfr,xcomplex::default_rnd);
#endif
  mpfr_clear(tmpfr);
  return tmp;
}

xreal norm(const xcomplex newz){
  xreal tmp;
  mpfr_t tmpfr;
  mpfr_init2(tmpfr,xcomplex::default_prec);
  mpc_norm(tmpfr,newz.z,xcomplex::default_rnd);
#ifdef MPFR_REALS
  mpfr_get_f(tmp.get_mpf_t(),tmpfr,xcomplex::default_rnd);
#else
  tmp = mpfr_get_d(tmpfr,xcomplex::default_rnd);
#endif
  mpfr_clear(tmpfr);
  return tmp;
}

long int ilogbl(const xcomplex newz){
  long e = INT_MIN;
  if (mpfr_cmp_si(mpc_realref(newz.z),0) != 0)
    e = mpfr_get_exp(mpc_realref(newz.z));
  if (mpfr_cmp_si(mpc_imagref(newz.z),0) != 0)
    e = max(e,mpfr_get_exp(mpc_imagref(newz.z)));
  return e;
}

void xcomplex::cscalbln(long int a){
  if (a > 0) {
    mpfr_mul_2ui(mpc_realref(z),mpc_realref(z),a,xcomplex::default_rnd);
    mpfr_mul_2ui(mpc_imagref(z),mpc_imagref(z),a,xcomplex::default_rnd);
  } else {
    mpfr_mul_2ui(mpc_realref(z),mpc_realref(z),-a,xcomplex::default_rnd);
    mpfr_mul_2ui(mpc_imagref(z),mpc_imagref(z),-a,xcomplex::default_rnd);
  }
}

long int PrecisionInt(const xcomplex newz) { return(mpc_get_prec(newz.z)); }

xreal Precision(const xcomplex newz)
{
#ifdef MPFR_REALS
  xreal tmp = 1;
  mpf_mul_2exp(tmp.get_mpf_t(),tmp.get_mpf_t(),-PrecisionInt(newz));
  return tmp;
#else
  return scalbln(1.0,-PrecisionInt(newz));
#endif
}
