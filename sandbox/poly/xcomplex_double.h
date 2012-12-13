typedef double xreal;

class xcomplex { 
 public:
  std::complex<double> z;

  // constants
  static const long int MAX_EXP = DBL_MAX_EXP;
  static const long int MIN_EXP = DBL_MIN_EXP;
  static const double INFIN = DBL_MAX;
  static const double ZERO = 0.0;

  // constructor
  xcomplex () { };
//xcomplex ( const int a ){ z = std::complex<double>(a); };
  xcomplex ( const double a ) { z = std::complex<double>( a ); };
  xcomplex ( const double a, const double b ) { z = std::complex<double>( a, b ); };
  xcomplex ( const std::complex<double> newz ){ z = newz; };

  // operations
//xcomplex operator + ( ){ return( xcomplex( z ) ); };
  xcomplex operator - ( ){ return( xcomplex( -z ) ); };
  xcomplex operator + ( const xcomplex newz ) const { return( xcomplex( z + newz.z ) ); };
  xcomplex operator - ( const xcomplex newz ) const { return( xcomplex( z - newz.z ) ); };
  xcomplex operator * ( const xcomplex newz ) const { return( xcomplex( z * newz.z ) ); };
  xcomplex operator / ( const xcomplex newz ) const { return( xcomplex( z / newz.z ) ); };

  void operator += ( const xcomplex newz ){ z += newz.z; };
  void operator -= ( const xcomplex newz ){ z -= newz.z; };
  void operator *= ( const xcomplex newz ){ z *= newz.z; };
//void operator /= ( const xcomplex newz ){ z /= newz.z; };

//void operator = ( const xcomplex newz ) { z = newz.z; };
//void operator = ( const std::complex<double> newz ) { z = newz; };
//void operator = ( const double newz ) { z = (std::complex<double>)(newz); };
  
  unsigned int operator == ( const xcomplex newz ) const { return( z == newz.z ); };
  unsigned int operator != ( const xcomplex newz ) const { return( z != newz.z ); };

  friend double real ( const xcomplex );
  friend double imag ( const xcomplex );
  friend double abs ( const xcomplex );
  friend long int ilogbl( const xcomplex );
  friend std::ostream& operator << (std::ostream &s, xcomplex &newz);
  void xscalbln(long int e){ 
    z = std::complex<double>(scalbln( real( z ), e ), scalbln( imag( z ), e));
  }
}; 

xcomplex operator / ( const double a, const xcomplex newz ){ return( xcomplex( a / newz.z ) ); };
xcomplex operator * ( const double a, const xcomplex newz ){ return( xcomplex( a * newz.z ) ); };

std::ostream& operator << ( std::ostream &s, xcomplex &newz ){
  s << "(" << real( newz ) << "," << imag( newz ) << ")";
  return( s );
}

double real ( const xcomplex newz ){ return( real( newz.z ) ); };

double imag ( const xcomplex newz ){ return( imag(newz.z) ); };

double abs ( const xcomplex newz ){ return( abs( newz.z ) ); };

double norm ( const xcomplex newz ){ return( norm( newz.z ) ); };

long int ilogbl( const xcomplex newz ){ return( ilogbl( abs( newz ) ) ); };

double Precision( const xcomplex z ){ return( DBL_EPSILON ); };

long int PrecisionInt( const xcomplex z ){ return( ilogbl( Precision( z ) ) ); };
