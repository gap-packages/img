class xcomplex { 
 public:
  std::complex<__float128> z;

  // constants
  static const long int MAX_EXP = LDBL_MAX_EXP;
  static const long int MIN_EXP = LDBL_MIN_EXP;
  static const long double INFIN = LDBL_MAX;

  // constructor
  xcomplex () { };
  xcomplex ( const int a ){ z = a; };
  xcomplex ( const long double a ) { z = a; };
  xcomplex ( const double a ) { z = a; };                      // needed for constants 1.0
  xcomplex ( const __float128 a, const __float128 b ) { z = std::complex<__float128>(a,b); };
  xcomplex ( std::complex<__float128> newz ){ z = newz; }

  // operations
  xcomplex operator + ( ){ return( xcomplex( z ) ); };
  xcomplex operator - ( ){ return( xcomplex( -z ) ); };
  xcomplex operator + ( const xcomplex newz ){ return( xcomplex( z + newz.z ) ); };
  xcomplex operator - ( const xcomplex newz ){ return( xcomplex( z - newz.z ) ); };
  xcomplex operator * ( const xcomplex newz ){ return( xcomplex( z * newz.z ) ); };
  xcomplex operator / ( const xcomplex newz ){ return( xcomplex( z / newz.z ) ); };

  void operator += ( const xcomplex newz ){ z += newz.z; };
  void operator -= ( const xcomplex newz ){ z -= newz.z; };
  void operator *= ( const xcomplex newz ){ z *= newz.z; };
  void operator /= ( const xcomplex newz ){ z /= newz.z; };

  void operator = ( const xcomplex newz ) { z = newz.z; };
  void operator = ( std::complex<__float128> newz ) { z = newz; };
  void operator = ( const long double newz ) { z = newz; };
  
  unsigned int operator == ( const xcomplex newz ){ return( z == newz.z ); };
  unsigned int operator != ( const xcomplex newz ){ return( z != newz.z ); };

  friend long double real ( const xcomplex );
  friend long double imag ( const xcomplex );
  friend long double abs ( const xcomplex );
  friend long int ilogbl( const xcomplex );
  friend xcomplex scalbln( const xcomplex, long int );
  friend std::ostream& operator << (std::ostream &s, xcomplex &newz);
};

xcomplex operator / ( const __float128 a, const xcomplex newz ){ return( xcomplex( a / newz.z ) ); };
xcomplex operator * ( const __float128 a, const xcomplex newz ){ return( xcomplex( a * newz.z ) ); };

std::ostream& operator << ( std::ostream &s, xcomplex &newz ){
  s << "(" << real( newz ) << "," << imag( newz ) << ")";
  return( s );
}

long double real ( const xcomplex newz ){ return( real( newz.z ) ); };

long double imag ( const xcomplex newz ){ return( imag(newz.z) ); };

long double abs ( const xcomplex newz ){ return( abs( std::complex<long double>(real(newz.z),imag(newz.z)) ) ); };

long int ilogbl( const xcomplex newz ){ return( ilogbl( abs( newz ) ) ); };

xcomplex scalbln( const xcomplex z, long int e ){ 
  xcomplex b ( scalbln( real( z ), e ), scalbln( imag( z ), e ) );
  return( b );
}; 

long int PrecisionInt( const xcomplex z ){ return( ilogbl( LDBL_EPSILON ) ); };

long double Precision( const xcomplex z ){ return( LDBL_EPSILON ); };
