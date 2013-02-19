typedef long double xreal;

struct xcomplex { 
  std::complex<xreal> z;

  // constants
  static const long int MAX_EXP = LDBL_MAX_EXP;
  static const long int MIN_EXP = LDBL_MIN_EXP;
  static const xreal INFIN = LDBL_MAX;
  static const ZERO = 0.0;

  // constructor
  xcomplex () { };
  xcomplex ( const int a ){ z = std::complex<xreal>(a); };
  xcomplex ( const double a ) { z = std::complex<xreal>( a ); };                      // needed for constants 1.0
  xcomplex ( const xreal a ) { z = std::complex<xreal>( a ); };
  xcomplex ( const xreal a, const xreal b ) { z = std::complex<xreal>( a, b ); };
  xcomplex ( const std::complex<xreal> newz ){ z = newz; };

  // operations
  xcomplex operator + ( ){ return( xcomplex( z ) ); };
  xcomplex operator - ( ){ return( xcomplex( -z ) ); };

  void operator += ( const xcomplex newz ){ z += newz.z; };
  void operator -= ( const xcomplex newz ){ z -= newz.z; };
  void operator *= ( const xcomplex newz ){ z *= newz.z; };
  void operator /= ( const xcomplex newz ){ z /= newz.z; };

  void operator = ( const xcomplex newz ) { z = newz.z; };
  void operator = ( const std::complex<xreal> newz ) { z = newz; };
  void operator = ( const xreal newz ) { z = (std::complex<xreal>)(newz); };
  
  unsigned int operator == ( const xcomplex newz ){ return( z == newz.z ); };
  unsigned int operator != ( const xcomplex newz ){ return( z != newz.z ); };

  friend xreal real ( const xcomplex );
  friend xreal imag ( const xcomplex );
  friend xreal abs ( const xcomplex );
  friend long int ilogbl( const xcomplex );
  friend xcomplex scalbln( const xcomplex, long int );
  friend std::ostream& operator << (std::ostream &s, xcomplex &newz);
};

xcomplex operator + ( const xcomplex a, const xcomplex b ){ return( xcomplex( a.z + b.z ) ); };
xcomplex operator - ( const xcomplex a, const xcomplex b ){ return( xcomplex( a.z - b.z ) ); };
xcomplex operator * ( const xcomplex a, const xcomplex b ){ return( xcomplex( a.z * b.z ) ); };
xcomplex operator / ( const xcomplex a, const xcomplex b ){ return( xcomplex( a.z / b.z ) ); };
xcomplex operator * ( const xreal a, const xcomplex b ){ return( xcomplex( a * b.z ) ); };
xcomplex operator / ( const xreal a, const xcomplex b ){ return( xcomplex( a / b.z ) ); };

std::ostream& operator << ( std::ostream &s, xcomplex &newz ){
  s << "(" << real( newz ) << "," << imag( newz ) << ")";
  return( s );
}

xreal real ( const xcomplex newz ){ return( real( newz.z ) ); };

xreal imag ( const xcomplex newz ){ return( imag(newz.z) ); };

xreal abs ( const xcomplex newz ){ return( abs( newz.z ) ); };

long int ilogbl( const xcomplex newz ){ return( ilogbl( abs( newz ) ) ); };

xcomplex scalbln( const xcomplex z, long int e ){ 
  xcomplex b ( scalbln( real( z ), e ), scalbln( imag( z ), e ) );
  return( b );
}; 

long int PrecisionInt( const xcomplex z ){ return( ilogbl( LDBL_EPSILON ) ); };

xreal Precision( const xcomplex z ){ return( LDBL_EPSILON ); };
