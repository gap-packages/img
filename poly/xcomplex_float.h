class xcomplex { 
 public:
  std::complex<float> z;

  // constants
  static const long int MAX_EXP = FLT_MAX_EXP;
  static const long int MIN_EXP = FLT_MIN_EXP;
  static const float INFIN = FLT_MAX;

  // constructor
  xcomplex () { };
  xcomplex ( const int a ){ z = std::complex<float>(a); };
  xcomplex ( const float a ) { z = std::complex<float>( a ); };
  xcomplex ( const double a ) { z = std::complex<float>( a ); };                      // needed for constants 1.0
  xcomplex ( const float a, const float b ) { z = std::complex<float>( a, b ); };
  xcomplex ( const std::complex<float> newz ){ z = newz; };

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
  void operator = ( const std::complex<float> newz ) { z = newz; };
  void operator = ( const float newz ) { z = (std::complex<float>)(newz); };
  
  unsigned int operator == ( const xcomplex newz ){ return( z == newz.z ); };
  unsigned int operator != ( const xcomplex newz ){ return( z != newz.z ); };

  friend float real ( const xcomplex );
  friend float imag ( const xcomplex );
  friend float abs ( const xcomplex );
  friend long int ilogbl( const xcomplex );
  friend xcomplex scalbln( const xcomplex, long int );
  friend std::ostream& operator << (std::ostream &s, xcomplex &newz);
};

xcomplex operator / ( const float a, const xcomplex newz ){ return( xcomplex( a / newz.z ) ); };
xcomplex operator * ( const float a, const xcomplex newz ){ return( xcomplex( a * newz.z ) ); };

std::ostream& operator << ( std::ostream &s, xcomplex &newz ){
  s << "(" << real( newz ) << "," << imag( newz ) << ")";
  return( s );
}

float real ( const xcomplex newz ){ return( real( newz.z ) ); };

float imag ( const xcomplex newz ){ return( imag(newz.z) ); };

float abs ( const xcomplex newz ){ return( abs( newz.z ) ); };

long int ilogbl( const xcomplex newz ){ return( ilogbl( abs( newz ) ) ); };

xcomplex scalbln( const xcomplex z, long int e ){ 
  xcomplex b ( scalbln( real( z ), e ), scalbln( imag( z ), e ) );
  return( b );
}; 

long int PrecisionInt( const xcomplex z ){ return( ilogbl( FLT_EPSILON ) ); };

float Precision( const xcomplex z ){ return( FLT_EPSILON ); };
