



//----------------------basicNumber----------------------------------------






template <typename TScalar>
 const basicNumber<TScalar> basicNumber<TScalar>::Zero = basicNumber<TScalar>((TScalar)0);

template <typename TScalar>
 const basicNumber<TScalar> basicNumber<TScalar>::One = basicNumber<TScalar>((TScalar)1);









//-------------------------Constructors--------------------------------------------------------

template <typename TScalar>
inline   basicNumber<TScalar>::basicNumber()  : x(0), eps(0) {};


/// constructs a basic number {s + 0*EPS }
template <typename TScalar>
inline basicNumber<TScalar>::basicNumber (TScalar s) : x(s), eps(0)	 {	 };


/// constructs a basic number, epsPrecision is ignored and assumed =1
template <typename TScalar>
inline basicNumber<TScalar>::basicNumber (TScalar epsPrecision, std::string s ) :x(0), eps(0)
{
	if (epsPrecision>1)
	{
 		std::cerr << " allowed epsPrecision is { 0, 1 }" << std::endl;
		throw(" allowed epsPrecision is { 0, 1 }");
	}
	//	assert(epsPrecision<2);
};

/// copy constructor
template <typename TScalar>
inline basicNumber<TScalar>::basicNumber (const  basicNumber<TScalar> & z) 	: x(z.x), eps(z.eps)	 {};


	
/// constructs a basic number {s + t*eps }
template <typename TScalar>
inline basicNumber<TScalar>::basicNumber (TScalar s, TScalar t) 	: x(s), eps(t)	 	 {};


template <typename TScalar>
 inline bool  basicNumber<TScalar>::wellDefined(unsigned int characteristic)
{
	//std::cerr << "wellDefined: characteristic" << characteristic << std::endl;

	assert(sizeof(long long)>sizeof(TScalar));	//otherwise the following test does not work

	long long 	test  = pow( 2.0, (int)sizeof(TScalar)*8 )-1;
	TScalar test2 = test;
	test = test2;

	if (test<0)
	{
		
		//std::cerr << "TScalar is signed!" << std::endl;
		assert(pow( 2.0, (int)sizeof(TScalar)*8-1 )>characteristic);
	}
	else
	{
		std::cerr << "warning: TScalar is unsigned!" << std::endl;
		assert(pow( 2.0, (int)sizeof(TScalar)*8 )>characteristic);
	}
	return true;
}




//-------------------------Properties--------------------------------------------------------
/// returns maximum possible epsPrecision (here: 1)
template <typename TScalar>
inline unsigned short basicNumber<TScalar>::getEpsPrecision() const
{
	return 1;
}


template <typename TScalar>
inline short basicNumber<TScalar>::getEpsDegree() const
{
	if ( getEps() !=0 ) 
		return 1;
	else if ( getX() !=0 ) 
		return 0;
	else return -1;
}

//-----------------------------------data Access-------------------------------------------------------

template <typename TScalar>
inline   TScalar   basicNumber<TScalar>::getX() const 
{
	return x;
}

template <typename TScalar>
inline  TScalar   basicNumber<TScalar>::getEps() const
{
	return (eps);
}

template <typename TScalar>
inline TScalar& basicNumber<TScalar>::operator[](int i)
{
	#ifdef SAFE
		assert(i==0||i==1);
	#endif

	if (i==0)
		return x;
	if (i==1)
		return eps;
	
}




/// set coeff of e^epsPrecision to 'coeff'. TODO: 
template <typename TScalar>
inline  void   basicNumber<TScalar>::setValue(TScalar epsExponent, TScalar coeff) 
{
	
	assert(epsExponent==0||epsExponent==1);
	
	if ( epsExponent==0 )
		setX(coeff);
	else if ( epsExponent==1 )
		setEps(coeff);
	else
		assert(true==false);
	
}


template <typename TScalar>
TScalar 	 basicNumber<TScalar>::getValue(TScalar epsExponent)  const
{
	if ( epsExponent==0 )
		return getX();
	else if ( epsExponent==1 )
		return getEps();
	else
		return 0;
}


template <typename TScalar>
inline void basicNumber<TScalar>::setX(TScalar xxx)
{
	x=xxx;
}

template <typename TScalar>
inline void basicNumber<TScalar>::setEps(TScalar _eps) 
{
	eps=_eps;
}




//-------------------------Operators--------------------------------------------------------

template <typename TScalar>
inline  bool basicNumber<TScalar>::isZero() const
{
	return (x==0 && eps==0);
};


template <typename TScalar>
inline  bool basicNumber<TScalar>::isNotZero() const
{
	return (x!=0 || eps!=0);
};

///comparison; ignores EPS components.
template <typename TScalar>
inline int basicNumber<TScalar>::nearlyEqual(const basicNumber z) const
{
	return  (x==z.x);
};

template <typename TScalar>
inline bool basicNumber<TScalar>::operator==( const basicNumber z) const
{	
	return (x==z.getX() && eps==z.getEps());
}

template <typename TScalar>
inline bool basicNumber<TScalar>::operator!=( const basicNumber z)  const
{	
	return (x!=z.getX() || eps!=z.getEps());
}


/// @todo wieso hast du den Operator += geschrieben? Das birgt Fehleranfälligkeit!
/*
template <typename TScalar>
inline void  basicNumber<TScalar>::operator+=(const basicNumber &z)	
{
	//#ifdef SAFE
		std::cerr <<"Warning: using unsafe +=operator!" << std::endl;
	//#endif
	x+=z.x;
	eps+=z.eps;
}
*/



template <typename TScalar>
inline void  basicNumber<TScalar>::printMultSecure(std::ostream &os)	const
{
	os <<  "(" << *this << ")";
	 
}



//-------------------------Index Computation--------------------------------------------------------




/// @note all index functions must me inter-coordinated
template <typename TScalar>
size_t 	
basicNumber<TScalar>::getPairIndex     (	const basicNumber a, 
							const basicNumber b, 
							const TScalar 	characteristic)
{
	//return ( (a.eps*characteristic+a.x)*characteristic + b.eps )*characteristic + b.x  ;
	return ( (a.eps*characteristic+b.eps) * characteristic + a.x ) * characteristic + b.x  ;
}

template <typename TScalar>
size_t 	
basicNumber<TScalar>::getPairIndexByRef(	const basicNumber & a, 
							const basicNumber & b, 
							const TScalar 	& characteristic)
{
	//return ( (a.eps*characteristic+a.x)*characteristic + b.eps )*characteristic + b.x  ;
	return ( (a.eps*characteristic+b.eps) * characteristic + a.x ) * characteristic + b.x  ;
}

template <typename TScalar>
 size_t 	basicNumber<TScalar>::getSingleIndex(const basicNumber b, 
								 const TScalar  	characteristic)
{
	//return  b.eps *characteristic + b.x  ;
	return  b.eps *characteristic*characteristic + b.x  ;
}

template <typename TScalar>
 size_t 	basicNumber<TScalar>::getSingleIndexByRef(const basicNumber & b, 
									const TScalar 	& characteristic)
{
	//return  b.eps *characteristic + b.x  ;
	return  b.eps *characteristic*characteristic + b.x  ;
}

template <typename TScalar>
 size_t   basicNumber<TScalar>::getMaxSingleIndex(const TScalar characteristic)
{
	return (characteristic-1)*characteristic*characteristic + characteristic-1;
}

template <typename TScalar>
 size_t   basicNumber<TScalar>::getMaxPairIndex  (const TScalar characteristic)
{
	return characteristic*characteristic*characteristic*characteristic-1;
}




//-------------------------IO--------------------------------------------------------


template <typename TScalar>
std::ostream &  operator<<(std::ostream & out, const basicNumber<TScalar>& z)  
{	
	out << (int)z.getX();

	if (!z.getEps()==0)
	{
		out << " +  ";
		out << (int)z.getEps() << "*eps^1 " ;
	}
	return out;
}

































//----------------------fieldScalar----------------------------------------




























template <typename TScalar>
 const fieldScalar<TScalar> fieldScalar<TScalar>::Zero = fieldScalar<TScalar>(0);

template <typename TScalar>
 const fieldScalar<TScalar> fieldScalar<TScalar>::One = fieldScalar<TScalar>(1);



//------------------------------Constructors----------------------------------------
template <typename TScalar>
inline   fieldScalar<TScalar>::fieldScalar()  : x(0) {};

// x initialisierung vergessen...
template <typename TScalar>
inline   fieldScalar<TScalar>:: fieldScalar(std::stringstream & sstream)
{
	x=0;
	sstream >> x;
	#ifdef DEBUG
	std::cerr << "parse fieldScalar from stream" << x << std::endl;	
	#endif
}

/// constructs afieldScalar {s  }
template <typename TScalar>
inline fieldScalar<TScalar>::fieldScalar (TScalar s) : x(s)	 {	 };


/// constructs a fieldScalar, epsPrecision is ignored and assumed =0
template <typename TScalar>
inline fieldScalar<TScalar>::fieldScalar (TScalar epsPrecision, std::string s ) :x(0)
{
	if (epsPrecision!=0)
	{
		std::cerr << " allowed epsPrecision is { 0 }" << std::endl;
		throw(" allowed epsPrecision is { 0 }");
	}
	//assert(epsPrecision==0);
	
};

/// @note pow call is constructed to be compilable with gcc 3.x
template <typename TScalar>
inline bool fieldScalar<TScalar>::wellDefined(unsigned int characteristic)
{

	// std::cerr << "wellDefined: characteristic" << characteristic << std::endl;

	assert(sizeof(long long)>sizeof(TScalar));	//otherwise the following test does not work

	long long 	test  = pow( 2.0, (int)sizeof(TScalar)*8 )-1;
	TScalar test2 = test;
	test = test2;

	if (test<0)
	{
		//std::cerr << "TScalar is signed!" << std::endl;
		assert(pow( 2.0, (int)sizeof(TScalar)*8-1 )>characteristic);
	}
	else
	{
		std::cerr << "warning: TScalar is unsigned!" << std::endl;
		assert(pow( 2.0, (int)sizeof(TScalar)*8 )>characteristic);
	}
	return true;
}



/// copy constructor
template <typename TScalar>
inline fieldScalar<TScalar>::fieldScalar (const  fieldScalar<TScalar> & z) 	: x(z.x)	 {};


	
/// constructs a basic number {s  }
template <typename TScalar>
inline fieldScalar<TScalar>::fieldScalar (TScalar s, TScalar t) 	: x(s) 	 {};


/////////////////////////////////-Properties-------------------------------------------------
/// returns maximum possible epsPrecision
template <typename TScalar>
inline unsigned short fieldScalar<TScalar>::getEpsPrecision() const
{
	return 0;
}

template <typename TScalar>
inline short fieldScalar<TScalar>::getEpsDegree() const
{
	if ( getX() !=0 ) 
		return 0;
	else return -1;
}


//-----------------------------Value Access---------------------------------------------------

template <typename TScalar>
inline   TScalar   fieldScalar<TScalar>::getX() const 
{
	return x;
}


/// implementation only for compatibility
template <typename TScalar>
inline  TScalar   fieldScalar<TScalar>::getEps() const
{
	return (0);
}


template <typename TScalar>
inline void fieldScalar<TScalar>::setX(TScalar xxx)
{
	x=xxx;
}


template <typename TScalar>
inline void fieldScalar<TScalar>::setEps(TScalar _eps) 
{
	if (_eps!=0 ) 
		assert(true==false);
}


/// set coeff of e^epsPrecision to 'coeff'. 
///TODO: error checking @todo setValue-Operation nur über den Körper/ Ring laufen lassen?
template <typename TScalar>
inline  void   fieldScalar<TScalar>::setValue(TScalar epsPrecision, TScalar coeff) 
{
	assert(epsPrecision==0);
	
	if (epsPrecision==0)
		setX(coeff);
	
}

template <typename TScalar>
inline TScalar 	 fieldScalar<TScalar>::getValue(TScalar epsPrecision) const
{
	if (epsPrecision==0)
		return x;
	return 0;
}


template <typename TScalar>
inline TScalar& fieldScalar<TScalar>::operator[](int i)
{
	#ifdef SAFE
		assert(i==0);
	#endif

	if (i==0)
		return x;
	return 0;
	
}


//-----------------------------------Operators-------------------------------------------------

template <typename TScalar>
inline  bool fieldScalar<TScalar>::isZero() const
{
	return (x==0);
};


template <typename TScalar>
inline  bool fieldScalar<TScalar>::isNotZero() const
{
	return (x!=0 );
};


/// @todo wo zum Teufel hast du diesen Operator eingesetzt? Todo: testen, ob nicht aus Versehen zwei 'fieldScalar'-Objekte unbedarft addiert werden können.
/*
template <typename TScalar>
inline void  fieldScalar<TScalar>::operator+=(const fieldScalar &z)	
{
	#ifdef SAFE
		std::cerr <<"Warning: using unsafe +=operator!" << std::endl;
	#endif
	x+=z.x;
}*/



template <typename TScalar>
inline bool fieldScalar<TScalar>::operator==( const fieldScalar z)  const
{	
	return ( x==z.x );
}

template <typename TScalar>
inline bool fieldScalar<TScalar>::operator!=( const fieldScalar z)  const
{	
	return ( x!=z.x );
}


/// nearlyEqual-comparison; ignores EPS components.
template <typename TScalar>
inline int fieldScalar<TScalar>::nearlyEqual(const fieldScalar z) const
{
	return  ( x==z.x );
};


template <typename TScalar>
inline void 	fieldScalar<TScalar>::printMultSecure(std::ostream &os)	const
{
	os <<   (int)getX() ;
	 
}




//--------------------Index computation-------------------------------------------------

template <typename TScalar>
inline size_t 	
fieldScalar<TScalar>::getPairIndex     (	const fieldScalar a, 
							const fieldScalar b, 
							const TScalar 	characteristic)
{
	return a.x*characteristic + b.x;
}

template <typename TScalar>
inline  size_t 	
fieldScalar<TScalar>::getPairIndexByRef(	const fieldScalar & a,
							const fieldScalar & b, 
							const TScalar 	& characteristic)
{
	return a.x*characteristic + b.x;
}

template <typename TScalar>
inline size_t 	
fieldScalar<TScalar>::getSingleIndex	(const fieldScalar b, const TScalar characteristic)
{
	return getSingleIndex(b);
}


template <typename TScalar>
inline size_t 	
fieldScalar<TScalar>::getSingleIndexByRef	(const fieldScalar & b, const TScalar & characteristic)
{
	return getSingleIndexByRef(b);
}


template <typename TScalar>
inline size_t 	fieldScalar<TScalar>::getSingleIndex		(const fieldScalar b)
{
	return b.x;
}


template <typename TScalar>
inline size_t 	fieldScalar<TScalar>::getSingleIndexByRef	(const fieldScalar & b)
{
	return b.x;
}


template <typename TScalar>
 size_t   fieldScalar<TScalar>::getMaxSingleIndex(const TScalar characteristic)
{
	return characteristic - 1;
}

template <typename TScalar>
 size_t  fieldScalar<TScalar>::getMaxPairIndex  (const TScalar characteristic)
{
	return characteristic*characteristic - 1;
}




template <typename TScalar>
std::ostream &  operator<<(std::ostream & out, const fieldScalar<TScalar>& z)  
{	
	out << (int)z.getX();

	return out;
}



