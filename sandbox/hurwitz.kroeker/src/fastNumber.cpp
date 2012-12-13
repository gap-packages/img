
/*
#ifdef INTEL

#include "fastNumber.h




template <int CHAR,class TScalar> 
const number_eps0<CHAR,TScalar> number_eps0<CHAR,TScalar>::Zero(0);

template <int CHAR,class TScalar> 
const number_eps0<CHAR,TScalar> number_eps0<CHAR,TScalar>::One(1);

template <int CHAR,class TScalar> 
const int number_eps0<CHAR,TScalar>::bitsize(needbits<CHAR>::value);

template <int CHAR>
const number_eps1<CHAR> number_eps1<CHAR>::Zero=number_eps1<CHAR>(0);

template <int CHAR>
const number_eps1<CHAR> number_eps1<CHAR>::One=number_eps1<CHAR>(1);
*/







//-----------------------------------------------------number_eps1-----------------------------------


	template <int CHAR, typename TScalar, typename TScalarPair>
	const number_eps1<CHAR, TScalar, TScalarPair>
	number_eps1<CHAR, TScalar, TScalarPair>::Zero = number_eps1<CHAR, TScalar, TScalarPair>(0);
	
	template <int CHAR, typename TScalar, typename TScalarPair>
	const number_eps1<CHAR, TScalar, TScalarPair> 
	number_eps1<CHAR, TScalar, TScalarPair>::One = number_eps1<CHAR, TScalar, TScalarPair>(1);




//------------------------------------------ Constructors -------------------------------


	template <int CHAR, typename TScalar, typename TScalarPair>
	inline number_eps1<CHAR, TScalar, TScalarPair>::number_eps1() : epsx(0)	{	};



	/**
	@brief einheitlicher Konstruktor fuer eine Zahl mit beliebigen EpsPrecision
	*
	* @todo Konstruktor fuer eine Zahl mit belibigen EpsPrecision:
	* Problem: eine 0 wird als leerer String interpretiert, bzw, ein char als ein Integer-Wert
	* Insgesamt noch keine gute Loesung
	*/
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline number_eps1<CHAR, TScalar, TScalarPair>::
	number_eps1(int epsPrec, std::string dummy)		:	epsx(0)
	{

		if (epsPrec!=1)
		{
			std::cerr << " allowed epsPrecision is {  1 }" << std::endl;
			throw(" allowed epsPrecision is   {  1 }");
			//assert(epsPrec==1);
		}
	}

	///constructs (x,0) from x
	/// @todo wieso scheitert dieser Test bei der Kombination char/short und Charakteristik =29 ? - nein, scheitert nicht für Charakteristik=29 sondern für Charakteristik=197.
	template <int CHAR, typename TScalar, typename TScalarPair>
	 inline number_eps1<CHAR, TScalar, TScalarPair>::number_eps1 (TScalar x) :	epsx(x)
	 {
		#ifdef SAFE
			//cerr << " <CHAR> "<<  CHAR << endl;
			//cerr << "needbits<CHAR>::doubledvalue="<< needbits<CHAR>::doubledvalue << endl;
			//cerr << "sizeof(TScalarPair)*8="<< sizeof(TScalarPair)*8 << endl;
			 assert (needbits<CHAR>::doubledvalue<sizeof(TScalarPair)*8);
		#endif
	 };
	
	/// @todo:   bereits waehrend des Compilevorgangs prüfen, 
	/// ob die Datenstrukturen gross genug ausgelegt sind	 -  , es reicht aber auch, vor der ersten Verwendung die Tests Durchzuführen.
	///
	///@brief Konstruiert number_eps1 (x,eps) aus number_eps1(x,eps)
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline number_eps1<CHAR, TScalar, TScalarPair>::
	number_eps1 (const  number_eps1<CHAR, TScalar, TScalarPair> & z_x) :	epsx( z_x.epsx )
	 {
		#ifdef SAFE
			 assert (needbits<CHAR>::doubledvalue<sizeof(TScalarPair)*8);
		#endif
	 };
	 
	///Konstruiert number_eps1(x,eps) aus (x,eps)
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline number_eps1<CHAR, TScalar, TScalarPair>::
	number_eps1 (TScalar x, TScalar eps):  epsx((eps << bitsize)|x)//:epsx(s),eps(t)
	 {
		#ifdef SAFE
		  // todo: sollte bereits waehrend des Compilevorgangs geprft werden
			 assert(fullbitsize<sizeof(TScalarPair)*8);
		#endif
	 };

	template <int CHAR, typename TScalar, typename TScalarPair>
	inline  short	number_eps1<CHAR, TScalar, TScalarPair>::getEpsDegree() const
		{
			if ( getEps()!=0 )
				return 1;
			else if ( getX() !=0 ) 
				return 0;
			else return -1;
		}

//------------------------------------------ Safety-------------------------------


	template <int CHAR, typename TScalar, typename TScalarPair>
	inline bool number_eps1<CHAR, TScalar, TScalarPair>::wellDefined(unsigned int characteristic)
		{

		//	std::cerr << "wellDefined: characteristic" << characteristic << std::endl;
		return wellDefined();
	}
		
	/// @todo Grenzen fuer enums ueberpruefen
	template	 <int CHAR, typename TScalar, typename TScalarPair>
	 inline bool number_eps1<CHAR, TScalar, TScalarPair>::wellDefined()
	{
		assert(sizeof(long long)>sizeof(TScalar));	//otherwise the following test does not work

		/// trage in Tscalar test2 an jede zulässige Bitstelle Einsen ein.
		long long 	test  = pow( 2, sizeof(TScalar)*8 )-1;
		TScalar test2 = test;
		test = test2;

		if (test<0)
		{
		//	std::cerr << "TScalar is signed!" << std::endl;
			assert( sizeof(TScalar)*8-1>=bitsize);
		}
		else
		{
			std::cerr << "warning: TScalar is unsigned!" << std::endl;
			assert( sizeof(TScalar)*8>=bitsize);
		}

		assert(sizeof(long long)>sizeof(TScalarPair));	//otherwise the following test does not work

			test  = pow( 2, sizeof(TScalar)*8 )-1;
			test2 = test;
		test = test2;

		if (test<0)
		{
		//	std::cerr << "TScalarPair is signed!" << std::endl;
			assert( sizeof(TScalarPair)*8-1>=fullbitsize);
 			assert (needbits<CHAR>::doubledvalue<sizeof(TScalarPair)*8);
		}
		else
		{
			std::cerr << "warning: TScalarPair is unsigned!" << std::endl;
			assert( sizeof(TScalarPair)*8>=fullbitsize);
 			assert (needbits<CHAR>::doubledvalue<=sizeof(TScalarPair)*8);
		}
		return true;			
	}

//------------------------------------------ Data Access-------------------------------

	/// get x from (x,eps)
	template <int CHAR, typename TScalar, typename TScalarPair>

	inline TScalar number_eps1<CHAR, TScalar, TScalarPair>::getX() const 
	{
		return ( maskx&epsx );
	}

	/// get eps from (x,eps)
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline TScalar  	number_eps1<CHAR, TScalar, TScalarPair>::getEps() const
	{
		#ifdef COUNT
			bitwiseShift+=1;
		#endif	
		return ( epsx >> needbits<CHAR>::value );
	}

	
	/// set x in (x,eps)
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline void  	number_eps1<CHAR, TScalar, TScalarPair>::setX(TScalar x)  
	{
		#ifdef COUNT
			bitwiseOR+=1;
			bitwiseAND++;
		#endif	
		///todo: Assert, that x has not too much bits!
		epsx = ( epsx & maskeps ) | x ; 
	}

	/// setzt eps in (x,eps)
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline void 	number_eps1<CHAR, TScalar, TScalarPair>::setEps(TScalar eps) 
	{
		#ifdef COUNT
			bitwiseOR+=1;
			bitwiseShift+=1;
			bitwiseAND++;
		#endif	
		epsx = (epsx & maskx) | (TScalarPair) (eps << needbits<CHAR>::value ); 
	}
 
	
	/**
	* @brief setzt entweder .x oder .eps gleich 'val', je nach 'eps_exponent'-Wert
	*
	* @param _epsPrecision beinhaltet die epsilon--Potenz von val - mit val ist gemeint \f$ VAL=val*(eps^eps\_exponent) \f$
	* @param val der Koeffizient von  (\f$ eps^\_epsPrecision \f$)
	* @todo Funktionsbezeichnung ungluecklich, setMonom passt aber auch nicht.
	*/
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline void  	
	number_eps1<CHAR, TScalar, TScalarPair>::
	setValue(  unsigned short _epsPrecision, TScalar val)
	{
		if (_epsPrecision==0)
			setX(val);
		else if (_epsPrecision==1)
			setEps(val);
		else if ( val!=0 )
			assert(true==false);

			
	}

//------------------------------------------ Operators-------------------------------
	
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline  bool 	number_eps1<CHAR, TScalar, TScalarPair>::isZero() const	
	{	
		return (epsx==0);	
	};
	
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline  bool 	number_eps1<CHAR, TScalar, TScalarPair>::isNotZero() const	
	{
		return (epsx!=0);	
	};
	
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline bool 
	number_eps1<CHAR, TScalar, TScalarPair>::
	nearlyEqual(const number_eps1<CHAR, TScalar, TScalarPair> & z) const
	  {
		  return  ( getX()==z.getX() );
	  };

	
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline int  
	number_eps1<CHAR, TScalar, TScalarPair>::
	operator==(const  number_eps1<CHAR, TScalar, TScalarPair> & z) const
	{
		return epsx == z.epsx;
	}
	
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline int  
	number_eps1<CHAR, TScalar, TScalarPair>::
	operator!=(const  number_eps1<CHAR, TScalar, TScalarPair> & z) const
	{
		return epsx != z.epsx;
	}



//------------------------------------------ Index computing -------------------------------

	
	
	template <int CHAR, typename TScalar, typename TScalarPair>
	// inline size_t number_eps1<CHAR, TScalar, TScalarPair>::
	//getSingleIndex(const number_eps1<CHAR, TScalar, TScalarPair> b)
	inline size_t 	 
	number_eps1<CHAR, TScalar, TScalarPair>::getSingleIndex(const number_eps1<CHAR, TScalar, TScalarPair> b)
	{
		#ifdef COUNT
			bitwiseShift	+= 1;
			bitwiseOR	+= 1;
		#endif	
		// warum mal zwei? -> weil epsx den eos0 und eps1-Anteil speichert!
		return( (size_t) b.epsx);
	
	}

	
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline size_t 
	number_eps1<CHAR, TScalar, TScalarPair>::getSingleIndexByRef(
									const number_eps1<CHAR, TScalar, TScalarPair> & b )
	{
		#ifdef COUNT
			bitwiseShift	+= 1;
			bitwiseOR	+= 1;
		#endif	
		//  mal zwei weil epsx den eos0 und eps1-Anteil speichert!
		return( (size_t) b.epsx);
	}
	

	template <int CHAR, typename TScalar, typename TScalarPair>
	inline size_t
	number_eps1<CHAR, TScalar, TScalarPair>::getPairIndex(
								const number_eps1<CHAR, TScalar, TScalarPair> a,
								const number_eps1<CHAR, TScalar, TScalarPair> b)
	{
		#ifdef COUNT
			bitwiseShift	+= 1;
			bitwiseOR	+= 1;
		#endif	
		//  mal zwei, weil epsx den eps0 und eps1-Anteil speichert!
		return( ((size_t)a.epsx<< (needbits<CHAR>::doubledvalue)) | (size_t)b.epsx );
	}
	
	
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline size_t 
	number_eps1<CHAR, TScalar, TScalarPair>::getPairIndexByRef(
									const number_eps1<CHAR, TScalar, TScalarPair>& a,
									const number_eps1<CHAR, TScalar, TScalarPair>& b)
	{
		#ifdef COUNT
			bitwiseShift	+= 1;
			bitwiseOR	+= 1;
		#endif	
		//  mal zwei, weil epsx den eps0 und eps1-Anteil speichert!
		return( ((size_t)a.epsx<< (needbits<CHAR>::doubledvalue)) | (size_t) b.epsx );
	}


	template <int CHAR, typename TScalar, typename TScalarPair>
	inline size_t 
	number_eps1<CHAR, TScalar, TScalarPair>::
	getSingleIndex(const number_eps1<CHAR, TScalar, TScalarPair> b,
											const TScalar characteristic)
	{

		#ifdef SAFE
			assert(characteristic==CHAR);
		#endif 
		return number_eps1<CHAR, TScalar, TScalarPair>::getSingleIndex(b);
	}


	template <int CHAR, typename TScalar, typename TScalarPair>
	inline size_t 
	number_eps1<CHAR, TScalar, TScalarPair>::
	getSingleIndexByRef(		const number_eps1<CHAR, TScalar, TScalarPair>& b,
						const TScalar & characteristic)
	{

		#ifdef SAFE
			assert(characteristic==CHAR);
		#endif 
		return number_eps1<CHAR, TScalar, TScalarPair>::getSingleIndexByRef(b);
	}


	template <int CHAR, typename TScalar, typename TScalarPair>
	inline size_t 
	number_eps1<CHAR, TScalar, TScalarPair>::
	getPairIndex(			const number_eps1<CHAR, TScalar, TScalarPair> a, 
						const number_eps1<CHAR, TScalar, TScalarPair> b, 
						const TScalar  					characteristic	)
	{
		#ifdef SAFE
			assert(characteristic==CHAR);
		#endif 
		return number_eps1<CHAR, TScalar, TScalarPair>::getPairIndex(a,b);
	}


	template <int CHAR, typename TScalar, typename TScalarPair>
	inline size_t 
	number_eps1<CHAR, TScalar, TScalarPair>::
	getPairIndexByRef(const number_eps1<CHAR, TScalar, TScalarPair>	& a,
				const number_eps1<CHAR, TScalar, TScalarPair>	& b, 
				const TScalar 						& characteristic	)
	{
		#ifdef SAFE
			assert(characteristic==CHAR);
		#endif 
		return number_eps1<CHAR, TScalar, TScalarPair>::getPairIndexByRef(a, b);
	}




	template <int CHAR, typename TScalar, typename TScalarPair>
	inline size_t   
	number_eps1<CHAR, TScalar, TScalarPair>::
	getMaxSingleIndex(const TScalar characteristic)
	{
		//#ifdef SAFE
			assert(characteristic==CHAR);
		//#endif 

		number_eps1<CHAR, TScalar, TScalarPair> num (characteristic-1,characteristic-1);

		return number_eps1<CHAR, TScalar, TScalarPair>::getSingleIndex(num);
		//return nextpow2num<CHAR>::value*nextpow2num<CHAR>::value;
	}

	
	template <int CHAR, typename TScalar, typename TScalarPair>
	inline size_t   
	number_eps1<CHAR, TScalar, TScalarPair>::
	getMaxPairIndex  (const TScalar characteristic)
	{
		//#ifdef SAFE
		if (characteristic!=CHAR)
		{
			std::cerr << "characteristic" << characteristic <<  std::endl;
			std::cerr << "CHAR" << CHAR << std::endl;
		}
			assert(characteristic==CHAR);
		//#endif 
		number_eps1<CHAR, TScalar, TScalarPair> num1 (characteristic-1, characteristic-1);
		number_eps1<CHAR, TScalar, TScalarPair> num2 (characteristic-1, characteristic-1);
	
		return number_eps1<CHAR, TScalar, TScalarPair>::getPairIndex(num1, num2);
		//return nextpow2num<CHAR>::value*nextpow2num<CHAR>::value*nextpow2num<CHAR>::value*nextpow2num<CHAR>::value;
	}
	
//------------------------------------------ IO -------------------------------

	
	template <int CHAR, typename TScalar, typename TScalarPair>
	std::ostream &  
	operator<<(std::ostream & out, const number_eps1<CHAR, TScalar, TScalarPair>& z)
	{  
	
		out << (int)z.getX() ;
		if (z.getEps()!=0)
		{
			out << " +  ";
			out << (int)z.getEps()  << "*eps " ;
		}
		return out;
	} ;


























 //-----------------------------(number_eps0)--------------------------------------------



	template <int CHAR, class TScalar> 
	const number_eps0<CHAR, TScalar> 
	number_eps0<CHAR, TScalar>::Zero(0);
	
	
	template <int CHAR, class TScalar> 
	const number_eps0<CHAR, TScalar> 
	number_eps0<CHAR, TScalar>::One(1);




//-----------------------------------Constructors-----------------------------------------------------


	template <int CHAR, class TScalar> 
	inline 
	number_eps0< CHAR, TScalar >::number_eps0() :x(0)
	{
		//x=0;
	};
	
	/// @todo: es sollte waehrend des Compilevorgangs oder zumindest einmal während des Programmstarts geprueft werden,
	/// ob der interne Datentyp groß  genug ausgelegt ist.
	///
	/// create .x from (_x).
	template <int CHAR, class TScalar>  
	inline number_eps0< CHAR, TScalar >::number_eps0 (TScalar _x) :x(_x)
	{
	};

	
		///create .x from (_x,eps). eps is ignored
	template <int CHAR, class TScalar> 
	inline 
	number_eps0< CHAR, TScalar >::number_eps0 (TScalar _x, TScalar eps):  x(_x)
	{ 
	};
	
	
	/**
	@brief einheitlicher Konstruktor fuer eine Zahl mit belibigen EpsPrecision
	*
	* @todo Konstruktor fuer eine Zahl mit belibigen EpsPrecision:
	* Problem: eine 0 wird als leerer String interpretiert, bzw. ein char als ein Integer-Wert.
	* Insgesamt noch keine gute Loesung.
	*/
	template <int CHAR, class TScalar> 
	inline 	
	number_eps0< CHAR, TScalar >::number_eps0(int epsPrec, std::string dummy):x(0)
	{
		if (epsPrec!=0)
		{
			std::cerr << " allowed epsPrecision is { 0 }" << std::endl;
			throw(" allowed epsPrecision is { 0 }");
		}
		//assert(epsPrec==0);
	}
//--------------------------------------Safety-----------------------------------------------------
	
	
	template <int CHAR, class TScalar> 
	inline bool 
	number_eps0< CHAR, TScalar >::wellDefined()
	{
		assert(sizeof(long long)>sizeof(TScalar));	//otherwise the following test does not work
	
		long 	test  = pow( 2, sizeof(TScalar)*8 )-1;
		TScalar test2 = test;
		test = test2;
	
		if (test<0)
		{
		//	std::cerr << "TScalar is signed!" << std::endl;
			assert( sizeof(TScalar)*8-1>=bitsize);
		}
		else
		{
			std::cerr << "warning: TScalar is unsigned!" << std::endl;
			assert( sizeof(TScalar)*8>=bitsize);
		}
		return true;
	}
	
	
	template <int CHAR, class TScalar> 
	inline  short	
	number_eps0< CHAR, TScalar >::getEpsDegree() const
	{
		if ( getX() !=0 ) 
			return 0;
		else return -1;
	}

//--------------------------------------Data Access------------------------------------------------------
	
	template <int CHAR, class TScalar> 
	inline TScalar 	
	number_eps0< CHAR, TScalar >::getX() 	const 
	{
		return (x);	
	};

	template <int CHAR, class TScalar> 
	inline void 
	number_eps0< CHAR, TScalar >::setX(const TScalar _x)	
	{
		x = _x;	
	};

	/// returns 0
	template <int CHAR, class TScalar> 
	inline TScalar 	
	number_eps0< CHAR, TScalar >::getEps() 	const 
	{	
		  return 0; 
	}	


	template <int CHAR, class TScalar> 
	inline void 	
	number_eps0< CHAR, TScalar >::setEps(const TScalar eps) 	const 
	{
		 if (eps!=0)	
			assert(true==false); 	
	};


	template <int CHAR, class TScalar> 
	inline int 		
	number_eps0< CHAR, TScalar >::getValue( unsigned short _epsExp) const
	{
		if(_epsExp==0) 
			return x; 
		return 0;	
	};


	/**
	* @brief set  .x to 'coeff' if  _epsExponent==0, otherwise nothing
	*
	* @param _epsExponent beinhaltet die epsilon--Potenz von val -
	*  mit val ist gemeint VAL=coeff*(eps^eps_exponent)
	* @todo Funktionsbezeichnung ungluecklich, setMonom passt aber auch nicht
	*/
		
	template <int CHAR, class TScalar> 
	inline void number_eps0< CHAR, TScalar >::setValue( unsigned short _epsExponent, TScalar coeff)
	{
		if (_epsExponent==0)
			x = coeff;
		else  if ( coeff!=0 )
			assert(true==false);
	}
	


//------------------------------------ Index Computation ---------------------------------------------------
	
	template <int CHAR, class TScalar> 
	inline size_t 	
	number_eps0< CHAR, TScalar >::getSingleIndex(const number_eps0< CHAR, TScalar >  b)  
	{
		return( (size_t)b.x );
	}
	
	
	template <int CHAR, class TScalar> 
	inline  size_t 	
	number_eps0< CHAR, TScalar >::getSingleIndexByRef(const number_eps0< CHAR, TScalar > & b)  
	{
		return( (size_t)b.x );
	}
	
	template <int CHAR, class TScalar> 
	inline  size_t 	
	number_eps0< CHAR, TScalar >::getPairIndex(const number_eps0< CHAR, TScalar >  a, 
											const number_eps0< CHAR, TScalar > b)  
	{
		#ifdef SAFE
			size_t res = a.x;
			res = ( (res<<bitsize) | b.x );
			assert( res==( ((size_t)a.x<< bitsize) | b.x) );
		#endif
		return( ((size_t)a.x<< needbits<CHAR>::value) | (size_t)b.x );
	}
	
	
	template <int CHAR,class TScalar> 
	inline  size_t 	
	number_eps0< CHAR, TScalar >::getPairIndexByRef(const number_eps0< CHAR, TScalar > & a, 
									const number_eps0< CHAR, TScalar > & b)  
	{
		#ifdef SAFE
			size_t res = a.x;
			res = ( (res<<bitsize) | b.x );
			assert( res==( ((size_t)a.x<<bitsize) | b.x ) );
		#endif
		return( ( (size_t)a.x<< needbits<CHAR>::value ) | (size_t)b.x );
	}
	
	
	template <int CHAR,class TScalar> 
	inline size_t 
	number_eps0< CHAR, TScalar >::getSingleIndex	(const number_eps0< CHAR, TScalar > b,
									const TScalar 			characteristic)
	{
		return getSingleIndex(b);
	}
	
	
	template <int CHAR,class TScalar> 
	inline size_t 	
	number_eps0< CHAR, TScalar >::getSingleIndexByRef(const number_eps0< CHAR, TScalar > & b, 
									const TScalar 			& characteristic)
	{
		return getSingleIndexByRef(b);
	}
	
	
	template <int CHAR, class TScalar> 
	inline size_t 
	number_eps0< CHAR, TScalar >::getPairIndex(const number_eps0< CHAR, TScalar > a, 
								const number_eps0< CHAR, TScalar > b, 
								const TScalar 		characteristic)
	{
		return getPairIndex(a,b);
	}
	
	
	template <int CHAR, class TScalar> 
	inline size_t 	
	number_eps0< CHAR, TScalar >::getPairIndexByRef( const number_eps0< CHAR, TScalar > & a,
									const number_eps0< CHAR, TScalar > & b, 
									const TScalar 				& characteristic)
	{
		return getPairIndexByRef(a,b);
	}
	
	
	template <int CHAR, class TScalar> 
	inline size_t   
	number_eps0< CHAR, TScalar >::getMaxSingleIndex(const TScalar characteristic)
	{
		return( (size_t)(characteristic-1) );
	}
	
	
	template <int CHAR, class TScalar> 
	inline size_t   
	number_eps0< CHAR, TScalar >::getMaxPairIndex  (const TScalar characteristic)
	{
		return  ( ( (size_t)(characteristic-1)<<bitsize) | (size_t)(characteristic-1) );
	}

//----------------------- Operators ---------------------------------------------------

	template <int CHAR, class TScalar> 
 	inline  bool 	
	number_eps0< CHAR, TScalar >::isZero()    const		
	{
		return ( x==0 );     
	};

		
	template <int CHAR, class TScalar> 
	inline  bool 	
	number_eps0< CHAR, TScalar >::isNotZero() const 		
	{
		return ( x!=0 );     
	};

		/// compare only the .x-Component
	template <int CHAR, class TScalar> 
	inline bool 
	number_eps0< CHAR, TScalar >::nearlyEqual(const number_eps0 z) 	const 
	{
		return ( x == z.x );    
	};


	template <int CHAR, class TScalar> 
	inline int 	
	number_eps0< CHAR, TScalar >::operator==(const  number_eps0 & z)	const	
	{	
	 	return ( x == z.x );	
	}


	template <int CHAR, class TScalar> 
	inline int 	
	number_eps0< CHAR, TScalar >::operator!=(const  number_eps0 & z)	const	
	{
		return ( x != z.x );	
	}

//----------------------- IO (number_eps0) ---------------------------------------------------


	template <int CHAR,class TScalar> 
	std::ostream &  
	operator<<(std::ostream & out, const number_eps0<CHAR, TScalar>& z)
	{  
	
		out << (int)z.getX() ;
		if (z.getEps()!=0)
		{
			out << " +  ";
			out << (int)z.getEps()  << "*eps " ;
		}
		return out;
	} ;




/*

#endif // ifdef INTEL

*/
