#ifndef improved_zahl
#define improved_zahl



#include "CompileFunctions.h"
#include <math.h>
#define STR(X) #X
#define SW_STATUS(v1)  STR(v1)

/** \file fastNumber.h
*
* @brief  contains optimized (packed) datatypes for elements of  finite Field  F_q[] and finite Ring  F_q[epsilon]. <br>Field characteristic is static and have to be defined during compile time. See also basicNumber.h.
*
* number_eps1(epsPRecision=1) ,<br>
* number_eps0 (epsPRecision=0),<br>
* 
* @todo The interfaces of  number_eps1 nnumber_eps0 and epsZahl should be identical! 
* @note inheritance was not used to define an explizit interface, because this would slow down the program
* there exists some complicated template methods to define a interface without performance loss, but that is on the todo list
*
@todo:  Fehlertest für die ausreichende Dimensierung der Datentype (insb. TScalar) bereits zur Kompilezeit durchführen
*

* @note  implementation note: Error: cannot declare member function ‘static int Foo::bar()’ to have static linkage
if you declare a method to be static in your .cc file.
The reason is that static means something different inside .cc files than in class declarations 
It is really stupid, but the keyword static has three different meanings. In the .cc file, 
the static keyword means that the function isn't visible to any code outside of that particular file.
This means that you shouldn't use static in a .cc 
*/

using namespace std;



/**
* @brief  Class representing elements of F_q[epsilon];<br> 
*  compact datatype for (x,eps)-Pair; field characteristic is parametrized during compile time
*
*datatype for (x,eps)-Pair; 
-----------number_eps1 DATA LAYOUT<br>
// Data is stored compact in a 'short int' as follows: <br>
// first 'bitsize' bits for 'x' and next 'bitsize' bits for 'eps'<br>
// MSB ist the first bit of 'x'<br>
// LSB ist the last bit of 'eps'<br>
Der Template-Parameter CHAR ist die Charakteristik des Koerpers, dessen Werte in number_eps1.x , bzw number_eps1.eps<br>
dargestellbar sein sollen. CHAR ist auf 256 beschraenkt! <br>
*
* @todo security check: Template-Parameter CHAR <=sizeof(TScalar) and (CHAR*CHAR)<=sizeof(TScalarPair) !
*

* @todo statt (int getPairIndex) (size_t getPairIndex)?
*
* @todo herausfinden, warum und um wieviel uebergabe per referenz in manchen Fällen schneller ist.
* basierend auf diesem Wissen könnte  man sich die Funktionen getxxxByRef eventzell ersparen!
*/
template <int CHAR, typename TScalar, typename TScalarPair>
struct number_eps1
{
 public:

		typedef 	TScalar scalarType; 	

	/** @name  value representation
         * @{ */
		TScalarPair 	epsx; ///< data 
	/** @} */ 

	/** @name static data
         * @{ */
		static  const number_eps1 	Zero;//=number_eps1(0);
		static  const number_eps1 	One;//=number_eps1(1);
	
		enum { bitsize		= needbits<CHAR>::value		}; ///< number of reserved bits for x or eps   in epsx
		enum { fullbitsize	= needbits<CHAR>::doubledvalue	}; ///< number of reserved bits for x AND eps  in epsx
	
		enum { maskx	= pow2<needbits<CHAR>::value>::valueMinusOne	}; ///< bitmask for  x
		
		/** @brief  bitmask for  eps*/
		enum { maskeps	= pow2<needbits<CHAR>::doubledvalue>::valueMinusOne - ( pow2<needbits<CHAR>::value>::valueMinusOne)	};
	/** @} */ 

	/** @name safety
         * @{ */ 

		/// returns true, if it is allowed to initialise class objects with memset(0)
		static inline bool memsetClearAllowed() {	return true; }
		
		static inline bool wellDefined(unsigned int characteristic);
		
		/// @todo Grenzen fuer enums ueberpruefen
		static inline bool wellDefined();
	/** @} */

	/** @name Constructors
         * @{ */ 

		inline number_eps1();
	
		/**  @brief einheitlicher Konstruktor fuer eine Zahl mit beliebigen EpsPrecision - notwendig ???
		*
		* @todo allgemeiner Konstruktor fuer eine Zahl mit belibigen EpsPrecision:
		* Problem: eine 0 wird als leerer String interpretiert, bzw, ein char als ein Integer-Wert
		* Insgesamt noch keine gute Loesung	*/
		inline number_eps1(int epsPrec, string dummy) ;
	
		/** @brief constructs (x,0) from x */
		inline number_eps1 (scalarType x);
		
		/** @brief copy constructor */
		inline number_eps1 (const  number_eps1 & z_x) ;
		
		/** @brief constructs x + eps*EPS  */
		inline number_eps1 (scalarType x, scalarType eps);
 	/** @} */ 

	/** @name properties
         * @{ */
		inline unsigned short	getEpsPrecision() const 	{	return 1;	}; ///< returns 1

		inline void	setEpsPrecision(int epsPrecision) const 	{	assert( epsPrecision==1);	}; ///< returns 1

		/// returns the highest EpsExponent where the Coeffitient is not zero.
		inline  short	getEpsDegree() const;
		
	/** @} */ 


	/** @name data access
         * @{ */
		inline scalarType 	getX() const ;	///< get x from  (x,eps)
		inline scalarType 	getEps() const;	///< get eps from (x,eps)
	
		inline void 		setX(scalarType x)  ;///<  set x in (x,eps)
		inline void 		setEps(scalarType eps) ;///< set eps in (x,eps)
		inline void 		setValue( unsigned short _epsPrecision, scalarType val);

		inline int 		getValue( unsigned short _epsPrecision) 
		{	
			if(_epsPrecision==0) 		return getX(); 
			else if (_epsPrecision==1) 	return getEps();
			else	return 0;
		 };
	/** @} */ 


	/** @name operators
         * @{ */
		inline  bool 	isZero() const	;
		inline  bool 	isNotZero() const	;

		/** @brief 	compare  x-components, eps-parts are ignored */
	 	inline bool 	nearlyEqual(const number_eps1 & z) const;

		inline int 	operator==(const  number_eps1<CHAR, TScalar, TScalarPair> & z) const;
		inline int 	operator!=(const  number_eps1<CHAR, TScalar, TScalarPair> & z) const;
	/** @} */ 


	/** @name regular index computation interface
         * @{ */
		/** @brief 		returns  b.epsx	*/
		static 
		inline size_t 	getSingleIndex	   (   const number_eps1<CHAR, TScalar,TScalarPair>  b,
									 const TScalar 				characteristic ) ;
		
		/** @brief 		returns  b.epsx	*/
		static 
		inline size_t 	getSingleIndexByRef( 	const number_eps1<CHAR, TScalar,TScalarPair> & b,
									 const TScalar & characteristic ) ;

		/** @brief 		computes a*(2^(fullbitsize)) + b.epsx	*/
		static 
		inline size_t	getPairIndex	 (	const number_eps1<CHAR, TScalar,TScalarPair> a, 
									const number_eps1<CHAR, TScalar,TScalarPair> b, 
									const TScalar characteristic) ;
	
		/** @brief 		computes a*(2^(fullbitsize)) + b.epsx	*/
		static 
		inline size_t	getPairIndexByRef(	const number_eps1<CHAR, TScalar,TScalarPair>& a, 
									const number_eps1<CHAR, TScalar,TScalarPair>& b,
									const TScalar & characteristic) ;
	
		
		static inline size_t 	getMaxSingleIndex(const TScalar characteristic) ;
		static inline size_t 	getMaxPairIndex  (const TScalar characteristic) ;
 	/** @} */ 

	/** @name reduced index computation interface
         * @{ */
		/** @brief 		returns  b.epsx	*/
		static inline size_t	getSingleIndex		( const number_eps1<CHAR, TScalar,TScalarPair> b ) ;
		
		/** @brief 		returns  b.epsx	*/
		static inline size_t 	getSingleIndexByRef	( const number_eps1<CHAR, TScalar,TScalarPair>& b ) ;

		/** @brief 		computes a*(2^(fullbitsize)) + b.epsx	*/ 
		static inline size_t 	getPairIndex	 (	const number_eps1<CHAR, TScalar,TScalarPair> a,
										const number_eps1<CHAR, TScalar,TScalarPair> b) ;
	
		/** @brief 		computes a*(2^(fullbitsize)) + b.epsx	*/
		static inline size_t 	getPairIndexByRef(	const number_eps1<CHAR, TScalar,TScalarPair>& a, 
										const number_eps1<CHAR, TScalar,TScalarPair>& b) ;
	
	/** @} */ 
};





/**
* @brief Datatype for an element of an finite field with characteristik 'CHAR'.
*
*datatype for (x,); 
-----------number_eps0 DATA LAYOUT--------------------<br>
// Data is stored  in a 'TScalar' as follows: <br>
// first 'bitsize' bits for 'x' and no  'bits' for epsilon!

Template parameter CHAR determines the maximal possible regular value. <br>
The TScalar is the internal type for the stored value.

Der Template-Parameter CHAR ist die Charakteristik des Koerpers, dessen Werte in number_eps0.x
dargestellbar sein sollen.
*  maximal zulaessiger CHAR-Wert  haengt vom Template-Datentyp 'datatype' ab!
*
* @todo Sicherheitspruefung: maximal zulaessiger Template-Parameter CHAR <=sizeof(datatype)*8!
*
** @todo warnen, wenn als DataType ein unsigned ding verwendet wird 
*
* @todo einmalig warnen wenn  (bitsize>=sizeof(TScalar)*8) - 
*      dann kann nicht im Speicher subtrahiert werden

* @todo einmalig warnen wenn  (bitsize>=(sizeof(TScalar)-1)*8) - 
*      dann kann nicht im Speicher addiert werden (es kommt zum Überlauf)
*
* @todo size_t or NOT size_t as return type for index funtions?
*
 @todo in den Indexfunktionen eventuell bitsize durch  needbits<CHAR>::value ersetzen, weil schneller!
*
*/
template <int CHAR, class TScalar> 
struct number_eps0
{
 public:
	
	//operator int() {return x;}

	TScalar x; ///< number_eps0  data

	typedef TScalar scalarType; 	///< number datatype

	static  const number_eps0 	Zero;
	static  const number_eps0 	One;

	/// @todo 3rd template parameter with delivers reserved bits for stored data
	enum{  bitsize = needbits<CHAR>::value }; 


 	/** @name Constructors
         * @{ */ 
		inline number_eps0() ; ///< value is set to ::Zero;

		inline number_eps0 (TScalar _x) ;///<create  number_eps0, .x ist set to _x. 
		
	 	inline number_eps0 (TScalar _x, TScalar eps) ;///<create  number_eps0, .x ist set to _x, eps is ignored
	 
	 	inline number_eps0(int epsPrec, std::string dummy);///<create  number_eps0 with epsprecision=epsPrec. fails if epsPrec is not 0 
	/** @} */
	
	/** @name safety
         * @{ */ 
		
		/// returns true, if it is allowed to initialise class objects with memset(0)
		static inline bool memsetClearAllowed() {	return true; }

		static inline bool wellDefined();

		static inline bool wellDefined(unsigned int characteristic)
		{
			return wellDefined();
		}
	/** @} */
	

	/** @name properties
         * @{ */ 
		
		inline unsigned short 	getEpsPrecision() const { 	return 0;  	};

		/// returns the highest EpsExponent where the coeffitient is not zero.
		inline  short		getEpsDegree() 	const;

	/** @} */

	/*inline explicit operator int () 	{	return x; 	} *///...
	
	/** @name index computation regular interface
         * @{ */
		/// returns  a.x*(2^bitsize) + b.x 
		 static inline  size_t 	getPairIndex	 (	const number_eps0<CHAR, TScalar>  a,
								 		const number_eps0<CHAR, TScalar> b, 
										const TScalar characteristic) ;

		/// returns  a.x*(2^bitsize) + b.x
	 	 static inline size_t 	getPairIndexByRef(const number_eps0<CHAR, TScalar> & a,
									const number_eps0<CHAR, TScalar> & b,
									const TScalar &  characteristic) ;
	
		/// returns   b.x
		 static inline size_t 	getSingleIndex		(const number_eps0<CHAR, TScalar>  b,
								 		const TScalar characteristic)  ;
		/// returns  b.x
		 static inline size_t 	getSingleIndexByRef	(const number_eps0<CHAR, TScalar> & b,
								 		const TScalar &  characteristic)  ;

		static inline size_t   getMaxSingleIndex(const TScalar characteristic);
		static inline size_t   getMaxPairIndex  (const TScalar characteristic);
	/** @} */

	/** @name index computation reduced interface
         * @{ */
		
		/// returns  a.x*(2^bitsize) + b.x
		 static inline  size_t 	getPairIndex	 (	const number_eps0<CHAR, TScalar>  a,
								 		const number_eps0<CHAR, TScalar> b) ;

		/// returns  a.x*(2^bitsize) + b.x
	 	 static inline size_t 	getPairIndexByRef(	const number_eps0<CHAR, TScalar> & a,
										const number_eps0<CHAR, TScalar> & b) ;
	
		/// returns   b.x
		 static inline size_t 	getSingleIndex		(const number_eps0<CHAR, TScalar>  b)  ;

		/// returns  b.x 
		 static inline size_t 	getSingleIndexByRef	(const number_eps0<CHAR, TScalar> & b)  ;

	/** @} */

	/** @name data access
         * @{ */ 
		inline TScalar 	getX() 				const  ;
		inline void 	setX(const TScalar _x) ;

		inline TScalar 	getEps() 				const ;

		inline void 	setEps(const TScalar eps) 	const ;

		///	sets entweder .x to coeff if _epsExp==0, otherwise does nothing
		inline void 	setValue( unsigned short _epsExp, 
							TScalar 	  coeff    );

		inline int 		getValue( unsigned short _epsExp) const;
	/** @} */
	

	/** @name operators
         * @{ */ 
		 inline  bool 	isZero()    const	 ;
	 	 inline  bool 	isNotZero() const  ;

		/// compare only the .x-Component
		inline bool nearlyEqual(const number_eps0 z) 	const  ;

		inline int 	operator==(const  number_eps0 & z)	const	 ;
		inline int 	operator!=(const  number_eps0 & z)	const	 ;
	
	/** @} */
};


	#include "fastNumber.cpp"


#endif
