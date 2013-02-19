
#ifndef basicNumber_h
#define basicNumber_h



#include <vector>
#include <iostream>
#include <assert.h>


#include "typedefs.h"
#include <math.h>

#define STR(X) #X
#define SW_STATUS(v1)  STR(v1)





/** \file basicNumber.h 
*
* @brief  contains class representing elements of finite Field F_q and finite Ring F_q[epsilon]. 
		See also fastNumber.h
*
* @author Martin Cremer, redesign Jakob Kröker
*
* @todo besseren (sprechenden) Namen einfallen lassen.
* @todo rename basicNumber (to what?)?
*
*/





/**
*
* @brief basicNumber represents elements of finite Field F_q and  finite Ring F_q[epsilon].  
* <br> 
* Limitation: field element representatives assumed as positive and only values 
		between {0...characterristic-1}
* 
@note if data design changes, please update memsetClearAllowed()-function
*

*
*
* @todo implizite TNum Initialisierung nur über den dazugehörigen Ring - mal sehen 
*
* @todo Skalar Datentyp abhängig von der Characteristik wählen - ist prinzipiell moeglich
* 
* @todo folgende Frage klaeren:<cstdint>) and use the 'int_fastX_t' 
*
* @todo als TScalar nur primitive Datentypen zulassen und dies prüfen. Wieso? -> 
*/
template <typename TScalar>
class basicNumber
{

protected:

	  TScalar x; 	///< represents <b>x*EPS^0</b>
	  TScalar eps; 	///< represents <b>eps</b>*EPS^1

 public:

 	/** @name static data
         * @{ */
		typedef TScalar scalarType; // ja? dann sollte man auch fieldScalar als TScalar uebergeben!
	
		/// represents zero (<b>x</b>=0 and <b>eps</b>=0)
		static  const 	basicNumber Zero;//=basicNumber(0);
		
		/// represents one (<b>x</b>=1 and <b>eps</b>=0)
		static  const 	basicNumber One;//=basicNumber(1);
	/** @} */
	
 	/** @name safety
         * @{ */

		/// returns true, if it is allowed to initialise class objects with memset(0)
		static inline bool memsetClearAllowed() {	return true; }
		/** @brief checks, if the datatype TScalar is large enough (in bits) to store 
				field element representatives {0...characterristic-1}

			in the case, that TScalar is signed, sizeof(TScalar) must be one bit greater 
			than nessesary to store positive field element representatives {0...characterristic-1}
		*/
		static inline bool wellDefined(unsigned int characteristic);
	/** @} */

	 /** @name Constructors
         * @{ */
		/// constructs a basic number {0 + 0*EPS }
		inline   basicNumber();
	
		/// constructs a basic number {s + 0*EPS }
		inline basicNumber (TScalar s);
	
		/// constructs a basic number, epsPrecision is ignored and assumed =1
		inline basicNumber (TScalar epsPrecision, std::string s );
	
		/// constructs a basic number {s + t*EPS }
		inline basicNumber (TScalar s, TScalar t);
		
		/// copy constructor
		inline basicNumber (const  basicNumber & z);
 	
	/** @} */

	 /** @name properties
         * @{ */
		inline 	unsigned short	getEpsPrecision() const;
		inline 	void	setEpsPrecision(int epsPrecision) const {	assert(false); 	};
		/// returns the highest EpsExponent where the coeffitient is not zero.
		inline  short	getEpsDegree() const;
	/** @} */


	/** @name data Access
         * @{ */		
		inline  TScalar	 	getX()  	const;	///< returns  <b>x</b>
		inline  TScalar		getEps()  	const;	///< returns  <b>eps</b>

		inline  void 		setX  (TScalar xxx);	///< set <b>x</b>= <b>xxx</b>
		inline  void 		setEps(TScalar _eps) ;///< set <b>eps</b>= <b>_eps</b>

		/** @brief  get  coefficient with defined EPSPrecision to  {<b>val</b> + <b>coeff</b>*epsPrcision } */
		inline  TScalar 	getValue(TScalar epsPrecision) const ;

		/** @brief 	set  coefficient with defined EPSPrecision to  {<b>val</b> + <b>coeff</b>*epsPrcision } */
		inline  void 		setValue(TScalar epsPrecision, TScalar coeff) ;
	/** @} */


	/** @name getset
         * @{ */
		inline TScalar& 	operator[](int i);

		// inline void 		operator+=(const basicNumber &z);

	/** @} */


	/** @name operators
         * @{ */
		/// returs true if (<b>z</b> . x == this -> x);    eps is ignored!
		inline  int 		nearlyEqual(const basicNumber z) const;
	
		inline  bool		isZero() 	const;
		inline  bool		isNotZero()  	const;

		inline  bool 		operator==( const basicNumber z) const;
		inline  bool 		operator!=( const basicNumber z) const;
		
	/** @} */


	/** @name index computation
         * @{ */
		static size_t 	getPairIndex     (const basicNumber a,
								const basicNumber b,
								const TScalar 	characteristic);
	
		static size_t 	getPairIndexByRef(const basicNumber & a, 
								const basicNumber & b,
								const TScalar 	& characteristic);
	
		static size_t 	getSingleIndex	(const basicNumber b,
								 const TScalar  	 characteristic);

		static size_t 	getSingleIndexByRef(const basicNumber & 	b,
								  const TScalar 	  & characteristic);

		static size_t   	getMaxSingleIndex(  const TScalar 	characteristic);
		static size_t   	getMaxPairIndex  (  const TScalar 	characteristic);
	/** @} */

		void printMultSecure(std::ostream &os)	const;

};





/// @todo bei zusammengesetzten Typen bei getX() und getEPS keinen primitiven Datentyp zurueckliefern, 
/// sondern einen mit Memberfunktionen wie z.B. isZero.

/// @todo wenn eine Initialisierung von TNum unter Umgehung des Rings geschehen soll, 
///muss das über einen expliziten Konstruktor-Aufruf	erfolgen, der Fehlerfreiheit zuliebe!

/**  @brief class representing  elements of F_q. 
<br> Limitation: 
field element representatives assumed as positive and only values between {0...characterristic-1}
*/
template <typename TScalar>
class fieldScalar
{

protected:

	  TScalar x; 		///< represents <b>x</b>*EPS^0

 public:
	inline operator int() {return x;}
	//inline operator uint16_t() {return x;}
	//inline operator const uint16_t() const {return x;}
	
	//operator uint16_t() {return x;}
	//operator uint32_t() {return x;}

	 /** @name static data
         * @{ */
		typedef TScalar scalarType;
	
		/// represents zero (<b>x</b>=0 and <b>eps</b>=0)
		static  const 	fieldScalar Zero;//=fieldScalar(0);
		
		/// represents one (<b>x</b>=1 and <b>eps</b>=0)
		static  const 	fieldScalar One;//=fieldScalar(1);
	/** @} */


		/// returns true, if it is allowed to initialise class objects with memset(0)
		static inline bool memsetClearAllowed() {	return true; }

		/** @brief checks, if the datatype TScalar is large enough (in bits) to store field element 
				representatives {0...characterristic-1}

			in the case, that TScalar is signed, sizeof(TScalar) must be one bit greater 
			than nessesary to store positive field element representatives {0...characterristic-1}

			@todo 
		*/
		static inline bool wellDefined(unsigned int characteristic);
	

	 /** @name Constructors
         * @{ */

		inline   fieldScalar(std::stringstream  & sstream);
		/// constructs a basic number {0 + 0*EPS }
		inline   fieldScalar();
		//inline   fieldScalar(int bla):x((TScalar)bla) {};
	
		/// constructs a basic number {s + 0*EPS }
		inline fieldScalar (TScalar s);
		inline fieldScalar (TScalar s,TScalar eps) ; //{ s + 0*EPS }
	
		/// constructs a basic number, epsPrecision is ignored and assumed =1
		inline fieldScalar (TScalar epsPrecision, std::string s );
		
		/// copy constructor
		inline fieldScalar (const  fieldScalar & z);
 	
	/** @} */

	 /** @name properties
         * @{ */
	
		inline 	unsigned short	getEpsPrecision() const;

		/// returns the highest EpsExponent where the coeffitient is not zero.
		inline  short	getEpsDegree() const;
	/** @} */


	/** @name data Access
         * @{ */
		
		inline  TScalar	 	getX()   const;///< returns  <b>x</b>
		inline  TScalar		getEps()  const;///< returns  <b>0</b>

		inline  void 		setX  (TScalar xxx);	///< set <b>x</b>= <b>xxx</b>
		inline  void 		setEps(TScalar _eps) ;///<does nothing @TODO : throw Exception !


		/// get  coefficient with defined EPSPrecision to  {<b>val</b> + <b>coeff</b>*epsPrcision }
		inline  TScalar 	getValue(TScalar epsPrecision) const;

		/// set  value to <b>coeff</b> if _epsExponent==0, otherwise does nothing
		inline  void 		setValue(TScalar _epsExponent, TScalar coeff) ;
	/** @} */


	/** @name getset
         * @{ */
		inline TScalar& 	operator[](int i);
		//inline void 		operator+=(const fieldScalar &z);
	/** @} */


	/** @name operators
         * @{ */
		/// returs true if (<b>z</b> . x == this -> x);    eps is ignored!
		inline  int 		nearlyEqual(const fieldScalar z) const;
	
		inline  bool		isZero() 	 const;
		inline  bool		isNotZero()  	const;

		inline  bool 		operator==( const fieldScalar z) const;
		inline  bool 		operator!=( const fieldScalar z) const;
	/** @} */


	/** @name index computation
         * @{ */
		/// @todo const Correctness
		inline static size_t 	getSingleIndex		(const fieldScalar b,    
										 const TScalar 	characteristic  	);

		inline static size_t 	getSingleIndexByRef	(const fieldScalar & b,  
										const TScalar 	& characteristic	);

		inline static size_t 	getSingleIndex		(const fieldScalar b  );

		inline static size_t 	getSingleIndexByRef	(const fieldScalar & b);

		inline static size_t 	getPairIndex     ( const fieldScalar a,  
									 const fieldScalar b, 
									 const TScalar 	characteristic	);

		inline static size_t 	getPairIndexByRef(const fieldScalar & a,
									const fieldScalar & b,
									const TScalar 	&characteristic	);

		static size_t   getMaxSingleIndex(const TScalar characteristic	);

		static size_t   getMaxPairIndex  (const TScalar characteristic	);
	/** @} */

	void printMultSecure(std::ostream &os)	const;
};


	#include "basicNumber.cpp"

#endif // #ifndef basicNumber_h
