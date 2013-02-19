
#ifndef FAST_RING_H
#define FAST_RING_H


#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "fastNumber.h"




/** @file fast_Ring.h
*
* @brief (optimized) implementation of a finite field F_q and finite ring F_q[epsilon] template for small characteristic.<br>

        Characteristic of the field/ring may be defined at compile time.
	*TODO Ring mit Characteristic = 0!
	
*/

/** @file fast_Ring.cpp
*
* @brief 
*/

using namespace std;

// zum Design:
// bei Koerper fehlt istElement(zahl)
// generic_Ring sollte sowas bieten, wie getZahl


/**
* @brief Implements some fast basic operations for elemens of a ring  F_q[epsilon] or field F_q with small characteristic; <br>
	 For all operations the user must be sure that the operators are regular represantants of ring elements, in doubt by calling 
	 the Convert() function. Currently all represantant values must lie between 0 and (characteristic-1).
	 Allowed 'eps' precisions are 0 (is a field) and 1 (is a ring) <br>
         Parametrized with element Type (TNum) and characteristic of the field the ring is based on during compile time

 Main idea: table lookup instead of computation. 
For big characteristics this method costs performance, because the operation lookup tables does not fit in L1 or L2 cache
But for small characteristics it works very well.

To avoid copying data, functions with reference parameters were implemented, because it is impossible to enforce inlining.
But that costs implementation overhead and is dangerous, because is is possible to pass a reference to a temporary variable.
and woult result in an error.
*
* Die Schnelligkeit der Berechnungen beruht auf 
* Verwendung von Operationstabellen und von der Annahme, dass alle Funktionsaufrufe mit gültigen Parametern ausgeführt werden.
* Die Operationstabellen sind so ausgelegt, dass mit Hilfe von Bit Shift Operationen das Ergebnis
* einer Rechenoperation abgelesen werden kann:
* statt 
*
* ArrayIndex(a,b)=a*charakteristik+b  
*
* ArrayIndex(a,b)=a*Zweierpotenz+b benutzt wird.; Zweierpotenz>=charakteristik,
* Somit kann die Indexfunktion mit Bit Schifts und bitweisem 'OR'  berechnet werden
*  Der Nachteil ist, dass mehr Speicherplatz, als benötigt belegt wird.
*
* Die Berechnungen werden wieder langsamer, wenn die Tabellen nicht in den Cache passen.
* Bei Charakteristik 29 und dem Zahlentype 'number_eps1' beträgt der Speicherplatzbedarf
* 29^4*2 Bytes * 2 Tabellen (Multiplikation, Addition) = 4 MByte
* und bei Charakteristik 23 1119364 = 1 MByte
*
*@note Komisch, dass die 23-er Berechnung nicht wesentlich schneller ist  - wieso sollte sie? es gibt nur 1 strudelgroesse weniger!
*
* @note in the optimized and error-free version it is not tested, if an operation result is really a legal ring/field representative.
* Moreover, in a correct frommer algorithm implementation this test is not neccessary.
*
*
* @todo multiplyOnCPU einfuehren bzw multiplyByHand
*
* @todo Idee: shift index bei generischer Version - hätte vielleicht Vorteile, da schneller addiert werden kann. 
* Hm, warum wird die Optimierte Version ohne Addition nicht schneller? Schliesslich wird der Platz für den Pointer
*  auf die Additionstabelle frei.


@todo ueberpruefen, ob man nicht ohne die funktionen mit parametern 'passed by reference' auskommen kann,
d.h. ob dies nicht wesentlich die performance verschlechtert, und wie man bei 'pass by value'-parametern 
und Rückgabewerten unnötiges Kopieren verhindern kann.

\ingroup Algebra
\ingroup FieldRing


*/
template <class TNum, class kdefs>
class fast_Ring 
{
public:
	typedef 	struct sqrtInf 
			{
				short solutions;
				TNum	sqrt;
				sqrtInf(short _sol,TNum _sqrt): solutions(_sol), sqrt(_sqrt)
				{ }
			sqrtInf(): solutions(0), sqrt(TNum::Zero)
				{ }
	} sqrtInf_t;


	

	typedef 	TNum 					ElementType;

	typedef 	fast_Ring< TNum,kdefs  >		FieldType;

	typedef 	typename FieldType::ElementType		ScalarType;



private:

	const unsigned short characteristic;

	long 	moduloTableSize_m;

	bool  bContainsImagNum_m;

	

	/// @todo Variable umbenenen in epsPrecision.
	const unsigned short epsilon;
	
	TNum 	imagNum_m;///@todo sollte eigentlich 'const' sein

	const TNum 	generator;
	const TNum* 	elementsToExponentsTab;
	const TNum* 	exponentsToElementTab;

	const TNum* 	additiveInverseTable;
	const TNum*  	multiplicationTable;
	const TNum*  	additionTable;
	const TNum*  	multiplicativeInverseTable;
	const TNum* 	fastAdditionTable;

	const ScalarType * moduloTable;
	const sqrtInf_t* sqrtTable;


protected:
	inline void	init() ;

 public:
	inline TNum getZero() const {  return TNum::Zero; };

    inline  TNum getOne() const {  return TNum::One; };
    
    
	/// containsImag  wird  in  createMultiplicationTable() initialisiert
	bool 		containsImagNum()	const		{	return 	bContainsImagNum_m;	}

	/// imagNum_m	 wird  in  createMultiplicationTable() initialisiert
	TNum 		getImagNum()	const		{	assert(containsImagNum() );	return 	imagNum_m;	}


	bool 		wellDefined();

	/// @todo alternativ Field asls const Membervariable.
	inline const FieldType * getField() const
	{
		return this;
	}

	inline const FieldType & getFieldRef() const
	{
		return *this;
	}

	/** @name Constructors / Destructors
	 @{ */
		fast_Ring(unsigned short _characteristic, unsigned short epsPrec); 

        fast_Ring(unsigned short _characteristic, unsigned short epsPrec, short generator ); 
		
		~fast_Ring();
	/** @} */
	

	/** @name properties
	 @{ */

		inline unsigned short  getCharacteristic() const	{	return characteristic;	}

        inline unsigned short  getCardinality() const    
        {
              assert(epsilon==0);
              return characteristic;  
        }

		inline const  unsigned short & getCharacteristicRef() const	{	return characteristic;	}
	
		inline unsigned short  getEpsPrecision() const 		{	return epsilon;		}

		inline void  setEpsPrecision(int epsPrecision) const 		
		{
			assert(epsPrecision<=epsilon);
			//if (epsPrecision<2)
			//return epsilon;		
			return;
		}

	/** @} */

	
	/** @name safety
	 @{ */
		inline bool 	 isValid(TNum a ) const		{ assert(a ==Convert(a) );	return true;	}

        bool isGenerator(const TNum & _generator) const;
	
/** @} */
	//-----------Additionsfunktionen
	

	/** @name add
         * @{ */
		inline 		TNum 	add		( const TNum a,	const  TNum b) 	const;	
		inline 		TNum	addRef	( const TNum &a,const  TNum& b) const;
	/** @} */

	/** @name add in place
	 @{ */
		inline void	addInPlace   		(TNum& a, const TNum  b) 	const;
		inline void	addInPlaceRef		(TNum& a, const TNum & b) 	const;
 	/** @} */

	//-----------Multiplikationsfunktionen

	/** @name additive inverse
	 @{ */
		inline  TNum		addInv	 ( const TNum   a ) const;
		inline  TNum		addInvRef( const TNum & a ) const;

		inline  void		addInvInPlace(  TNum & a ) const;
	/** @} */

	/** @name accMult
	 @{ */
		inline void 	accMult( TNum& a ,const TNum b , const TNum c)  const;
		inline void 	accMultRef( TNum& a ,const TNum& b , const TNum& c)  const;

		inline void 	accMult( TNum* a ,const TNum b , const TNum c)  const;
		inline void 	accMultAddr( TNum* a ,const TNum* b , const TNum* c) const;

		inline void 	accMultSpec( TNum* const a ,const TNum b , const TNum * const c) const;
	
	/** @} */

	/** @name multiplication
	 @{ */
		inline TNum    		multiply(const TNum a, const TNum b) const;
		inline TNum 		multiplyRef(const TNum& a, const TNum& b) const;

	/** @} */


	/** @name multiply in place
         * @{ */
		inline void 	multiplyInPlace( TNum& a, const TNum b) const;
		inline void 	multiplyInPlaceRef( TNum& a, const TNum& b) const;
	/** @} */


	/** @name scalar multiplication
	 @{ */
		inline TNum    		scalarMultiply   (const FieldType::ElementType a, const TNum b) const;
		inline TNum 		scalarMultiplyRef(const FieldType::ElementType& a, const TNum& b) const;

	/** @} */


	/** @name scalar multiply in place
         * @{ */
		inline void 	scalarMultiplyInPlace   (const FieldType::ElementType a,  TNum & b) const;
		inline void 	scalarMultiplyInPlaceRef(const FieldType::ElementType& a,  TNum& b) const;
	/** @} */

	/** @name multiplication by Exponents
	 @{ */
		inline TNum const   	multByExp   (const TNum a, const TNum b) const;
		inline TNum const  	multByExpRef(const TNum & a, const TNum & b) const;
		inline void  		multByExpInPlace( TNum & a, const TNum  b) const;
		inline void  		multByExpInPlaceRef( TNum & a, const TNum & b) const;
	/** @} */


	/** @name multiplicative inverse
         * @{ */
		inline  TNum  		multInv   (const TNum a) const;
		inline  TNum 		multInvRef(const TNum &a) const;
		inline void		multInvInPlace( TNum &a) const;
	/** @} */


	/** @name power
         * @{ */
		inline TNum 	pow	  ( TNum const a ,unsigned int exp) const;
		inline void 	powInPlace( TNum &     a ,unsigned int exp) const;
	/** @} */

	/** @name sqrd
         * @{ */
		inline 		sqrtInf 	sqrt	( const TNum a)  const;	
		inline 		sqrtInf		sqrtRef	( const TNum &a) const;
	/** @} */

	/** @name Conversion
	 @{ */

        // todo: es fehlt eigentlich RepToInt...

        inline int  repToInt(const  TNum &)    const;

		inline  typename FieldType::ElementType  lookupModuloTable(int convertee) const;

		inline int getLookupModuloTableSize() const;

		template <class TConvNum >
		inline TNum	Convert(const TConvNum a) const;

		//template <class TConvNum =double >
		///@TODO Convert eventuell umbenennen ConvertScalarToRingElement 
		inline TNum	Convert(const double a) const;
		///@TODO Convert eventuell umbenennen ConvertScalarToRingElement 
		inline TNum	Convert(const int a) const;
        inline TNum Convert(const short a) const;
        inline  TNum Convert(const unsigned long  a) const;

		template <class TConvNum >
		inline void	ConvertInPlace( TConvNum & a) const;
		
		///@TODO einfuehren ConvertScalarToFieldElement 
		inline int	ConvertScalar(const int a) const;

		inline int 	ConvertScalarSpec(const int  a) const;
		
		inline int 	FastConvertScalar(const int  a) const;

        inline int elemToGeneratorExponent(const TNum z1) const
        {
            return ( elementsToExponentsTab[z1.getX() ] ).getX() ;
        }

        // todo: ueberarbeiten!!!
        inline int generatorExponentToElem(const TNum z1) const
        {
            return ( exponentsToElementTab[ z1.getX() ] ).getX() ;
        }

      
	/** @} */

	/** @name Table index computation
	 @{ */
		inline  size_t 	getMaxSingleIndex() const;
		inline  size_t 	getMaxPairIndex() const;
		inline  size_t 	getSingleIndex(const TNum z1) const;
		inline  size_t 	getSingleIndexByRef(const TNum &z1) const;

		inline  size_t 	getPairIndexByRef(const TNum &z1, const TNum& z2) const;
		inline  size_t 	getPairIndex(const TNum z1, const TNum z2) const;
	/** @} */

	/// only for the case, that table pointers are static
	//inline static void  clear() { }; /* gibt den Index in der Multiplikations- bzw. Additionstabelle zurck */


 protected:
	
	TNum* 	createAdditionTable(); /* initialisiert Additionstafel */

   
	TNum*	createMultiplicationTable(); /* initialisiert Multiplikationstafel */	

	TNum* 	createSubtractionTable(); /* initialisiert SubtraAdditionstafel */
	TNum* 	createAdditiveInverseTable();
	
	TNum* 	createMultiplicativeInverseTable(); /* initialisiert Tabelle fr multiplikative Inverse */

	/// Returns a generator of the Ring\{0} if exists, otherwise TNum::Zero
 	TNum 		getGenerator();

	TNum* 	initElementsToExponentsTab(TNum erzeuger);
	TNum* 	initExponentsToElementTab(TNum erzeuger);

	TNum* 	createFastAdditionTable();

	typename FieldType::ElementType* 	createModuloTable(); /* initialisiert Additionstafel */

	sqrtInf_t* createSqrtTable();

};



/// @todo Optimization: Addition with tables is faster on a intel Processor for epsprecision 1. 0 ?
/// Auch auf hoech ist auf einmal die Addition mit Tabellen schneller... 

/// @todo diese 'assert(Convert(a)==a)' können wieder raus, wenn ein TNum nur über den Ring->convert(int,int) erzeugt wird
/// @todo den !=NULL-Test fuer die Tabellen braucht man eingentlich nicht - das wuerde bei einem valgrind-Lauf sofort auffliegen

/*
template <class TNum>
inline int generic_Ring<TNum >::ConvertScalar(const int a) const
{
	int res = a;
	while (res<0)
	{
		res += kdefs::charakteristik;
	}
	if ( res >= kdefs::charakteristik ) 
	{	
		res %= kdefs::charakteristik;
	}
	return res;
}*/




	#include "fast_Ring.cpp"

#endif
