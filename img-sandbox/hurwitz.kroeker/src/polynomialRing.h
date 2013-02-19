#pragma once

/*
class expIterator
{
	short degree_m;
	short currDegree_m;
	pair<short>	xyexp;
	
	expIterator(polynomDegree)
	
	bool end()

	
	bool next()
	{
		if (yexp<currDegree)
	{
		yexp++;
		xexp--;
	}
	else 
	{
		currDegree++;
		if (currDegree>polynomialDegree)
			return false;
		else
		{
			xexp=currDegree;
			yexp=0;
		}
	}
	}

};*/


/// @note Polynomring in zwei Variablen - naja, fast ein Polynomring, Multiplikation und Division wird erst implementiert, wenn diese gebraucht wird.
/// TODO Variablentyp für Maximalen Polynomgrad muss TPolynomXY vorgeben!
/// TODO typedef tring und tpolynomxy ueberall gleich durchfuehren!
template<class TPolynomXY, class TRing>
class PolynomialRing
{
	
public :

	typedef 	TPolynomXY		PolynomXY;

	typedef 	TRing		RingType;

	
	const RingType & 	ring_ref_m;
	

	PolynomialRing(const RingType & ring);

	/** @name additive Inverse
	  @{ */

		TPolynomXY 		addInv		(const TPolynomXY  & polynom) const;
	
		TPolynomXY*		addInvReturnPtr	(const TPolynomXY  & polynom) const;
		
		void 			addInvInPlace	(	TPolynomXY & polynom) const;

	/** @} */


	/** @name add polynoms
	  @{ */
		TPolynomXY 		add(		const TPolynomXY & polynom1,
									const TPolynomXY & polynom2 ) const;
	
		TPolynomXY* 	addReturnPtr( 	const TPolynomXY & polynom1,
									const TPolynomXY & polynom2	) const;
	
		inline void 			addInPlace( 		TPolynomXY & polynom1,
									const 	TPolynomXY & polynom2) const;
	/** @} */


	/** @name scalar multiply
	  @{ */
		TPolynomXY 		scalarMultiply(	const typename TPolynomXY::CoefficientType 	scalar,
								const TPolynomXY 				& polynom			) const;
		
		inline void 		scalarMultiplyInPlace(const 	typename TPolynomXY::CoefficientType	 scalar,
									TPolynomXY 					& polynom		) const;
	
		TPolynomXY*		scalarMultiplyRetPtr(	const	typename TPolynomXY::CoefficientType 	scalar,
									const TPolynomXY 					& polynom	) const;
	/** @} */

	


	/** @name convert
	  @{ */

		//TPolynomXY	convert(	const typename TPolynomXY pxy ) const ;

		template <class TPolynomXY_Type>
		void 	convertInPlace( TPolynomXY_Type & pxy) const ;

		
		template <class TPolynomXY_SRC_Type, class TPolynomXY_DEST_Type>
		static void 	 copyPolynomWithGivenEpsPrecision(const TPolynomXY_SRC_Type  & srcPol,
								 TPolynomXY_DEST_Type  	&	 destPol, int epsPrecision);

	/** @} */
	
};

template<class TPolynomX, class TRing>
class UnivariatePolynomialRing
{
	
public :

	typedef 	TPolynomX		Element;
	typedef 	TPolynomX		ElementType;

	typedef 	typename TRing::ElementType	CoeffType;


	typedef 	TRing		RingType;

    typedef     TRing       CoeffRingType;
 

	const RingType & 	ring_ref_m;

	const TRing	&	getRing()	const	{	return ring_ref_m;};

    const TRing &   getCoeffRing()   const   {   return ring_ref_m;};

   inline std::string  coeffToGeneratorExponentGAPStr(const CoeffType z1) const;


	
	UnivariatePolynomialRing(const RingType & ring);

    void    normalizeInPlace(TPolynomX   &  ) const;


	TPolynomX 		addInv		(const TPolynomX  & polynom) const;

	TPolynomX 		add(		const TPolynomX & polynom1,
						const TPolynomX & polynom2 ) const;


	TPolynomX 		scalarMultiply(	const typename TPolynomX::CoefficientType 	scalar,
							const TPolynomX 				& polynom			) const;

 
	TPolynomX 		multiply(			const TPolynomX 				& polynom1,
								const TPolynomX 				& polynom2			) const;
								
	TPolynomX* 		multiplyPtr(			const TPolynomX 				& polynom1,
								const TPolynomX 				& polynom2			) const;
								
								
	TPolynomX 		square(				const TPolynomX 				& polynom1) const;

	// Returns (result,remainder)-pair
	pair<TPolynomX,TPolynomX> 		divide(	const TPolynomX 				& dividend,
							const TPolynomX 				& divisor
							) const;
	TPolynomX 		remainder(	const TPolynomX 				& dividend,
							const TPolynomX 				& divisor
							) const;


	TPolynomX 		gcd(				const TPolynomX 				& polynom1,
								const TPolynomX 				& polynom2			) const;
	TPolynomX 		fastgcd(				const TPolynomX 				& polynom1,
								const TPolynomX 				& polynom2			) const;

	bool 		fastgcdspec(				TPolynomX 				& polynom1,
								const TPolynomX 				& polynom2			,
								int minDegree) const;


	TPolynomX 		diff(				const TPolynomX 				& polynom1 			) const;

	void 			diffInPlace(	 TPolynomX 				& polynom1 					) const;

	TPolynomX 		pow(				const TPolynomX 				& polynom, unsigned int exponent 			) const;

    TPolynomX       subst(                const TPolynomX                 & polynom,          const TPolynomX                 & substPol          ) const;
	 

	CoeffType 		evalAt(const TPolynomX & polynom1 , CoeffType  el) const;

};


#include "polynomialRing.cpp"
