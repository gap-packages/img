
#ifndef POLYNOM_H4D45E74F99FC
#define POLYNOM_H4D45E74F99FC

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <fstream>
#include <iostream>
#include <assert.h>

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>



/** \file polynom.h
*
* @brief contains polynomial classes in one and two varibles with template coefficients
*/

//using namespace std;
using  std::cout;
using  std::cerr;
using  std::endl;
using  std::string;
using  std::vector;
using  std::ostream;
using  std::stringstream;

#include "parseTools.h"

template <class PolynomXY_Type, class RingType, class _istream>
PolynomXY_Type  createFromStream (_istream & _polynomStream, RingType & ringRef )
{

	
//	std::cerr << "_polynomStream: "<<  _polynomStream.str() << std::endl;

	std::string 		polynomStr = extractNextBracedData( _polynomStream );

	std::stringstream  	polynomStream( polynomStr);

	int 				mononGroupNum = countSubGroups( polynomStream );

	//std::cerr << "momonroups: "<<  mononGroupNum << std::endl;

	
	int 	MaxMonomDegree = getMaxMonomDegree( polynomStr );
	vector<int> degreeTracking(MaxMonomDegree + 1);
	
	PolynomXY_Type res(MaxMonomDegree);
	
	
	polynomStream.clear();
	polynomStream.seekg(0);
	
	extractChar(polynomStream, '{');
	//std::cerr << "matrixStream" << matrixStream.str() << std::endl;

	for ( int currMonomGroupNum = 1;  currMonomGroupNum <= mononGroupNum ;  currMonomGroupNum ++ )
	{
		std::string  monomGroup   = extractNextBracedData(polynomStream);
		std::stringstream ssMonomGroup(monomGroup);

		int currDegree = countElements(ssMonomGroup) -1;
		//std::cerr << "currDegree" << currDegree << std::endl;	
		ssMonomGroup.clear();
		ssMonomGroup.seekg(0);
		//safety:
		assert(degreeTracking[currDegree] == 0);
		degreeTracking[currDegree]  = 1 ;

		//std::cerr << "monomGroup" << monomGroup << std::endl;	
			
		extractChar(ssMonomGroup, '{');

		for ( int y_exp=0;  y_exp <= currDegree ;  y_exp++)
		{
			// extract next unbraced Element on current brace Level (finished by ',' or '}')
			string nextElement=extractNextData(ssMonomGroup);

		//	std::cerr << "nextElement" << nextElement << std::endl;	

			typename PolynomXY_Type::CoefficientType  coeff = parseNumber( &ringRef, nextElement.c_str(), ringRef.getEpsPrecision() );
		
			res.setCoeff( currDegree - y_exp, y_exp, coeff );

			if (y_exp !=currDegree ) 
				extractChar( ssMonomGroup, ',' );
			else extractChar( ssMonomGroup, '}' );
		}
		if (currMonomGroupNum !=mononGroupNum ) 
				extractChar( polynomStream, ',' );
		else extractChar( polynomStream, '}' );
	}
	//std::cerr << "polynom initialization OK" << std::endl;
	return res;
}

template <class PolynomXY_Type, class RingType>
PolynomXY_Type  createFromString (string & _polynomString, RingType & ringRef )
{
	stringstream strstream   (   _polynomString);
	return createFromStream<PolynomXY_Type>( strstream, ringRef);
}

/** @brief class used to store the coefficients of a polynom in (x,y)
*
* @TODO eventuell auch einen Ring hier speichern, da im Programm nicht allzuviele Polynome verwendet werden.
* @TODO soweit Verallgemeinern, dass die Klasse mit einer getIndex-Fkt parametrisiert werden kann. Dann kannste fast_polynom wegwerfen.
*/
template <class TNum, class TIndex>
class polynomXY  
{
private:
	 short 	maxDegree; ///< max possible degree of a contained (x,y)-monom
	 short 	maxDegreePlusOne;
	TIndex   idefs;
	 short 	size; 	///< monom coefficients count 
	TNum *	koeff;	///<koeff[x_exp*(maxDegree+1)+y_exp] = value;

	//TNum * *	offsets;
	string 	name;	/// object name
	

public:

	typedef TNum 	CoefficientType;

	/** @name Contsructors / Destructors
	 @{ */
		polynomXY();
		polynomXY(const short _maxDegree);
	
		polynomXY(string _name, const short _maxDegree);
	
		/// copy constructor
		polynomXY(const polynomXY&);

		~polynomXY();
 	/** @} */


	/** @name Initialization
	 @{ */
		inline void 	setDegree(short maxDegree);

		inline void 	clear(short maxDegree);
		inline void 	clear();
	/** @} */


	/** @name Data access
	 @{ */
		inline void 		setCoeff(const short x_exp,const short y_exp, const TNum value);
	
		inline TNum 		getCoeff(const short x_exp,const short y_exp) 	const;
	
		inline TNum const 	getCoeffConst(const short x_exp,const short y_exp) const;
	
		inline const TNum &  	getCoeffConstRef(const int x_exp,const int y_exp) const;
	
		inline TNum& 		getCoeffRef(const short x_exp, const short y_exp);

		inline TNum const * 	getCoeffConstAddr(const short x_exp, const short y_exp) const;
		inline TNum * 		getCoeffAddr(const short x_exp, const short y_exp);
	/** @} */


	/** @name Properties
	 @{ */
		inline  short 	getDegree() const;
		
		inline string 	getName() const {	return name;	};

		/// return max monom degree capacity
		inline short 	getMaxDegree() const;
	/** @} */

	


	/** @name Safety
	 @{ */
		inline void testBounds(const short x_exp, const short y_exp) const;
	/** @} */
	
	/** @name operators
	 @{ */
		/// assignment operator
		polynomXY& operator=(const polynomXY&);
		bool operator==(const polynomXY&) const;
	/** @} */
	//void update(std::string s);
	//void update(std::string s);

	/** @name IO
	 @{ */
		void OutputPureCoefficients( ostream &datei, int maxDegree, bool mitKomma) const;
		void output( std::ostream& os) const;
		void printInMacaulayStyle( std::ostream& os) const;
		void outputMatrix() const;
		void print( std::ostream& os) const;
	/** @} */
	private:
		inline int  getIndex(const short x_exp, const short y_exp) const;

};





/** \class polynomx
*
* @brief represents a polynom in one variable
*
*
*

*
*/
template <class TNum>
class polynomx
{

    public:
        static const polynomx   One;
        static const polynomx   Zero;

	private:
		 int 	maxDegree;
		 int 	size;
		TNum *		koeff; ///< data; koeff[x_exp] = value;

        
	protected:
		/** @name Safety
		@{ */
			inline void 		testbounds(const int x_exp) const;
		/** @} */
	public:
	
		typedef TNum 	CoefficientType;


        std::string     getVariableName() const  {   return "x"; };

        std::string      getStringRep() const;

		/** @name Constructors / Destructors
		@{ */
			polynomx();
			polynomx(const int gr);

            polynomx(const int gr, bool monic);
			inline  polynomx(const polynomx & fpx);
			inline  polynomx(const std::vector<TNum> & vecx);

			virtual ~polynomx();

		/** @} */

			inline polynomx & operator=(const polynomx & fpx);

			inline bool  operator==(const polynomx & fpx) const;
			inline bool  operator!=(const polynomx & fpx) const;

            inline const TNum  & operator[](const int x_exp) const
            {
                assert ( x_exp>=0 && x_exp <= maxDegree );
                return koeff[ x_exp ];
            }

	
		/** @name init
		@{ */
			void 	clear(int _degree);
			void 	clear();
		/** @} */
		

		/** @name Data access
		@{ */
                
			inline TNum 		getCoeff(const int x_exp)  const;	

			inline TNum 		getSafeCoeff(const int x_exp)  const;
		
			inline  TNum& 		getCoeffRef(const int x_exp) ;

			inline  const TNum& 	getCoeffConstRef(const int x_exp) const;

			inline  TNum const 	getCoeffConst(const int x_exp)  const;

			inline  TNum const 	getSafeCoeffConst(const int x_exp)  const;
		
			inline void 		setCoeff(const int x_exp, const TNum& value);
	
		/** @} */


		/** @name Properties
		 @{ */
			inline const int 	getDegree() const ;
            inline const int    getInitialDegree() const { return getDegree() ; }
			inline const int 	getExactDegree() const ;

			inline const bool isZero() const ;

			inline const bool isOne() const ;

			inline const bool isConstant() const ;

            // das problem lag in der include-reihenfolge  !!! warum ist inline gefährlich??
            static polynomx<TNum>  getOne()   
            {
                #ifdef DEBUG
                std::cerr << "polynomx<TNum>  getOne()  " << std::endl;
                #endif
                polynomx pol= polynomx(0);
                // manchmal nicht korrekt 
                pol.setCoeff(0, polynomx<TNum>::CoefficientType::One);
                assert(1 == pol.getCoeff(0).getX() );

                return pol;
            }

            static polynomx<TNum>  getZero()   
            {
                polynomx pol= polynomx(0);
                
                pol.setCoeff(0,TNum::Zero);
            
                return pol;
            }

		/** @} */

		void print( std::ostream& os) const;
		
	
       template <class TCoeffRing>
       inline bool nextInPlace(const TCoeffRing& coeffRing,  bool ignoreHighestCoeff=false )  
        {

            polynomx *   result = this;
            
            TNum   currCoeff = TNum::Zero;
            int degree = getInitialDegree();
            if ( ignoreHighestCoeff )
                degree--;
            
            int pos=0;
            while (currCoeff==TNum::Zero && pos<=degree)
            {                
                coeffRing.addInPlace( getCoeffRef(pos), TNum::One );
                currCoeff = getCoeff( pos );
                pos++;
            };
            if (pos>degree && currCoeff==TNum::Zero)
                return false;
            return true;           
        }
	

};



template <typename TNum, typename TIndex>
std::ostream &  operator<<(std::ostream & out, const polynomXY<TNum, TIndex>& z)  
{	
	z.print(out);
	return out;
}


template <typename TNum>
std::ostream &  operator<<(std::ostream & out, const polynomx<TNum>& z)  
{	
	z.print(out);
	return out;
}


#include "polynom.cpp"

#endif // ifndef POLYNOM_H4D45E74F99FC
