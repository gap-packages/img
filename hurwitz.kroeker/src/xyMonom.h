
#pragma once

#include <sstream>
#include <assert.h>
#include <iostream>
#include <algorithm>

#include "parseTools.h"

/** @file xyMonom.h 
*
*   @brief contains contains xyExpPair and xyMonom -  monom template representation in two variables 
*
*/

/// exponents of monom in (x,y)-> negative exponets not allowed!
struct xyExpPair
{
	/** @name Constructors
 	@{ */
		inline xyExpPair();

		inline xyExpPair( int _x_exp,  int _y_exp );

		xyExpPair(std::stringstream&  sstream);

	/** @} */

	/** @name Access to Exponents
 	@{ */
	
		inline void 		set( int _x_exp,  int _y_exp);
		inline unsigned int 	getXExp() 	const;
		inline unsigned int 	getYExp() 	const;
		
		/// returns ( x_exp + y_exp )
		inline unsigned int 	getDegree() const;

	/** @} */

	protected:
		unsigned int 	x_exp; ///< x - exponent 
		unsigned int 	y_exp; ///< y - exponent 

		unsigned int 	extractExplicitExponent(std::stringstream&  sstream);
};




/** @brief Represents a  monom template in two variables. Therefore negative exponents are not allowed
*
*/
template <class ccoeff>
struct xyMonom
{

	/** @name Constructors
 	@{ */
		xyMonom();

		xyMonom(ccoeff _coeff, int _x_exp, int _y_exp );
	
		/// constructs a monom from a input stringstream , see detailed description for format
		xyMonom(std::stringstream&  sstream);

		xyMonom(std::string str);
	/** @} */


	void createFromStream(std::stringstream&  monomListStream);

	/** @name Monom exponent access
 	@{ */
		inline unsigned int 	getDegree() const;
		
		inline unsigned int	getXExp() 	const;
		inline unsigned int 	getYExp() 	const;
		inline void 		setExp(int _x_exp, int _y_exp);
	/** @} */

	/** @name Monom coefficient access
 	@{ */
		inline ccoeff 	getCoeff() const;
		inline void 	setCoeff(ccoeff _coeff);
	/** @} */

	/** @name evaluate xyMonom 
 	@{ */
		template <class IMultiplyCoeffWithVariable, class IMultiplyVariables, class ResultType >
		inline ResultType 	substitute(	const typename IMultiplyVariables::ElementType & x, 
								const typename IMultiplyVariables::ElementType & y,
								 const IMultiplyCoeffWithVariable & imultCoeffXVariable,
								 const IMultiplyVariables & imultVariableXVariable
							 ) const;

		template <class IMultiplyRing >
		typename IMultiplyRing::ElementType substitute2(			
								const typename IMultiplyRing::ElementType 		& x, 
								const typename IMultiplyRing::ElementType 		& y,
								const 	IMultiplyRing	 				& imult
							 ) const;

	/** @} */

	protected:
		ccoeff 		coeff;	///< monom coeffitient
	
		xyExpPair 	exponents;///< monom exponents of x and y
	
		std::string 	monomString; ///< monom-representing string

};






template <class ccoeff>
std::ostream &  operator<<(std::ostream & out, const xyMonom<ccoeff>& xyMonomObj);

/// @todo Problem: Klammer kann noch nicht eingelesen werden. sollte auch nicht ausgegeben werden. 
/// Nachteil: Ohne Klammerverwendung ist die Lesbarkeit nicht gut.
/// @todo Herausfinden, wo es ein Problem bei der Ausgabe geben kann, wenn nicht geklammert wird.
template <class ccoeff>
std::ostream &  operator<<(std::ostream & out, const xyMonom<ccoeff>& xyMonomObj)
{  
	//out << xyMonomObj.getCoeff() << "*x^" << xyMonomObj.getXExp() << "*y^" << xyMonomObj.getYExp() << std::endl ;
	xyMonomObj.getCoeff().printMultSecure(out);
	out  << "*x^" << xyMonomObj.getXExp() << "*y^" << xyMonomObj.getYExp() << std::endl ;
	return out;
} ;






inline xyExpPair::xyExpPair()
{
	x_exp=0;
	y_exp=0;
}

inline xyExpPair::xyExpPair(int _x_exp, int _y_exp )
{
	assert(_x_exp>=0 && _y_exp>=0);

	x_exp=_x_exp;
	y_exp=_y_exp;
}

inline void xyExpPair::set(int _x_exp, int _y_exp )
{
	assert(_x_exp>=0 &&_y_exp>=0);
	x_exp=_x_exp;
	y_exp=_y_exp;
}

inline unsigned int xyExpPair::getDegree() const 
{
	return x_exp + y_exp;
}

inline unsigned int xyExpPair::getXExp() const 
{
	return x_exp ;
}

inline unsigned int xyExpPair::getYExp() const 
{
	return y_exp ;
}






template <class ccoeff>
xyMonom<ccoeff>::xyMonom()
{
	exponents.set(0,0);
	coeff=0;
};

template <class ccoeff>
inline unsigned int xyMonom<ccoeff>::getDegree() const
{
	return exponents.getDegree();
}


template <class ccoeff>
inline unsigned int xyMonom<ccoeff>::getXExp() const 
{
	return exponents.getXExp() ;
}



template <class ccoeff>
inline unsigned int xyMonom<ccoeff>::getYExp() const 
{
	return exponents.getYExp() ;
}



template <class ccoeff>
inline ccoeff xyMonom<ccoeff>::getCoeff() const 
{
	return coeff;
}

template <class ccoeff>
inline void  xyMonom<ccoeff>::setExp(int _x_exp, int _y_exp)
{
	assert(_x_exp>=0 && _y_exp>=0);
	exponents.set(_x_exp,_y_exp);
}

template <class ccoeff>
inline void xyMonom<ccoeff>::setCoeff(ccoeff _coeff)
{
	coeff=_coeff;
}

template <class ccoeff>
void xyMonom<ccoeff>::createFromStream(std::stringstream&  monomListStream)
{
	#ifdef DEBUG
		std::cerr << "xyMonom:: " << std::endl;
		std::cerr << "monomListStream= '" << monomListStream.str() << "'";
	#endif
	ccoeff parsedCoeff;
	bool failed = false;

	std::streampos pos = monomListStream.tellp ( );
	
	bool negative = false;

	// eigentlich muss man sich bis x oder y vorarbeiten. dann ist dieser Teil der Koeffizient.
	// Wenn es nur ein '+' oder ein '-' ist, den Koeffizienten mit 1 oder -1 initialisieren. 
	try{
		#ifdef DEBUG
			std::cerr << std::endl <<"parsedCoeff = ccoeff(monomListStream);" << std::endl;
		#endif
		monomListStream.peek();
		assert( monomListStream.good() );
		if (monomListStream.peek()=='-')
		{
			negative=true;
			extractChar('-',monomListStream);
		}
		monomListStream.peek();

		if ( !monomListStream.eof() )
		{
			parsedCoeff = ccoeff(monomListStream);
	
			if (monomListStream.fail() )
			{
				// hilft alles nix!
				monomListStream.clear();
				//std::cerr   <<"pos  "<< pos << std::endl;
				monomListStream.seekp(pos);
				assert(! monomListStream.eof() );
			
				failed = true;
				throw "failed to get the monom coefficient";
			}
		}
		else
		{
			parsedCoeff=-1;
			negative=false;
		}

	}
	catch(char const * error)
	{
		// hier muss noch geprueft werden, ob ein 
		#ifdef DEBUG
			std::cerr << "char const * error: failed=true" << std::endl;
		#endif
		failed = true;
	}
	catch(std::bad_exception &e) 
	{
		// hier muss noch geprueft werden, ob ein 
		#ifdef DEBUG
			std::cerr << "bad_exception: failed=true" << std::endl;
		#endif
		failed = true;
	}
	catch(...)
	{
		// hier muss noch geprueft werden, ob ein 
		#ifdef DEBUG
			std::cerr << "exception (...) : failed=true" << std::endl;
		#endif
		failed = true;
	}
	if (failed)
	{

		#ifdef DEBUG
		std::cerr << " xyMonom: get coeff from stream failed" << std::endl;
		#endif
	
		//std::cerr << "monomListStream" <<  monomListStream.str() << std::endl;
		assert(! monomListStream.eof() );
		assert(! monomListStream.fail() );
		assert(! monomListStream.bad() );
		assert( monomListStream.good() );
		monomListStream.seekp(pos);	// hilft alles nix!
		//std::cerr << " monomListStream.peek()" <<  (char)monomListStream.peek() << "'" << std::endl;
		if (!monomListStream.eof() && (monomListStream.peek()=='x' || monomListStream.peek()=='y' ||  monomListStream.peek()=='-') )
		{
			#ifdef DEBUG
				std::cerr << "parsedCoeff=1;"<< std::endl;
			#endif
			parsedCoeff=1;
			if (monomListStream.peek()=='-')
			{
				parsedCoeff=-1;
				extractChar('-',monomListStream);
			}
		}
		else
			throw "failed to get  monom coefficient";
	}

	coeff = parsedCoeff;

	if (negative)
	{
		coeff= -parsedCoeff;
	}
	
	#ifdef DEBUG	
		std::cerr << "monomListStream= '" << monomListStream.str() << "'";
	#endif
	
	xyExpPair parsedxyExpPair(monomListStream);

	exponents = parsedxyExpPair;
}


template <class ccoeff>
xyMonom<ccoeff>::xyMonom(std::string   monom )
{
	std::stringstream monomstrstream(monom);

	createFromStream(monomstrstream);
}




/**  @brief reads xy-monom from input stream. <br>
/// For coefficient format conventions see class ccoeff and <br>
/// for variable exponent format conventions xyExpPair(stringtream)
@pre: coeff kann negiert werden (implementiert den 'operator-'
*/
template <class ccoeff>
xyMonom<ccoeff>::xyMonom(std::stringstream&  monomListStream)
{
	createFromStream(monomListStream);
};



template <class ccoeff>
xyMonom<ccoeff>::xyMonom(ccoeff _coeff, int _x_exp, int _y_exp )
{
	coeff=_coeff;
	exponents.set(_x_exp,_y_exp);
}


/// @todo strenggenommen müsste das Zwischenergebnis immer wissen, zu welchem Ring es gehört. Dies ist aber nicht
/// optimierungsfreundlich 
template <class ccoeff>
template <class IMultiplyCoeffWithVariable, class IMultiplyVariables, class ResultType >
ResultType 	xyMonom<ccoeff>::substitute(			
								const typename IMultiplyVariables::ElementType 		& x, 
								const typename IMultiplyVariables::ElementType 		& y,
								const 	IMultiplyCoeffWithVariable 			& imultCoeffXVariable,
								const 	IMultiplyVariables 				& imultVariableXVariable
							 ) const
{

	//assert()

	typename IMultiplyVariables::MultiplicationResultType	res=	IMultiplyVariables::ElementType::One;

	for (int currXexp=1;	currXexp<=(int)exponents.getXExp(); currXexp++)
	{
		res	=	imultVariableXVariable.multiply(res,x);
	}

	for (int currYexp=1;	currYexp<=(int)exponents.getYExp(); currYexp++)
	{
		res	=	imultVariableXVariable.multiply(res,y);
	}

	for (int currXexp=-1;	currXexp>=exponents.getXExp(); currXexp--)
	{
		res	=	imultVariableXVariable.multiply(res, imultVariableXVariable.multInv(x) );
	}

	for (int currYexp=-1;	currYexp>=exponents.getYExp(); currYexp--)
	{
		res	=	imultVariableXVariable.multiply(res,imultVariableXVariable.multInv(y) );
	}

	return imultCoeffXVariable.scalarMultiply(coeff, res);
};


///
/// es ist nur der Sonderfall implementiert, in welchem coeff*x^i*y^j wieder in IMultiplyRing::ElementType landet. Dies ist nicht immer der Fall.
template <class ccoeff>
template <class IMultiplyRing  >
typename IMultiplyRing::ElementType 	xyMonom<ccoeff>::substitute2(			
								const typename IMultiplyRing::ElementType 		& x, 
								const typename IMultiplyRing::ElementType 		& y,
								const 	IMultiplyRing 				& imult
							 ) const
{

	
	/*typename IMultiplyRing::ElementType	res=IMultiplyRing::ElementType::One;

	res.setEpsPrecision( imult.getEpsPrecision() );*/

	typename IMultiplyRing::ElementType	res(imult.getEpsPrecision() ,std::string("") );
	

	assert(imult.getEpsPrecision()>= x.getEpsPrecision() );
	assert(imult.getEpsPrecision()>= y.getEpsPrecision() );
	res.setValue(0,1);
	
	for (int currXexp=1;	currXexp<=(int)exponents.getXExp(); currXexp++)
	{
		imult.multiplyInPlaceRef(res,x);
	}

	for (int currYexp=1;	currYexp<=(int)exponents.getYExp(); currYexp++)
	{
		 	imult.multiplyInPlaceRef(res,y);
			
	}
	
	assert(exponents.getYExp()>=0);
	assert(exponents.getYExp()>=0);
	
	///  @todo merke: was negatives ist meistens größer als ein unsigned int!
	/*for (int currXexp=-1;	currXexp>=(int)exponents.getXExp(); currXexp--)
	{
		 	imult.multiplyInPlace(res,  imult.multInv(x) );
	}

	for (int currYexp=-1;	currYexp>=(int)exponents.getYExp(); currYexp--)
	{
		 	imult.multiplyInPlace(res, imult.multInv(y) );
	}*/

	return imult.scalarMultiply(coeff, res);
	///res= imult.multiply(coeff, res);
	//res= imult.multiply( res, coeff);
	//return res;
};


template <class xyMonomType>
class xyOneFormTerm
{
	public:
		 enum TermType
	  	 {
      		DXTERM, DYTERM
		} ;

		xyOneFormTerm(const xyMonomType & mon, TermType whichForm);
		xyOneFormTerm(const std::string & str);
		xyOneFormTerm(std::stringstream&  monomListStream);
	
		xyMonomType	getMonom() const	{	return xyMonom_m;	}
		int getDegree()	const	{	return xyMonom_m.getDegree();	}
		bool	isDxTerm() const	{	return whichDForm_m==xyOneFormTerm::DXTERM;	};
		bool	isDyTerm() const	{	return whichDForm_m==xyOneFormTerm::DYTERM;	};

	private:
		 TermType 		whichDForm_m;
		xyMonomType		xyMonom_m;
};


template <class xyMonomType>
xyOneFormTerm<xyMonomType>::xyOneFormTerm(	const xyMonomType & mon,
					 TermType whichForm	):	xyMonom_m(mon), 
										whichDForm_m(whichForm)
{
	
}

/// @todo '+-' sollte nicht erlaubt sein,
/// @todo beim Einlesen der Differentialform kann diese theoretisch = 0 sein - dann geht das Einlesen wohl schief.
template <class xyMonomType>
xyOneFormTerm<xyMonomType>::xyOneFormTerm(const std::string & paramstr)
{

	std::string error("Error reading xyOneFormTerm from stream: xyOneFormTerm only supports 'xyMonom'*dx and 'xyMonom'*dy -formed terms !\n");
	
//	str.erase(std::remove_if(str.begin(), str.end(), std::isspace), str.end() );

	std::string str=eatWS(paramstr);

	size_t pos = str.find("d");

	#ifdef DEBUG
	std::cerr << "paramstr = " << paramstr << std::endl;
	#endif
	// '-' und '+' sind vorerst nur als führendes Zeichen erlaubt! - TODO: Ja und was machst du, wenn es mehrere '-' und '+' gibt?
	size_t posPlus = str.find("+");
	if (posPlus!=str.npos)
	{
		assert(posPlus==0);
	}

	posPlus = str.find("-");
	if (posPlus!=str.npos)
	{
		assert(posPlus==0);
	}

	if (pos==str.npos)
	{
		std::cerr << error<<	std::endl;;
		assert (pos!=str.npos);
	}
	
	str.substr(pos);

	/// also: erstmal alle Whitespaces wegräumen.
	/// wenn nur dx auftaucht, dann ist pos=0 und es handelt sich um  ein  1-Monom.
	if (pos==0 && str.length()>0 )
	{
		xyMonom_m = xyMonomType(1, 0, 0 );
	}
	else
	{

		assert(pos>0);	// wenn ein Monom am Anfang steht, gibt es entweder ein '...*'dy  oder '-dy' !.
		xyMonom_m = xyMonomType(1, 0, 0 );
		std::string monomString ;
		if (pos==1)
		{
			
			monomString= str.substr(0, pos);
			assert(monomString.at(0)=='-' || monomString.at(0)=='+');
			if ( monomString.at(0)=='+')
				monomString= monomString.substr(1);
		}
		else
		{
			assert(str.at(pos-1)=='*');
			monomString= str.substr(0, pos-1);		
			
		}
		//std::cerr << "monomString" << monomString <<std::endl;
		std::stringstream sstream(monomString);
		if (monomString.length()==0)
			xyMonom_m = xyMonomType(1, 0, 0 );
		else
			xyMonom_m = xyMonomType(sstream);
	
	};
	
	//std::string differentialFormString 	std.substring(pos)
	std::stringstream 	differentialFormString(	str.substr(pos) );
	//std::cerr << "differentialFormString = '" << differentialFormString.str() << "'" << std::endl;
	extractChar( 'd', differentialFormString );
	


	if (differentialFormString.eof() )
	{
		std::cerr << error << std::endl;;
		assert( ! differentialFormString.eof() );
	}

	if (differentialFormString.peek()=='x')
	{
		extractChar('x', differentialFormString);
		whichDForm_m=DXTERM;
	}
	else if (differentialFormString.peek()=='y')
	{
		extractChar('y', differentialFormString);
		whichDForm_m=DYTERM;
	}
	else
	{
		std::cerr << error ;
		assert(false);
	}
	differentialFormString.peek();
	if ( ! differentialFormString.eof() )
	{
		std::cerr << error << std::endl;
		std::cerr << "differentialFormString" << differentialFormString.str() << std::endl;
		assert(  differentialFormString.eof() );
	}
}

//string::npos


/*
template <class xyMonomType>
xyOneFormTerm<xyMonomType>::xyOneFormTermNew(const std::string & paramstr)
{
	//1. Zerlege String in Einzelteile.
}
*/



template <class PolynomXYType>
class xyOneForm 
{
	public:

		xyOneForm(const PolynomXYType & pdx, const PolynomXYType & pdy);
		xyOneForm(const std::string  & str);
		xyOneForm(std::stringstream&  monomListStream);
	
		PolynomXYType	getDxFormPart() const	{	return polynomDX_m;	}
		PolynomXYType	getDyFormPart() const	{	return polynomDY_m;	}


	private:

		 PolynomXYType 		polynomDX_m;
		PolynomXYType 		polynomDY_m;
		
};


template <class PolynomXYType>
xyOneForm<PolynomXYType>::xyOneForm(	const PolynomXYType & pdx,
				 const PolynomXYType & pdy): 	polynomDX_m(pdx),
									polynomDY_m(pdy)
{

}

/// TODO: kein Term darf doppelt vorkommen. 
template <class PolynomXYType>
xyOneForm<PolynomXYType>::xyOneForm(const std::string  & str)
{

	std::string tmpString = eatWS(str);
	// 1. zerlege Zeichenkette in Summanden und parse jede einzelnen Summanden.
	typedef xyMonom < typename PolynomXYType::CoefficientType> xyMonomType;

	std::vector<xyOneFormTerm <xyMonomType > > monomTerms;

	size_t pos = tmpString.npos;

	// todo: zugelassene Tokens sind nur 'dx','dy' 'x','y','x^exponent','y^exponent','+','-', '(0...9)*'. , Exponent ist eine Zahl ohne Vorzeichen. Keine Klammern etc.!
	size_t posKlammer = tmpString.find_last_of('(');
	assert(   posKlammer==tmpString.npos);
	posKlammer = tmpString.find_last_of(')');

	
	assert(   posKlammer==tmpString.npos);

	//suche von Ende der Zeichenkette aus nach Summanden  und trage diese in die monomTerms-Liste ein. TODO: ->Fehleranfaellig, da Klammerung vorkommen kann, etc. 
	
	while (pos !=0)
	{
		size_t posPlus = tmpString.find_last_of('+');
		size_t posMinus = tmpString.find_last_of('-');
	
		assert( posPlus==tmpString.npos || posPlus<tmpString.npos);
		assert(posMinus==tmpString.npos || posMinus<tmpString.npos);
	
		if (posPlus==tmpString.npos)
			pos=posMinus;
		else if (posMinus==tmpString.npos)
			pos=posPlus;
		else
			pos= std::max(posPlus, posMinus);

		if (pos==tmpString.npos)
			pos=0;

		// wenn das erste Zeichen ein Plus(+) oder Minus(-) ist, suche bis zum naechsten Zeichen
		// alternativ: 
		assert( tmpString.npos != 0 );
		if ( pos != tmpString.npos )
		{
			std::string monomString = tmpString.substr(pos);
			xyOneFormTerm <xyMonomType >  oneFormTerm(monomString);
			monomTerms.push_back(oneFormTerm);
			tmpString= tmpString.substr(0,pos);
		}
	}
	
	// 2. ermittle höchsten vorkommenden Grad der Monome

	int maxDegree=0;

	for (size_t i=0;i < monomTerms.size(); i++)
	{
		maxDegree = max( monomTerms[i].getDegree(), maxDegree);
	}

	// 3. lege polynomDX_m und polynomDY_m an - fertig!


	polynomDX_m.clear();
	polynomDX_m.setDegree(maxDegree);
	polynomDY_m.clear();
	polynomDY_m.setDegree(maxDegree);

	
	PolynomXYType pObserver(maxDegree);
	PolynomXYType qObserver(maxDegree);


	for (size_t i=0;i < monomTerms.size(); i++)
	{

		xyMonomType mon = monomTerms[i].getMonom();
		if (monomTerms[i].isDxTerm() )
		{
			//sicherstellen, dass kein Koeffizient in der Eingabe doppelt vorkommt
			assert( pObserver.getCoeff(mon.getXExp(), mon.getYExp())==typename PolynomXYType::CoefficientType(0) );
			polynomDX_m.setCoeff(mon.getXExp(), mon.getYExp(), mon.getCoeff() );
			pObserver.setCoeff(mon.getXExp(), mon.getYExp(), 1 );
		}
		else if (monomTerms[i].isDyTerm() )
		{
			//sicherstellen, dass kein Koeffizient in der Eingabe doppelt vorkommt 
			assert( qObserver.getCoeff(mon.getXExp(), mon.getYExp())==typename PolynomXYType::CoefficientType(0) );
			polynomDY_m.setCoeff( mon.getXExp(), mon.getYExp(), mon.getCoeff() );
			qObserver.setCoeff(mon.getXExp(), mon.getYExp(), 1 );
		}
		else
			assert(false);

	}

	
}

