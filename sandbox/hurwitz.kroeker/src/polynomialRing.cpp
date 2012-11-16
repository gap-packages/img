
/// @todo hier kann man ganz viel code sparen...

#include <algorithm>
#include <iostream>
#include <sstream>

template<class TPolynomXY, class TRing>
PolynomialRing<TPolynomXY, TRing>::PolynomialRing(const TRing & ring): ring_ref_m(ring)
	{
	}


/// @pre der Ring ring 1 hat die Elemente vom Typ  PolynomXY_Type::CoefficientType.
template<class TPolynomXY, class TRing>
TPolynomXY 		PolynomialRing<TPolynomXY, TRing>::addInv	 (const	TPolynomXY 	& _polynom_ref) const
{

	TPolynomXY  res(_polynom_ref.getDegree() );
	
	for (short  i=0; i<= _polynom_ref.getDegree(); i++)
		{
			for(short j=0; j<=i; j++)
			{
				res.setCoeff (i-j, j , ring_ref_m.addInv( _polynom_ref.getCoeffConst(i-j, j) ) );
				#ifdef COUNT
					addCount+=1;
				#endif
			}
		}
	return res;
}




/// @todo das Durchlaufen der Schleifen wiederholt sich, das kann man doch irgendwie abstrahieren!
template<class TPolynomXY, class TRing>
TPolynomXY*		PolynomialRing<TPolynomXY, TRing>::addInvReturnPtr(const	TPolynomXY 	& _polynom_ref	) const
{
	TPolynomXY * res= new TPolynomXY(_polynom_ref.getDegree() );
	
	for (short  i=0; i<= _polynom_ref->getDegree(); i++)
		{
			for(short j=0; j<=i; j++)
			{
				res->setCoeff (i-j, j , ring_ref_m.addInv( _polynom_ref.getCoeffConst(i-j, j) ) );
				#ifdef COUNT
					addCount+=1;
				#endif
			}
		}
	return res;
}


template<class TPolynomXY, class TRing>
	void 			PolynomialRing<TPolynomXY, TRing>::addInvInPlace(	TPolynomXY 	& _polynom_ref	) const
{

	for (short  i=0; i<= _polynom_ref->getDegree(); i++)
		{
			for(short j=0; j<=i; j++)
			{
				ring_ref_m.addInvInPlace( _polynom_ref.getCoeffRef(i-j, j) );
				#ifdef COUNT
					addCount+=1;
				#endif
			}
		}
	return ;
}


template<class TPolynomXY, class TRing>
TPolynomXY 		PolynomialRing<TPolynomXY, TRing>::add(const TPolynomXY & _polynom_1_ref,
										const TPolynomXY & _polynom_2_ref) const
{

	// add polynoms for different Degree currently not implemented.
	assert(_polynom_1_ref.getDegree()==_polynom_2_ref.getDegree() );

    std::cerr << "_polynom_1_ref.getDegree() " << _polynom_1_ref.getDegree()<< std::endl;
    std::cerr << "_polynom_1_ref.getDegree() " << _polynom_2_ref.getDegree()<< std::endl;

	const short resDegree= max(_polynom_1_ref.getDegree(),_polynom_2_ref.getDegree() );
	TPolynomXY  res (resDegree );
	//i-j is the x-exp and j is the y-exp., i is the currentMonomDegree
	for (short  currMonomDegree=0; currMonomDegree<= resDegree; currMonomDegree++)
		{
			for(short yExp=0; yExp<=currMonomDegree; yExp++)
			{
				short xExp = currMonomDegree - yExp;

				res.setCoeff (xExp, yExp , ring_ref_m.add(	_polynom_1_ref.getCoeffConst(xExp, yExp ),
										_polynom_2_ref.getCoeffConst(xExp, yExp ) 
								)
							);
				#ifdef COUNT
					addCount+=1;
				#endif
			}
		}
	return res;
}

template<class TPolynomXY, class TRing>
TPolynomXY* 	PolynomialRing<TPolynomXY, TRing>::addReturnPtr( const	TPolynomXY 	 & _polynom_1_ref,
											const	TPolynomXY 	 & _polynom_2_ref) const
{
// add polynoms for different Degree currently not implemented.
	assert(_polynom_1_ref.getDegree()==_polynom_2_ref.getDegree() );

	const short resDegree= max(_polynom_1_ref.getDegree(),_polynom_2_ref.getDegree);
	TPolynomXY *  res = new TPolynomXY(resDegree );
	//i-j is the x-exp and j is the y-exp., i is the currentMonomDegree
	for (short  currMonomDegree=0; currMonomDegree<= resDegree; currMonomDegree++)
		{
			for(short yExp=0; yExp<=currMonomDegree; yExp++)
			{
				short xExp = currMonomDegree - yExp;

				res->setCoeff (xExp, yExp , ring_ref_m.add(	_polynom_1_ref.getCoeffConst(xExp, yExp ),
										_polynom_2_ref.getCoeffConst(xExp, yExp ) 
								)
							);
				#ifdef COUNT
					addCount+=1;
				#endif
			}
		}
	return res;
}


template<class TPolynomXY, class TRing>
inline void 			PolynomialRing<TPolynomXY, TRing>::addInPlace( 	TPolynomXY 			& _polynom_1_ref,
											const TPolynomXY 		& _polynom_2_ref) const
{
	// add polynoms for different Degree currently not implemented.
	assert(_polynom_1_ref.getDegree() == _polynom_2_ref.getDegree() );

	//i-j is the x-exp and j is the y-exp., i is the currentMonomDegree
	for (short  currMonomDegree=0; currMonomDegree<= _polynom_1_ref.getDegree(); currMonomDegree++)
		{
			for(short yExp=0; yExp<=currMonomDegree; yExp++)
			{
				short xExp = currMonomDegree - yExp;

					ring_ref_m.addInPlace(	_polynom_1_ref.getCoeffRef(xExp, yExp ),
							_polynom_2_ref.getCoeffConst(xExp, yExp ) 
							);
				#ifdef COUNT
					addCount+=1;
				#endif
			}
		}
	return ;
}


template<class TPolynomXY, class TRing>
TPolynomXY 		PolynomialRing<TPolynomXY, TRing>::scalarMultiply(const typename TPolynomXY::CoefficientType 	scalar,
										 const	TPolynomXY 			& _polynom_ref) const
{
	TPolynomXY  res( _polynom_ref.getDegree() );
	//i-j is the x-exp and j is the y-exp., i is the currentMonomDegree
	for (short  currMonomDegree=0; currMonomDegree<= _polynom_ref.getDegree(); currMonomDegree++)
		{
			for(short yExp=0; yExp <= currMonomDegree; yExp++)
			{
				short xExp = currMonomDegree - yExp;

				res.setCoeff (xExp, yExp , ring_ref_m.scalarMultiply(	scalar,
											_polynom_ref.getCoeffConst( xExp, yExp ) 
							)
					);
				#ifdef COUNT
					addCount+=1;
				#endif
			}
		}
 	//short polynomDegree = _polynom_ref.getDegree()
	//short xexp,yexp, currDegree;
	// begin(xexp,yexp,currDegree, polynomDegree);
	// while (next (xexp,yexp,polynomDegree))
	//for (begin(xexp,yexp,polynomDegree); !end(currDegree, polynomDegree);  next (xexp, yexp, currDegree, polynomDegree))
	//	res.setCoeff (xExp, yExp , ring_ref_m.scalarMultiply(	scalar,
	//										_polynom_ref.getCoeffConst( xExp, yExp ) 
	//						)
	//				);
	return res;
}


/*
// so ungefähr könnte man es mit dem Durchlaufen der Monome machen: begin() und next()
// da brauchst du aber strenggenommen eine const-iterator, einen Ref-Iterator und einen normalen Iterator
inline bool begin(short & xexp, short & yexp , short currDegree, int polynomialDegree)
{ 
	currDegree = 0;
	xexp = 0;
	yexp = 0;
	return ( polynomialDegree>=0 );
}

inline bool end(  short & currDegree, int polynomialDegree)
{ 
	return  (currDegree<=polynomialDegree);
}

// so ungefähr könnte man ex machen.
inline bool nextMonomExp(short & xexp, short & yexp, short &currDegree, int polynomialDegree)
{
	//short currDegree = xexp + yexp;
	if (yexp<currDegree)
	{
		yexp++;
		xexp--;
	}
	else 
	{
		currDegree++;
		xexp=currDegree;
		yexp=0;
		if (currDegree>polynomialDegree)
			return false;
	}
	return true;

}*/

template<class TPolynomXY, class TRing>
void 		PolynomialRing<TPolynomXY, TRing>::scalarMultiplyInPlace(const typename TPolynomXY::CoefficientType 	scalar,
										TPolynomXY 				& _polynom_ref	) const
{
 
	//i-j is the x-exp and j is the y-exp., i is the currentMonomDegree
	for (short  currMonomDegree=0; currMonomDegree<= _polynom_ref.getDegree(); currMonomDegree++)
		{
			for(short yExp=0; yExp <= currMonomDegree; yExp++)
			{
				short xExp = currMonomDegree - yExp;

				/*_polynom_ref.setCoeff (xExp, yExp , ring_ref_m.scalarMultiply(	scalar,
												_polynom_ref.getCoeffConst( xExp, yExp ) 
							)
					);*/
					 ring_ref_m.scalarMultiplyInPlace(scalar, _polynom_ref.getCoeffRef (xExp, yExp) );

				#ifdef COUNT
					addCount+=1;
				#endif
			}
		}
}

template<class TPolynomXY, class TRing>
TPolynomXY*	PolynomialRing<TPolynomXY, TRing>::scalarMultiplyRetPtr(const	typename TPolynomXY::CoefficientType 	scalar,
										const	TPolynomXY 			& _polynom_ref) const
{
	
	TPolynomXY*  res = new TPolynomXY( _polynom_ref.getDegree() );
	//i-j is the x-exp and j is the y-exp., i is the currentMonomDegree
	for (short  currMonomDegree=0; currMonomDegree<= _polynom_ref.getDegree(); currMonomDegree++)
		{
			for(short yExp=0; yExp <= currMonomDegree; yExp++)
			{
				short xExp = currMonomDegree - yExp;

				res->setCoeff (xExp, yExp , ring_ref_m.scalarMultiply(	scalar,
											_polynom_ref.getCoeffConst( xExp, yExp ) 
							)
					);
				#ifdef COUNT
					addCount+=1;
				#endif
			}
		}
 	//short polynomDegree = _polynom_ref.getDegree()
	//short xexp,yexp, currDegree;
	// begin(xexp,yexp,currDegree, polynomDegree);
	// while (next (xexp,yexp,polynomDegree))
	//for (begin(xexp,yexp,polynomDegree); !end(currDegree, polynomDegree);  next (xexp, yexp, currDegree, polynomDegree))
	//	res.setCoeff (xExp, yExp , ring_ref_m.scalarMultiply(	scalar,
	//										_polynom_ref.getCoeffConst( xExp, yExp ) 
	//						)
	//				);
	return res;
}

/// @pre: srcPol!=destPol
template<class TPolynomXY, class TRing>
template <class TPolynomXY_SRC_Type, class TPolynomXY_DEST_Type>
void 	PolynomialRing<TPolynomXY, TRing>::copyPolynomWithGivenEpsPrecision(const TPolynomXY_SRC_Type  & srcPol,
									 TPolynomXY_DEST_Type  	&	 destPol, int epsPrecision)
{	
	//assert((&srcPol) != (&destPol)) ;

	for (int deg = srcPol.getDegree(); deg>=0; deg--)
	{
		for (int x_exp=deg; x_exp>=0; x_exp--)
		{
			int MaxEpsPrecision= min(epsPrecision, (int)srcPol.getCoeff(x_exp, deg-x_exp).getEpsPrecision() );

			for (int prec=0; prec<= MaxEpsPrecision ; prec++)
			{
				destPol.getCoeffRef(	x_exp, deg-x_exp).setValue(prec, 0);
			}

			
			for (int prec=0; prec<= MaxEpsPrecision ; prec++)
			{
				destPol.getCoeffRef(	x_exp, deg-x_exp).setValue(prec,
											srcPol.getCoeff(x_exp, deg-x_exp).getValue(prec) );
			}
		}
	}
}



template<class TPolynomXY, class TRing>
template <class TPolynomXY_Type>
void 	PolynomialRing<TPolynomXY, TRing>::convertInPlace( TPolynomXY_Type  & pol ) const
{
	for (int deg = pol.getDegree(); deg>=0; deg--)
	{
		for (int x_exp=deg; x_exp>=0; x_exp--)
		{
			int MaxEpsPrecision= (int)pol.getCoeff(x_exp, deg-x_exp).getEpsPrecision();
			for (int prec=0; prec<= MaxEpsPrecision ; prec++)
			{
				pol.getCoeffRef(	x_exp, deg-x_exp).setValue(prec,
											ring_ref_m.ConvertScalar(pol.getCoeff(x_exp, deg-x_exp).getValue(prec)) );
			}
		}
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////
template<class TPolynomX, class TRing>
UnivariatePolynomialRing<TPolynomX, TRing>::UnivariatePolynomialRing(const TRing & ring): ring_ref_m(ring)
	{
	}

template<class TPolynomX, class TRing>
TPolynomX 		UnivariatePolynomialRing<TPolynomX, TRing>::scalarMultiply(const typename TPolynomX::CoefficientType 	scalar,
											 const	TPolynomX 			& _polynom_ref) const
{
	int polynomialDegree =_polynom_ref.getDegree();

	TPolynomX  res( polynomialDegree  );



	for(short xExp=0; xExp <= polynomialDegree; xExp++)
	{
				res.setCoeff (xExp,  ring_ref_m.scalarMultiply(	scalar,	_polynom_ref.getCoeffConst( xExp ) 
							)
					);
	}

	return res;

}



template<class TPolynomX, class TRing>
TPolynomX 		UnivariatePolynomialRing<TPolynomX, TRing>::addInv( const	TPolynomX 			& _polynom_ref) const
{
	int polynomialDegree =_polynom_ref.getDegree();

	TPolynomX  res( polynomialDegree  );


	for(short xExp=0; xExp <= polynomialDegree; xExp++)
	{
				res.setCoeff (xExp,  ring_ref_m.addInv(	_polynom_ref.getCoeffConst( xExp ) 	)	);
	}

	return res;

}

template<class TPolynomX, class TRing>
TPolynomX 		UnivariatePolynomialRing<TPolynomX, TRing>::add( const	TPolynomX 			& _polynom_ref1,  
										const	TPolynomX 			& _polynom_ref2) const
{
	int polynomialDegree = std::max(_polynom_ref1.getDegree(),_polynom_ref2.getDegree() );

	TPolynomX  res( polynomialDegree  );


	for(short xExp=0; xExp <= polynomialDegree; xExp++)
	{
				res.setCoeff (xExp,  ring_ref_m.add(	_polynom_ref1.getSafeCoeffConst( xExp ) 	 ,
											_polynom_ref2.getSafeCoeffConst( xExp ) 	
									)
						);
	}

	return res;

}



template<class TPolynomX, class TRing>
TPolynomX 		UnivariatePolynomialRing<TPolynomX, TRing>::multiply( const	TPolynomX 			& _polynom_ref1,  
										const	TPolynomX 			& _polynom_ref2) const
{
	short int pol_1_deg= _polynom_ref1.getDegree();
	short int pol_2_deg= _polynom_ref2.getDegree();

	short int polynomialDegree =   pol_1_deg + pol_2_deg ;

	TPolynomX  res( polynomialDegree  );


	for(short xExp=0; xExp <= polynomialDegree; xExp++)
	{
		short int xfinish = std::min( xExp,    pol_1_deg     ); 
		short int xstart  = std::max( 0,    xExp - pol_2_deg );

		for ( short x=xfinish; x>= xstart; x--)
		{

			if (x <=pol_1_deg && (xExp-x) <=  pol_2_deg )
			{
				ring_ref_m.addInPlace(	res.getCoeffRef(xExp), 
											ring_ref_m.multiply( _polynom_ref1.getCoeffConstRef(x),
														 _polynom_ref2.getCoeffConstRef(xExp-x)
											
										)
							);
			}
		}
	}

	return res;
}



template<class TPolynomX, class TRing>
TPolynomX* 		UnivariatePolynomialRing<TPolynomX, TRing>::multiplyPtr( const	TPolynomX 			& _polynom_ref1,  
										const	TPolynomX 			& _polynom_ref2) const
{
	short int pol_1_deg= _polynom_ref1.getDegree();
	short int pol_2_deg= _polynom_ref2.getDegree();

	short int polynomialDegree =   pol_1_deg + pol_2_deg ;

	TPolynomX*   resPtr = new TPolynomX( polynomialDegree  );
	TPolynomX&   res(*resPtr);

	for(short xExp=0; xExp <= polynomialDegree; xExp++)
	{
		short int xfinish = std::min( xExp,    pol_1_deg     ); 
		short int xstart  = std::max( 0,    xExp - pol_2_deg );

		for ( short x=xfinish; x>= xstart; x--)
		{

			if (x <=pol_1_deg && (xExp-x) <=  pol_2_deg )
			{
				ring_ref_m.addInPlace(	res.getCoeffRef(xExp), 
											ring_ref_m.multiply( _polynom_ref1.getCoeffConstRef(x),
														 _polynom_ref2.getCoeffConstRef(xExp-x)
											
										)
							);
			}
		}
	}

	return resPtr;
}

template<class TPolynomX, class TRing>
    void    UnivariatePolynomialRing<TPolynomX, TRing>::normalizeInPlace(TPolynomX          & _polynom_ref1) const
{
    int degree=_polynom_ref1.getExactDegree();
    typename TRing::ElementType    correctionCoeff = ring_ref_m.multInv( _polynom_ref1.getCoeff(degree) ) ;
    for ( int currExponent= degree; currExponent>=0; currExponent--)
          ring_ref_m.multiplyInPlace( _polynom_ref1.getCoeffRef(currExponent), correctionCoeff );
    return;
}


template<class TPolynomX, class TRing>
TPolynomX 		UnivariatePolynomialRing<TPolynomX, TRing>::square( const	TPolynomX 			& _polynom_ref1 ) const
{
	short int pol_1_deg= _polynom_ref1.getDegree();

	short int polynomialDegree =   pol_1_deg + pol_1_deg ;

	TPolynomX  res( polynomialDegree  );


	for(short xExp=0; xExp <= polynomialDegree; xExp++)
	{
		short int xfinish = std::min( xExp,    pol_1_deg     ); 
		short int xstart  = std::max( 0,    xExp - pol_1_deg );

		for ( short x=xfinish; x>= xstart; x--)
		{

			if (x <=pol_1_deg && (xExp-x) <=  pol_1_deg )
			{
				ring_ref_m.addInPlace(	res.getCoeffRef(xExp), 
											ring_ref_m.multiply( _polynom_ref1.getCoeffConstRef(x),
														 _polynom_ref1.getCoeffConstRef(xExp-x)
											
										)
							);
			}
		}
	}

	return res;
}
 


template<class TPolynomX, class TRing>
pair<TPolynomX,TPolynomX> 		UnivariatePolynomialRing<TPolynomX, TRing>::divide( const	TPolynomX 			& dividend,  
										const	TPolynomX 			& divisor) const
{
/*
	FOR i = GradZ - GradN TO 0 STEP -1
		Quotient(i) = Zähler(i + GradN) / Nenner(GradN)
		FOR j = GradN TO 0 STEP -1
			Zähler(i + j) = Zähler(i + j) - Nenner(j) * Quotient(i)
		NEXT j
	NEXT i
	FOR j = GradN - 1 TO 0 STEP -1
    		Rest(j) = Zähler(j)
	NEXT j
	*/

	assert(! divisor.isZero() );
	TPolynomX dividend_copy(dividend);

	int dividend_deg= dividend.getExactDegree();
	int divisor_deg = divisor.getExactDegree();

    if (dividend_deg - divisor_deg <0)
    {
          //res(0 );
        return pair<TPolynomX,TPolynomX>( TPolynomX(0), dividend );
    }
    
	TPolynomX  res( dividend_deg - divisor_deg  );

	for (int i = dividend_deg - divisor_deg; i>=0; i--)
	{
		res.setCoeff(i, ring_ref_m.multiply( dividend_copy.getCoeffConstRef(i + divisor_deg),  
								  ring_ref_m.multInv(divisor.getCoeffConstRef( divisor_deg) ) 	
								)
				);

		for (int j=divisor_deg; j>=0; j--)
			 ring_ref_m.addInPlace(	dividend_copy.getCoeffRef( i + j ),
							ring_ref_m.addInv(ring_ref_m.multiply( res.getCoeffConst(i), divisor.getCoeffConst(j) ) )
						);
		
	}
	TPolynomX remainder(dividend_copy.getExactDegree() );

	for (int i = dividend_copy.getExactDegree(); i>=0; i--)
	{
		remainder.setCoeff(i, dividend_copy.getCoeffConst(i) );
	}

	return pair<TPolynomX,TPolynomX>(res, remainder);
}


template<class TPolynomX, class TRing>
TPolynomX		UnivariatePolynomialRing<TPolynomX, TRing>::remainder( const	TPolynomX 			& dividend,  
										const	TPolynomX 			& divisor) const
{
/*
	FOR i = GradZ - GradN TO 0 STEP -1
		Quotient(i) = Zähler(i + GradN) / Nenner(GradN)
		FOR j = GradN TO 0 STEP -1
			Zähler(i + j) = Zähler(i + j) - Nenner(j) * Quotient(i)
		NEXT j
	NEXT i
	FOR j = GradN - 1 TO 0 STEP -1
    		Rest(j) = Zähler(j)
	NEXT j
	*/

	assert(! divisor.isZero() );
	TPolynomX dividend_copy(dividend);

	int dividend_deg= dividend.getExactDegree();
	int divisor_deg = divisor.getExactDegree();

	for (int i = dividend_deg - divisor_deg; i>=0; i--)
	{
		typename TRing::ElementType el =ring_ref_m.multiply( dividend_copy.getCoeffConstRef(i + divisor_deg),  
								  ring_ref_m.multInv(divisor.getCoeffConstRef( divisor_deg) ) 	
								);

		for (int j=divisor_deg; j>=0; j--)
			 ring_ref_m.addInPlace(	dividend_copy.getCoeffRef( i + j ),
							ring_ref_m.addInv(ring_ref_m.multiply( el, divisor.getCoeffConst(j) ) )
						);
		
	}
	TPolynomX remainder(dividend_copy.getExactDegree() );

	for (int i = dividend_copy.getExactDegree(); i>=0; i--)
	{
		remainder.setCoeff(i, dividend_copy.getCoeffConst(i) );
	}

	return remainder;
}

template<class TPolynomX, class TRing>
TPolynomX 		UnivariatePolynomialRing<TPolynomX, TRing>::pow(	const TPolynomX 	& polynom,
										unsigned int 	exponent  ) const
{
	
	if (exponent==0)
	{
		TPolynomX	res(0);
		res.setCoeff(0,TRing::ElementType::One);
		return res;
	}
	TPolynomX	res(polynom);
	while (exponent>1)
	{
		res=multiply(res,polynom);
		exponent--;
	}
	return res;
}


template<class TPolynomX, class TRing>
TPolynomX		UnivariatePolynomialRing<TPolynomX, TRing>::gcd( const	TPolynomX 			& pol1_ref,  
										const	TPolynomX 			& pol2_ref) const
{
	TPolynomX 			 pol1 (  pol1_ref);  
	TPolynomX 			 pol2 (  pol2_ref);

	TPolynomX 			* pol1_p =  & pol1;  
	TPolynomX 			* pol2_p =  & pol2;

	if ( pol1_ref.getExactDegree() < pol2_ref.getExactDegree() )
	{
		pol1_p = & pol2;
		pol2_p = & pol1;
	}	

	assert(! (*pol2_p).isZero() );

	/*	rk−2(x) = qk(x) rk−1(x) + rk(x) */

	// Plan: zwei Polynome , und zwei Zeiger 'prev' und 'curr' die nach jedem Schritt wecseln. 
	//std::cerr << " *pol1_p: " << *pol1_p << std::endl;
	//std::cerr << " *pol2_p: " << *pol2_p << std::endl;
	TPolynomX	prev=*pol2_p;
	pair<TPolynomX,TPolynomX> divres = divide(*pol1_p,*pol2_p);

	while (! divres.second.isZero() )
	{
 		*pol1_p = prev;
		//pol2_p = &(divres.second);
		prev=divres.second;
		pol2_p= &prev;
		//std::cerr << "div result: " << divres.first << std::endl;
		//std::cerr << "div remainder: " << *pol2_p << std::endl;
		//assert(! (*pol2_p).isZero() );
		divres = divide( *pol1_p, *pol2_p );
		//std::cerr << "2. div result: " << divres.first << std::endl;
		//std::cerr << "2 div remainder: " <<divres.second << std::endl;
		//std::cerr << "*pol2_p: " <<*pol2_p << std::endl;
		
	}
	//std::cerr << " *pol1_p: " << *pol1_p << std::endl;
	//std::cerr << " *pol2_p: " << *pol2_p << std::endl;

	return prev;

}


template<class TPolynomX, class TRing>
TPolynomX		UnivariatePolynomialRing<TPolynomX, TRing>::fastgcd( const	TPolynomX 			& pol1_ref,  
										const	TPolynomX 			& pol2_ref) const
{
	TPolynomX 			 pol1 (  pol1_ref);  
	TPolynomX 			 pol2 (  pol2_ref);

	TPolynomX 			* pol1_p =  & pol1;  
	TPolynomX 			* pol2_p =  & pol2;

	if ( pol1_ref.getExactDegree() < pol2_ref.getExactDegree() )
	{
		pol1_p = & pol2;
		pol2_p = & pol1;
	}	

	assert(! (*pol2_p).isZero() );

	/*	rk−2(x) = qk(x) rk−1(x) + rk(x) */

	// Plan: zwei Polynome , und zwei Zeiger 'prev' und 'curr' die nach jedem Schritt wecseln. 
	//std::cerr << " *pol1_p: " << *pol1_p << std::endl;
	//std::cerr << " *pol2_p: " << *pol2_p << std::endl;
	TPolynomX	prev=*pol2_p;
	TPolynomX divrem = remainder(*pol1_p,*pol2_p);

	while (! divrem.isZero() )
	{
 		*pol1_p = prev;
		prev=divrem;
		divrem = remainder( *pol1_p, prev );
	}
	

	return prev;
}

template<class TPolynomX, class TRing>
bool		UnivariatePolynomialRing<TPolynomX, TRing>::fastgcdspec( TPolynomX 			& pol1_ref,  
												const	TPolynomX 		& pol2_ref,
												int 				minDeg) const
{

	//assert ( pol1_ref.getExactDegree() >= pol2_ref.getExactDegree() );
	assert(! (pol2_ref).isZero() );

	/*	rk−2(x) = qk(x) rk−1(x) + rk(x) */
	TPolynomX		rkm1= pol2_ref ;
	TPolynomX		rkm2= pol1_ref ;

	TPolynomX	*	rkm2_p= &rkm2;
	TPolynomX	*	rkm1_p= &rkm1 ;


	if ( pol1_ref.getExactDegree() < pol2_ref.getExactDegree() )
	{
		rkm2_p = & rkm1;
		rkm1_p = & rkm2;
	}	
 
	// um Kopieren zu vermeiden, brauchst du doch drei Polynome..damit diese nicht ständig erzeugt werden, sollten diese mit übergeben werden.
	TPolynomX	*	rk_p =  rkm2_p ;
	
	*rk_p  = remainder(*rkm2_p, *rkm1_p);

	while (! rk_p->isZero() )
	{
		if (rkm1_p->getExactDegree()<minDeg )
		{
			//std::cerr <<" (rk_p->getExactDegree()- minDeg )" << (rkm1_p->getExactDegree()- minDeg ) << std::endl;
			return false;
		}
		rkm2_p = rkm1_p;
		rkm1_p = rk_p;
		rk_p = rkm2_p;
		*rk_p = remainder( *rkm2_p, *rkm1_p );
	}
	pol1_ref= *rkm1_p;

	if (rkm1_p->getExactDegree()<minDeg )
		{
			//std::cerr <<" (rk_p->getExactDegree()- " << minDeg <<")" << (rkm1_p->getExactDegree()- minDeg ) << std::endl;
			return false;
		}

	return true;
}


template<class TPolynomX, class TRing>
TPolynomX		UnivariatePolynomialRing<TPolynomX, TRing>::diff( const	TPolynomX 			& pol)	const
{
	int eDeg =  pol.getExactDegree(); 
	TPolynomX	res( std::max(0, eDeg-1 ) );
	for (int x_exp= eDeg; x_exp>0; x_exp--)
		res.setCoeff( x_exp-1 ,    ring_ref_m.multiply( ring_ref_m.Convert(x_exp) , pol.getCoeff( x_exp ) ));
	
	return res;
};


  template<class TPolynomX, class TRing>
        inline std::string  UnivariatePolynomialRing<TPolynomX, TRing>::coeffToGeneratorExponentGAPStr(const CoeffType z1) const
        {
            if (z1.isZero() )
                return "z";
         
            std::stringstream oss;
           
             ( ring_ref_m.elemToGeneratorExponent(z1) ) >> oss ;
            
            return oss.str();   
        }



template<class TPolynomX, class TRing>
void		UnivariatePolynomialRing<TPolynomX, TRing>::diffInPlace( TPolynomX 			& pol)	const
{
	int eDeg =  pol.getExactDegree(); 
	 
	for (int x_exp= 0; x_exp<eDeg; x_exp++)
		pol.setCoeff( x_exp ,    ring_ref_m.multiply( ring_ref_m.Convert(x_exp+1) , pol.getCoeffConst( x_exp+1 ) ));

	pol.setCoeff( eDeg,TRing::ElementType::Zero);

	return ;
};


template<class TPolynomX, class TRing>
typename UnivariatePolynomialRing<TPolynomX, TRing>::CoeffType 		UnivariatePolynomialRing<TPolynomX, TRing>::evalAt(const TPolynomX & polynom1 ,  CoeffType  el) const
{
	typename TRing::ElementType res=0;
	int eDeg=polynom1.getDegree();
	for (int x_exp= 0; x_exp<=eDeg; x_exp++)
	{
 		ring_ref_m.addInPlace(res, ring_ref_m.multiply( polynom1.getCoeffConst(x_exp),ring_ref_m.pow(el,x_exp)));
	}
	return res;

}

template<class TPolynomX, class TRing> 
TPolynomX       UnivariatePolynomialRing<TPolynomX, TRing>::subst(  const TPolynomX     & polynom,    const TPolynomX      & substPol      ) const
{

    TPolynomX respol;
        int eDeg=polynom.getDegree();
     for (int x_exp= 0; x_exp<=eDeg; x_exp++)   
    {
        respol = add(respol, scalarMultiply(polynom.getCoeff(x_exp), pow(substPol,x_exp)));
    }
    return respol;
}

