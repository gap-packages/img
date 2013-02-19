// polynom.cpp: Implementierung der Klasse polynom.
//
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// Polynomxy
//////////////////////////////////////////////////////////////////////
template <class TNum, class TIndex>
polynomXY<TNum, TIndex>::polynomXY( ) : 	maxDegree(-1),
							maxDegreePlusOne( 0 ),
							idefs(-1),
							size( idefs.getSize() ),
						koeff(NULL)
{

	name	= std::string("");
   
}



/// create a polynom in (x,y). with maximal monom degree = <b>_maxDegree</b>
template <class TNum, class TIndex>
polynomXY<TNum, TIndex>::polynomXY(const short _maxDegree) : 	maxDegree(_maxDegree),
							maxDegreePlusOne( _maxDegree + 1 ),
							idefs( TIndex(_maxDegree)),
							size(idefs.getSize() )
									
{
	assert(_maxDegree>=0); // erstens das, und zweitens kann es theoretisch zum ueberlauf kommen, auch fuer maxDegreePlusOne

	name	= std::string("");
	koeff 	= new TNum[ idefs.getSize() ] ;
 
	for (short dim=idefs.getSize()-1; dim>=0; dim--)
		koeff[dim] = TNum::Zero;
}

/// create a polynom in (x,y). with maximal monom degree= <b>_maxDegree</b>
template <class TNum, class TIndex>
polynomXY<TNum, TIndex>::polynomXY(string _name, const short _maxDegree) :	maxDegree(_maxDegree),
													maxDegreePlusOne( _maxDegree+1 ),
													idefs( TIndex(_maxDegree) ),
													size(idefs.getSize() ),
													name( _name)
{
	assert(_maxDegree>=0); // erstens das, und zweitens kann es theoretisch zum ueberlauf kommen
	
	koeff 	= new TNum[ idefs.getSize() ] ;

	for (short dim = idefs.getSize()-1; dim >= 0 ; dim--)
		koeff[dim] = TNum::Zero;
}


template <class TNum, class TIndex>
polynomXY<TNum, TIndex>::polynomXY(const polynomXY<TNum, TIndex>& pxy) :	maxDegree(pxy.maxDegree),
											maxDegreePlusOne(pxy.maxDegree+1),
											idefs(pxy.idefs),
											size( pxy.size ),
											name(pxy.name)
{								
	koeff 	= new TNum[idefs.getSize() ] ;

	for (short dim=idefs.getSize()  - 1 ; dim>=0; dim--)
		koeff[dim] = pxy.koeff[dim];
}



template <class TNum, class TIndex>
inline int polynomXY<TNum, TIndex>::getIndex(const short x_exp, const short y_exp) const
{
	return idefs.getPairIndex( x_exp, y_exp );
}



template <class TNum, class TIndex>
inline TNum const * 	polynomXY<TNum, TIndex>::getCoeffConstAddr(const short x_exp, const short y_exp) const
{
	#ifdef SAFE
		testBounds(x_exp, y_exp);
	#endif
	return &koeff[ getIndex(x_exp, y_exp) ];
}

template <class TNum, class TIndex>
inline TNum * 		polynomXY<TNum, TIndex>::getCoeffAddr(const short x_exp, const short y_exp)
{
	#ifdef SAFE
		testBounds(x_exp, y_exp);
	#endif
	return &koeff[ getIndex(x_exp, y_exp) ];
}


/// reset max possible degree of a contained (x,y)-monom. All data is erased !
template <class TNum, class TIndex>
void polynomXY<TNum, TIndex>::setDegree(const short _maxDegree)
{
	assert(_maxDegree>=0); // erstens das, und zweitens kann es theoretisch zum ueberlauf kommen

	#ifdef DEBUG
		std::cerr << "polynomXY<TNum>::setDegree"  << std::endl;
		std::cerr << "setDegree"  << std::endl;
		std::cerr << "_maxDegree" << _maxDegree << std::endl;
	#endif

	delete[] koeff;

	maxDegree	 = _maxDegree;
	maxDegreePlusOne = _maxDegree+1;
	idefs=TIndex(_maxDegree);

	size  = idefs.getSize() ;

	koeff = new TNum[idefs.getSize() ] ;

	for (short dim=idefs.getSize() -1; dim>=0; dim--)
		koeff[dim]=TNum::Zero;

}


template <class TNum, class TIndex>
short polynomXY<TNum, TIndex>::getMaxDegree() const
{
	return maxDegree;
}

template <class TNum, class TIndex>
bool polynomXY<TNum, TIndex>::operator==(const polynomXY<TNum, TIndex>& pxy) const
{

	assert(&pxy!=NULL);
	for (short currDegree=0; currDegree <= maxDegree; currDegree++)
		for(short yexp=0; yexp <= currDegree; yexp++)
			if (  getCoeff( currDegree-yexp, yexp) != pxy.getCoeff( currDegree-yexp, yexp) )
				return false;
	return true;
	 
}

template <class TNum, class TIndex>
polynomXY<TNum, TIndex>& polynomXY<TNum, TIndex>::operator=(const polynomXY<TNum, TIndex>& pxy)
{
	if ( this!=&pxy )
	{

		assert(&pxy!=NULL);

		maxDegree	 = pxy.maxDegree;
		maxDegreePlusOne = pxy.maxDegreePlusOne;

		idefs=pxy.idefs;
		size = pxy.size;
		
		if (koeff!=NULL) delete[] koeff;

		koeff = new TNum[idefs.getSize() ] ;

		for (short dim = idefs.getSize()  - 1; dim >= 0; dim--)
			koeff[dim] = pxy.koeff[dim];
	}
	return *this;
}


template <class TNum, class TIndex>
polynomXY<TNum, TIndex>::~polynomXY()
{
	if (koeff!=NULL)
		delete[] koeff;
	koeff = NULL;
}



template <class TNum, class TIndex>
void polynomXY<TNum, TIndex>::output(std::ostream& os) const
{
	int i, j;
	TNum 	z;

	for (i=0; i<=maxDegree; i++)
		for(j=0; j<=i; j++)
		{
			z = getCoeff( i-j, j); 
			if (! z.isZero())
				os << "(" << z <<")" << "x^" << i-j << "y^" << j << " + ";
		}
}



/// output polynom in Macaulay-Style 
template <class TNum, class TIndex>
void polynomXY<TNum, TIndex>::print(std::ostream& os) const
{
	int i, j;

	TNum 	z;

	bool first = true;

	for (i=0; i<=maxDegree; i++)
		for( j=0; j<=i; j++)
		{
			z = getCoeff(i-j, j); 
			if ( z.isNotZero() )
			{
				
				if (first)
					first=false;
				else
				{
					os << " + ";
				}
				 
				os << "(" << z << ")*x^" << i-j << "*y^" << j ;
			}
		}
	//os << ";\n  ";
}


template <class TNum, class TIndex>
void polynomXY<TNum, TIndex>::printInMacaulayStyle(std::ostream& os) const
{
	 os << name << "=" ;
	print(os);
	 os <<  ";" ;
}


/// print object ( debug)
template <class TNum,class TIndex>
void polynomXY<TNum, TIndex>::outputMatrix() const
{
	int i, j;

	TNum 	z;

	std::cerr << std::endl <<"Name: " << name  << " degree: " << maxDegree << std::endl << "-----------------" << std::endl;
	for (i=0; i<=maxDegree; i++)
	{
		std::cerr << "i" << i << " ";
		for(j=0; j<=i; j++)
		{
			z =getCoeff(i-j, j); 
			if ( z.isNotZero() )
				std::cerr << z << "*x^" << i-j << "*y^" << j << " + " ; 
		}
		std::cerr << std::endl;
	}
}

/** @brief output coefficients in monomgroups with equal degree  
*/
template <class TNum, class TIndex>
void polynomXY<TNum, TIndex>::OutputPureCoefficients(ostream &OStream, int monomDegree, bool mitKomma) const
{
	int j;

	TNum	 z;

	OStream << "{";
	for(j=0; j<=monomDegree; j++)
		{
			z = getCoeff(monomDegree-j, j); 
			OStream << z << " ";
			if (j<monomDegree)
				OStream << ", ";
		}

	OStream << "} ";
	if ( mitKomma==true )
		 (OStream) << ", ";

	OStream << "	-- ";

	for(j=0; j<=monomDegree; j++)
	{
		z =getCoeff(monomDegree-j, j); 

		OStream << "(" << z <<")*x^" << monomDegree-j << "y^" << j;
		if (j<monomDegree)
			OStream << " + ";
	}


	OStream << std::endl;
}


/// Test, if x_exp and y_exp have legal values due do maxDegree !
template <class TNum, class TIndex>
inline void polynomXY<TNum, TIndex>::testBounds(const short x_exp,const short y_exp) const
{
	if ( !( getIndex(x_exp, y_exp)>=0) )
	{
		std::cerr <<"polynomXY::setCoeff() : Err!" << std::endl;
	}
	if ( !( (getIndex(x_exp, y_exp)) < ( maxDegreePlusOne*maxDegreePlusOne )) )
	{
		std::cerr <<"polynomXY::setCoeff() : Err!" << std::endl;
		std::cerr << "(x_exp*(maxDegree+1)+y_exp)=" << (x_exp*(maxDegreePlusOne)+y_exp) << std::endl;
		std::cerr << "(maxDegree+1)*(maxDegree+1)=" << (maxDegreePlusOne)*(maxDegreePlusOne);
	}
	assert( (x_exp + y_exp) <=maxDegree);
	assert( x_exp>=0 && y_exp>=0);
}


template <class TNum, class TIndex>
inline void polynomXY<TNum, TIndex>::setCoeff(const short x_exp,const short y_exp, const TNum value)
{
	#ifdef SAFE
		testBounds(x_exp,y_exp);
	#endif

	koeff[ getIndex(x_exp, y_exp) ] = value;
};


template <class TNum, class TIndex>
inline TNum polynomXY<TNum, TIndex>::getCoeff(const short x_exp, const short y_exp) const
{
	#ifdef SAFE
		testBounds(x_exp,y_exp);
	#endif
	return koeff[ getIndex(x_exp, y_exp) ];
}


template <class TNum, class TIndex>
inline TNum const polynomXY<TNum, TIndex>::getCoeffConst(const short x_exp,const short y_exp) const
{
	#ifdef SAFE
		testBounds(x_exp,y_exp);
	#endif
	return koeff[ getIndex(x_exp, y_exp) ];
}

template <class TNum, class TIndex>
inline const TNum & polynomXY<TNum, TIndex>::getCoeffConstRef(const int x_exp,const int y_exp) const
{
	#ifdef SAFE
		testBounds(x_exp,y_exp);
	#endif
	return koeff[ getIndex(x_exp, y_exp) ];
}


template <class TNum, class TIndex>
inline TNum& polynomXY<TNum, TIndex>::getCoeffRef(const short x_exp, const short y_exp)
{
	#ifdef SAFE
		testBounds(x_exp,y_exp);
	#endif
	return koeff[ getIndex(x_exp, y_exp) ];
}


template <class TNum, class TIndex>
inline  short polynomXY<TNum, TIndex>::getDegree() const
{
	return maxDegree;
};

template <class TNum, class TIndex>
inline void polynomXY<TNum, TIndex>::clear(short _grad)
{
	#ifdef SAFE
	assert (_grad<=maxDegree);
	#endif 
	for (int i=0; i<= _grad; i++)
		for(int j=0; j<=i; j++)
			koeff[ getIndex( i-j, j ) ] = TNum::Zero;


}

template <class TNum, class TIndex>
inline void polynomXY<TNum, TIndex>::clear()
{
	for (int dim=0; dim<idefs.getSize() ; dim++)
	{
		koeff[dim] = TNum::Zero;
	} 
}


  

///////////////////////////////////////////////////////////////////
// polynomx
///////////////////////////////////////////////////////////////////
 
template <class TNum>
const polynomx<TNum>   polynomx<TNum>::Zero ( polynomx<TNum>::getZero() );

template <class TNum>
const polynomx<TNum>   polynomx<TNum>::One ( polynomx<TNum>::getOne() );


template <class TNum>
polynomx<TNum >::polynomx():    maxDegree(-1), size(0)
{
	 
	koeff = NULL ;
}


template <class TNum>
polynomx<TNum >::polynomx(const int gr):    maxDegree(gr), size(gr+1)
{
	assert(gr>=0); // erstens das, und zweitens kann es theoretisch zum ueberlauf kommen

	koeff = new TNum[size] ;
	for (int dim=size-1; dim>=0; dim--)
		koeff[dim] = TNum::Zero;
}



template <class TNum>
polynomx< TNum >::polynomx(const polynomx<TNum> & fpx) : maxDegree(fpx.maxDegree),  size(fpx.maxDegree + 1 )
{
	koeff = new TNum[size] ;
	for (int dim=size-1; dim>=0; dim--)
		koeff[dim] = fpx.koeff[dim];
}


template <class TNum>
polynomx< TNum >::polynomx(const vector<TNum> & vecx) : maxDegree( vecx.size()-1 ),  size( vecx.size()  )
{
	koeff = new TNum[size] ;
	for (int dim=size-1; dim>=0; dim--)
		koeff[dim] = vecx[dim];
}


		
template <class TNum>
inline polynomx<TNum> & polynomx< TNum >::operator=(const polynomx<TNum> & fpx)
{
	if ( this != &fpx )
	{
		assert(&fpx!=NULL);
		//assert(maxDegree	 == fpx.maxDegree)
		maxDegree	 = fpx.maxDegree;
		size = fpx.size;
		
		if (koeff!=NULL) 
			delete[] koeff;

		koeff = new TNum[size] ;

		for (short dim=size-1; dim>=0; dim--)
			koeff[dim] = fpx.koeff[dim];
	}
	return *this;

}


template <class TNum>
inline void 	polynomx<TNum >::testbounds(const int x_exp) const
{
	if (x_exp<0 || x_exp>=maxDegree+1)
	{
			std::cerr <<"polynomx::getCoeff(): Err";
			std::cerr << "x_exp =" << x_exp;
			std::cerr << "maxDegree+1= " << maxDegree+1;
			std::cerr << fflush;
	}
	assert( x_exp <= maxDegree );
	assert( x_exp >= 0 );
}

 

template <class TNum>
inline TNum 	polynomx<TNum >::getCoeff(const int x_exp)  const
{
	#ifdef SAFE
		testbounds(x_exp);
	#endif

	return koeff[x_exp];
}



template <class TNum>
inline const TNum& 	polynomx<TNum >::getCoeffConstRef(const int x_exp)  const
{	
	#ifdef SAFE
		testbounds(x_exp);
	#endif

	return koeff[x_exp];
}


template <class TNum>
inline  TNum& 	polynomx<TNum >::getCoeffRef(const int x_exp)  
{	
	#ifdef SAFE
		testbounds(x_exp);
	#endif

	return koeff[x_exp];
}


template <class TNum>
inline  TNum const polynomx<TNum >::getCoeffConst(const int x_exp)  const
{	
	#ifdef SAFE
		testbounds(x_exp);
	#endif



	return koeff[x_exp];
}

template <class TNum>
inline  TNum polynomx<TNum >::getSafeCoeff(const int x_exp)  const
{	
	#ifdef SAFE
		testbounds(x_exp);
	#endif
	if (x_exp>maxDegree)
		return TNum::Zero;

	return koeff[x_exp];
}

template <class TNum>
inline  TNum const polynomx<TNum >::getSafeCoeffConst(const int x_exp)  const
{	
	#ifdef SAFE
		//testbounds(x_exp);
	#endif
	if (x_exp>maxDegree)
		return TNum::Zero;


	return koeff[x_exp];
}


template <class TNum>
inline void 	polynomx<TNum >::setCoeff(const int x_exp, const TNum& value)
{
	#ifdef SAFE
		testbounds(x_exp);
	#endif

	koeff[x_exp] = value;

}


template <class TNum>
inline const int 	polynomx<TNum >::getDegree() const 
{
	return(maxDegree);

};

template <class TNum>
inline const int 	polynomx<TNum >::getExactDegree() const 
{
	//return(maxDegree);

	for (int deg= maxDegree; deg>=0; deg--)
		if (koeff[deg] != TNum::Zero)
			return deg;

	return 0;
};



template <class TNum>
inline const bool 	polynomx<TNum >::isConstant() const 
{
	//return(maxDegree);

	for (int deg= maxDegree; deg>0; deg--)
		if (koeff[deg] != TNum::Zero)
			return false;

	return true;
};


template <class TNum>
inline const bool 	polynomx<TNum >::isZero() const 
{
	for (int x_exp= maxDegree; x_exp>=0; x_exp--)
		if (koeff[x_exp] != TNum::Zero)
			return false;

	return true;

};


template <class TNum>
bool polynomx<TNum>::operator==(const polynomx<TNum>& px) const
{
	assert(&px!=NULL);
	for (short currDegree=0; currDegree <= maxDegree; currDegree++)
		if (  getCoeff( currDegree ) != px.getCoeff( currDegree) )
				return false;
	return true;
	 
}

template <class TNum>
bool polynomx<TNum>::operator!=(const polynomx<TNum>& px) const
{
	assert(&px!=NULL);
	for (short currDegree=0; currDegree <= maxDegree; currDegree++)
		if (  getCoeff( currDegree ) != px.getCoeff( currDegree) )
				return true;
	return false;
	 
}

template <class TNum>
inline const bool 	polynomx<TNum >::isOne() const 
{
	for (int x_exp= maxDegree; x_exp>0; x_exp--)
		if (koeff[x_exp] != TNum::Zero)
			return false;

	return  (koeff[0] == TNum::One);

};




template <class TNum>
polynomx<TNum >::~polynomx()
{
	delete[] koeff;
}


template <class TNum>
void polynomx<TNum >::clear(int _grad)
{
	for (int dim=0; ( dim<size && dim<=_grad) ; dim++)
		koeff[dim] = TNum::Zero;
}



template <class TNum>
void polynomx<TNum >::clear()
{
	for (int dim=0; dim<size ; dim++)
		koeff[dim] = TNum::Zero;
}

template <class TNum>
std::string polynomx<TNum>::getStringRep() const
{
    int i;

    TNum    z;

    bool first = true;

    std::stringstream strstream;

    for (i=0; i<=maxDegree; i++)
    {
         
        z = getCoeff(i); 
        if ( z.isNotZero() )
        {
            if (first)
                first=false;
            else
            {
                strstream << " + ";
            }
                
            strstream << "(" << z << ")*x^" << i ;
        }
    }
    return strstream.str();
    //os << ";\n  ";
}
/// output polynom in Macaulay-Style 
template <class TNum>
void polynomx<TNum>::print(std::ostream& os) const
{
	int i;

	TNum 	z;

	bool first = true;

	for (i=0; i<=maxDegree; i++)
	{
		 
		z = getCoeff(i); 
		if ( z.isNotZero() )
		{
			if (first)
				first=false;
			else
			{
				os << " + ";
			}
				
			os << "(" << z << ")*x^" << i ;
		}
	}
		
	//os << ";\n  ";
}
