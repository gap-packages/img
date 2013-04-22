/** @file fast_Ring.cpp
*
* @brief contains fast_Ring implementation
*/

#include <cstdlib>
#include <exception>
#include <new>


/// @todo hier pruefen, ob tNum gross genug ausgelegt ist
template <class TNum, class kdefs>
fast_Ring<TNum,kdefs>::fast_Ring(	unsigned short _char, 
						unsigned short _epsPrec) : 	characteristic(_char),
											epsilon(_epsPrec), 
											generator(getGenerator())
{
	assert(TNum::wellDefined(characteristic));
	assert( epsilon <= 1 );
    assert( isGenerator(generator));
	init();
	assert( wellDefined() );
}

template <class TNum, class kdefs>
fast_Ring<TNum,kdefs>::fast_Ring(   unsigned short _char, 
                        unsigned short _epsPrec, short _generator) :  characteristic(_char),
                                            epsilon(_epsPrec), 
                                            generator( Convert(_generator) )
{
    assert( TNum::wellDefined(characteristic) );
    assert( epsilon <= 1 );
    assert( isGenerator(generator));
    init();
    assert( wellDefined() );
}

template <class TNum, class kdefs>
bool 		fast_Ring<TNum,kdefs>::wellDefined()
{
	TNum num;
	assert (num.wellDefined(characteristic) );
	return true;
}

/// @todo warum  dieser Init-Kram ? gibt sowieso nur einen Konstruktor, aber nun gut...
template <class TNum,class kdefs>
inline void	fast_Ring< TNum, kdefs>::init() 
{
	additiveInverseTable		=NULL;
	multiplicationTable		=NULL;
	additionTable			=NULL;
	multiplicativeInverseTable	=NULL;
	
	fastAdditionTable		=NULL;
	elementsToExponentsTab		=NULL;
	exponentsToElementTab		=NULL;

	sqrtTable		=NULL;
	
	moduloTable = NULL;

	bContainsImagNum_m=false;

	imagNum_m=0;

	try
	{
		if (additionTable==NULL)
			additionTable = createAdditionTable();
	
		if (additiveInverseTable==NULL)
			additiveInverseTable = createAdditiveInverseTable();
	
		if (multiplicationTable==NULL)
			multiplicationTable = createMultiplicationTable();
		
	
		if (multiplicativeInverseTable==NULL)
			multiplicativeInverseTable = createMultiplicativeInverseTable();
	
		if (fastAdditionTable==NULL)
			fastAdditionTable = createFastAdditionTable();
	
		if (elementsToExponentsTab==NULL)
			elementsToExponentsTab = initElementsToExponentsTab(generator);
	
		if (exponentsToElementTab==NULL)
			exponentsToElementTab = initExponentsToElementTab(generator);

		if (sqrtTable==NULL)
			sqrtTable = createSqrtTable( );

		if (moduloTable==NULL)
			moduloTable=createModuloTable();
	}
	catch(std::bad_alloc &e)
	{
		std::cerr << "Memory allocation error in fast_Ring::init() " << std::endl;
		exit(0);
	}
	catch(...)
	{
		std::cerr << "Error in creating ring operation Tables, " << std::endl;
		exit(0);
	}
}



template <class TNum,class kdefs>
fast_Ring<TNum,kdefs>::~fast_Ring()
{
	if (additiveInverseTable!=NULL) 		{ 	delete[] additiveInverseTable;	
									additiveInverseTable = NULL;	}

	if (multiplicationTable!=NULL) 		{	 delete[] multiplicationTable;	
									multiplicationTable = NULL;	}

	if (additionTable!=NULL) 			{ 	delete[] additionTable;	 
									additionTable = NULL;		}

	if (multiplicativeInverseTable!=NULL) 	{	 delete[] multiplicativeInverseTable;
									 multiplicativeInverseTable = NULL;	}
	
	if ( fastAdditionTable!=NULL   ) 		{ 	delete[] 	fastAdditionTable;	 
									fastAdditionTable = NULL;		}

	if ( elementsToExponentsTab!=NULL ) 	{ 	delete[] 	elementsToExponentsTab; 
									elementsToExponentsTab = NULL;	}

	if ( exponentsToElementTab!=NULL ) 		{ 	delete[] 	exponentsToElementTab;	 
									exponentsToElementTab = NULL;	}

	if ( sqrtTable!=NULL ) 			{ 	delete[] 	sqrtTable;	 
									sqrtTable = NULL;	}

	if ( moduloTable!=NULL ) 			{ 	delete[] 	moduloTable;	 
									moduloTable = NULL;	}

}





//----------------------Convert---------------------------------------------------


/**

@param src TConvNum object src  must implement getX() and getEps() -Interface .
 In general this can be replaced by getValue(epsPrecision)-Interface, 
but is not neccessary in this specialized class
*/

template <class TNum, class kdefs>
template <class TConvNum >
inline  TNum fast_Ring<TNum, kdefs>::Convert(const TConvNum  src ) const
{
	#ifdef SAFE
		assert( getEpsPrecision()<=1 );
	#endif
	
	TNum 	z;
	
	z.setX( ConvertScalar( src.getX() ) );
	//alternativ: 
	//int 	currEpsPrecision = 0;
	//z.setValue(currEpsPrecision, ConvertScalar( a.getValue(currEpsPrecision) ) );
	
	if ( getEpsPrecision()==1 )
	{
		z.setEps( ConvertScalar(src.getEps() ) );
		//alternativ: 
		//currEpsPrecision=1;
		//z.setValue(currEpsPrecision, ConvertScalar( a.getValue(currEpsPrecision) ) );
	}
	return z;
}

template <class TNum, class kdefs>
inline  TNum fast_Ring<TNum, kdefs>::Convert(const double  a) const
{
	#ifdef SAFE
		assert( getEpsPrecision()<=1 );
	#endif
	
	TNum 	z;
	
	z.setX( ConvertScalar( a ) );
	
	return z;
}

template <class TNum, class kdefs>
inline int  fast_Ring<TNum, kdefs>::repToInt(const  TNum & z )    const
{
    assert( getEpsPrecision()==0 );
    return z.getX();
}


template <class TNum, class kdefs>
inline  TNum fast_Ring<TNum, kdefs>::Convert(const int  a) const
{
	#ifdef SAFE
		assert( getEpsPrecision()<=1 );
	#endif
	
	TNum 	z;
	
	z.setX( ConvertScalar( a ) );
	
	return z;
}

template <class TNum, class kdefs>
inline  TNum fast_Ring<TNum, kdefs>::Convert(const short  a) const
{
    #ifdef SAFE
        assert( getEpsPrecision()<=1 );
    #endif
    
    TNum    z;
    
    z.setX( ConvertScalar( a ) );
    
    return z;
}

template <class TNum, class kdefs>
inline  TNum fast_Ring<TNum, kdefs>::Convert(const unsigned long  a) const
{
    #ifdef SAFE
        assert( getEpsPrecision()<=1 );
    #endif
    
    TNum    z;
    
    z.setX( ConvertScalar( a ) );
    
    return z;
}

/*
template <class TNum, class kdefs>
template < const int>
inline  TNum fast_Ring<TNum, kdefs>::Convert(const int  a) const
{
	#ifdef SAFE
		assert( getEpsPrecision()<=1 );
	#endif
	
	TNum 	z;
	
	z.setX( ConvertScalar( a ) );
	//alternativ: 
	//int 	currEpsPrecision = 0;
	//z.setValue(currEpsPrecision, ConvertScalar( a.getValue(currEpsPrecision) ) );
	
	return z;
}*/

/// riscy in some cases? does currently ignore the case where TConvNum.getEpsPrecision > 1.
/// @TODO 	assert( getEpsPrecision()<=1 ); gehoert in den Konstruktor dieser Klasse, 
///		da sie offenbar nicht mit EPSPrecision >1 arbeiten kann und muss. 
///		An anderen Stellen kann diese Abfrage entfernt werden. 
template <class TNum,class kdefs>
template <class TConvNum >
inline  void fast_Ring<TNum, kdefs>::ConvertInPlace( TConvNum & a) const
{

	assert( a.getEpsPrecision()<=1 );
	
	
	a.setX( ConvertScalar(a.getX() ) );
	
	if ( getEpsPrecision()==1 )
		a.setEps( ConvertScalar(a.getEps() ) );

	#ifdef SAFE
		assert(a.getX()   ==  ConvertScalar(a.getX() ) );
		assert(a.getEps() ==  ConvertScalar(a.getEps() ) );
	#endif

	return ;
}

/// @todo TNum::scalarType zurueckgeben oder bei int bleiben?
template <class TNum, class kdefs>
inline int	 fast_Ring<TNum,kdefs>::ConvertScalar(const int  a) const
{
	int res = a;
	while (res<0)
	{
		res += getCharacteristic();
	}
	if ( res >= getCharacteristic() ) 
	{	
		res %= getCharacteristic();
	}
	// Alternativ zur Modulo-Rechnung:
	//while (res >=getCharacteristic()<)
	//{
	//	res -= getCharacteristic();
	//}
	return res;
 }


/// convert a integer to Ring element. Assumption: <b>a </b> >=0
template <class TNum,class kdefs>
inline int fast_Ring<TNum,kdefs>::ConvertScalarSpec(const int a) const
{
	#ifdef SAFE
		assert(a>=0);
	#endif
	return 	a % getCharacteristic();
 }


 
/** @brief convert a integer to Ring Elemen. Assumption: <b>a </b> >=-getCharacteristic()
/// @pre assumption: a >=-getCharacteristic()
*/
 template <class TNum,class kdefs>
inline  int fast_Ring<TNum,kdefs>::FastConvertScalar(const int  a) const
{
	return ( a + getCharacteristic() ) % getCharacteristic();
 }








//----------------------Add---------------------------------------------------

template <class TNum,class kdefs>
inline TNum 	fast_Ring< TNum, kdefs>::add( const TNum a, const  TNum b) const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		assert(Convert(b)==b);
	#endif
	return additionTable[ getPairIndex(a, b) ];
}




template <class TNum,class kdefs>
inline TNum 	fast_Ring< TNum, kdefs>::addRef( const TNum &a, const  TNum &b) const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		assert(Convert(b)==b);
	#endif
	return additionTable[ getPairIndexByRef(a, b) ];
}





//----------------------Add in place---------------------------------------------------

template <class TNum,class kdefs>
inline void fast_Ring< TNum, kdefs>::addInPlace(TNum& a, const TNum b) const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		assert(Convert(b)==b);
	#endif

	a= this->additionTable[ getPairIndex(a, b) ];
	return;

	
}

	
template <class TNum,class kdefs>
inline void fast_Ring< TNum, kdefs>::addInPlaceRef(TNum& a, const TNum &b) const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		assert(Convert(b)==b);
	#endif

	a = this->additionTable[ getPairIndexByRef(a, b) ];
	return;
}




//----------------------Additive Inverse---------------------------------------------------


template <class TNum,class kdefs>
inline TNum  fast_Ring< TNum, kdefs>::addInv(const TNum a) const
{
	#ifdef SAFE
		assert(Convert(a)==a);
	#endif
	return additiveInverseTable[ getSingleIndex(a) ];
}




template <class TNum,class kdefs>
inline TNum  fast_Ring< TNum, kdefs>::addInvRef(const TNum& a) const
{
	#ifdef SAFE
		assert(Convert(a)==a);
	#endif
	return additiveInverseTable[ getSingleIndexByRef(a) ];
}



template <class TNum, class kdefs>
inline void	fast_Ring< TNum, kdefs>::addInvInPlace( TNum  & a) const
{
	#ifdef SAFE
		assert(Convert(a)==a);
	#endif
	a=additiveInverseTable[ getSingleIndexByRef(a) ];
	return;
}






//-------------------------------multiply-------------------------------------------





template <class TNum,class kdefs>
inline  TNum   fast_Ring< TNum, kdefs>::multiply(const TNum a, const TNum b) const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		assert(Convert(b)==b);
	#endif

	return multiplicationTable[getPairIndex(a, b)];
}



// inline bringt hier kaum Zeitvorteil.
///multipliziert zwei Zahlen
template <class TNum,class kdefs>
 inline  TNum  fast_Ring< TNum, kdefs>::multiplyRef(const TNum & a, const TNum & b) const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		assert(Convert(b)==b);
	#endif

	return multiplicationTable[getPairIndexByRef(a, b)];
}



//-------------------------------multiply in place-------------------------------------------



template <class TNum,class kdefs>
inline void  fast_Ring< TNum, kdefs>::multiplyInPlace(TNum &a, const TNum b) const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		assert(Convert(b)==b);
	#endif
	a = multiplicationTable[ getPairIndex(a, b) ];
}


template <class TNum,class kdefs>
inline void  fast_Ring< TNum, kdefs>::multiplyInPlaceRef(TNum &a, const TNum &  b) const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		assert(Convert(b)==b);
	#endif
	a = multiplicationTable[ getPairIndexByRef(a, b) ];
}

//-------------------------------scalarmultiply-------------------------------------------





template <class TNum,class kdefs>
inline  TNum   
fast_Ring< TNum, kdefs>::scalarMultiply(const FieldType::ElementType a, const TNum b) const
{
	#ifdef SAFE
		assert(Convert(TNum(a))==a);
		assert(Convert(b)==b);
	#endif

	return multiplicationTable[getPairIndex(a, b)];
}



// inline bringt hier kaum Zeitvorteil.
///multipliziert zwei Zahlen
template <class TNum,class kdefs>
 inline  TNum  
fast_Ring< TNum, kdefs>::scalarMultiplyRef(const FieldType::ElementType & a, const TNum & b) const
{
	#ifdef SAFE
		assert(Convert(TNum(a))==a);
		assert(Convert(b)==b);
	#endif

	return multiplicationTable[getPairIndexByRef(a, b)];
}



//-------------------------------scalarmultiply in place-------------------------------------------


template <class TNum,class kdefs>
inline void  
fast_Ring< TNum, kdefs>::scalarMultiplyInPlace(const FieldType::ElementType 	a,
												TNum & b ) const
{
	#ifdef SAFE
		assert(Convert(TNum(a))==a);
		assert(Convert(b)==b);
	#endif
	b = multiplicationTable[ getPairIndex(a, b) ];
}


template <class TNum,class kdefs>
inline void  
fast_Ring< TNum, kdefs>::scalarMultiplyInPlaceRef(const FieldType::ElementType & a,
												 TNum &  b 	) const
{
	#ifdef SAFE
		assert(Convert(TNum(a))==a);
		assert(Convert(b)==b);
	#endif
	b = multiplicationTable[ getPairIndexByRef(a, b) ];
}

//-------------------------------multiply by exponents-------------------------------------------


template <class TNum,class kdefs>
inline TNum const   fast_Ring< TNum, kdefs>::multByExp(const TNum a, const TNum b) const
{
	if (a.isZero() || b.isZero() )
		return TNum::Zero;

	register short res = a.getX()+b.getX();
	if ( res>=getCharacteristic() )
		return res-getCharacteristic();
	return res;
}



template <class TNum,class kdefs>
inline TNum const   fast_Ring< TNum, kdefs>::multByExpRef(const TNum & a, const TNum & b) const
{
	if (a.isZero() || b.isZero() )
		return TNum::Zero;

	register short res = a.getX()+b.getX();
	if (res>=getCharacteristic())
		return res-getCharacteristic();
	return res;
}



template <class TNum,class kdefs>
inline void  fast_Ring< TNum, kdefs>::multByExpInPlace( TNum & a, const TNum  b) const
{
	if (a.isZero() || b.isZero() )
		 a=TNum::Zero;

	register short res = a.getX()+b.getX();
	if (res>=getCharacteristic())
		a=res-getCharacteristic();
	a=res;
}



template <class TNum,class kdefs>
inline void  fast_Ring< TNum, kdefs>::multByExpInPlaceRef( TNum & a, const TNum & b) const
{
	if (a.isZero() ||b.isZero())
		 a=TNum::Zero;

	register unsigned short res = a.getX() + b.getX();
	if (res>getCharacteristic()-1)
		a= res-getCharacteristic();
	return;
}
//-------------------------------multiplicative inverse-------------------------------------------

// invers: sollte nur skalare invertieren und nicht zahlen. -
//  no das ist ja ein superring und kann auch mit basicNumber rechnen :-)
template <class TNum,class kdefs>
inline  TNum  fast_Ring< TNum, kdefs>::multInv(const TNum  a) const
{
	#ifdef SAFE
		assert( a==Convert(a) );
	#endif
	TNum const &  res = multiplicativeInverseTable[ getSingleIndex( a)];
	if (res.isNotZero() ) 
	{
		return res;
	}
	else 
	{
		std::cerr << "Multiplicative inverse does not exist!" << std::endl;
		throw "Multiplicative inverse does not exist!" ;
	}
}

// invers: sollte nur skalare invertieren und nicht zahlen.
template <class TNum,class kdefs>
inline  TNum  fast_Ring< TNum, kdefs>::multInvRef(const TNum & a) const
{
	#ifdef SAFE
		assert( a==Convert(a) );
	#endif
	TNum const &	res = multiplicativeInverseTable[ getSingleIndexByRef( a) ];
	if (res.isNotZero() ) 
	{
		return res;
	}
	else 
	{
		std::cerr << "Multiplicative inverse does not exist!" << std::endl;
		throw "Multiplicative inverse does not exist!" ;
	}
}







// invers: sollte nur skalare invertieren und nicht zahlen.
template <class TNum,class kdefs>
inline void fast_Ring< TNum, kdefs>::multInvInPlace( TNum & a) const
{
	#ifdef SAFE
		assert( a==Convert(a) );
	#endif
	TNum const &	res = multiplicativeInverseTable[getSingleIndex( a)];
	if (res.isNotZero() ) 
	{
		a = res;
	}
	else 
	{
		std::cerr << "Multiplicative inverse does not exist!" << std::endl;
		throw "Multiplicative inverse does not exist!" ;
	}
}




//-------------------------------accMmult-------------------------------------------


template <class TNum,class kdefs>
inline void fast_Ring<TNum,kdefs>::accMult( TNum& a ,const TNum b , const TNum c)  const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		assert(Convert(b)==b);
		assert(Convert(c)==c);
	#endif
	#ifdef COUNT
		accMultCount=accMultCount+1;
	#endif
	a=additionTable[ getPairIndex (a, multiplicationTable[ getPairIndex(b, c)] ) ];
}


template <class TNum,class kdefs>
inline void fast_Ring<TNum,kdefs>::accMult( TNum* a ,const TNum b , const TNum c)  const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		assert(Convert(b)==b);
		assert(Convert(c)==c);
	#endif
	#ifdef COUNT
		accMultCount=accMultCount+1;
	#endif

	*a = additionTable[ getPairIndex( *a,multiplicationTable[ getPairIndex(b, c)] ) ];
}



template <class TNum,class kdefs>
inline void fast_Ring<TNum,kdefs>::accMultRef( TNum& a ,const TNum& b , const TNum& c)  const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		assert(Convert(b)==b);
		assert(Convert(c)==c);
	#endif
	#ifdef COUNT
		accMultCount=accMultCount+1;
	#endif

	a=additionTable[ getPairIndexByRef(a, multiplicationTable[ getPairIndex(b, c) ]) ];
}



/// Optimierung : Eigenen Überlegungen nach sollte die Multiplikation an dieser Stelle per Tabelle
/// und die Addition auf der CPU erledigt werden, weil ein Multiplikator fest ist und nur 
/// der andere variiert  ( characteristic*characteristic mögliche Werte)
/// Dagegen variieren bei der Addition beide Zahlen, was zu erheblichen Cache misses führen müsste.
/// Allerdings ist die Laufzeit auf einem Intel-Rechner katastrophal, wenn man die Addition auf
/// der CPU durchführt. Es gab mal eine Konfiguration, wo die Addition auf der CPU schneller war,
///  jetzt ist die Laufzeit in etwa gleich,  oder schlechter mit der Addition auf der CPU
template <class TNum,class kdefs>
inline void 
fast_Ring<TNum,kdefs>::accMultSpec( TNum* const a ,const TNum b , const TNum * const c)  const
{
	#ifdef SAFE
		assert(Convert(*a)==*a);
		assert(Convert(b)==b);
		assert(Convert(*c)==*c);
	#endif

	#ifdef COUNT
		accMultCount=accMultCount+1;
	#endif 

	*a=additionTable[ getPairIndex(*a, multiplicationTable[ getPairIndex(b, *c)]) ];

	/*
	#if EPSPRECISION==1
	
		register TNum tmp=multiplicationTable[getPairIndex(b, *c)];
		register short a2=tmp.getX()+a->getX();
		if (a2>=getCharacteristic())
			a->setX(a2-getCharacteristic());
		else
			a->setX(a2);
		a2=tmp.getEps()+a->getEps();
		if (a2>=getCharacteristic())
			a->setEps(a2-getCharacteristic());
		else
			a->setEps(a2);
	
	#else
		//  on Pentiums following code has catastrophic performance: but on hoech its fast
		int tmp=(*a).getX()+multiplicationTable[getPairIndex(b, *c)].getX();
		if (tmp>=getCharacteristic())
			(a)->setX(tmp-getCharacteristic());
		else
			(a)->setX(tmp);
	#endif*/

}


template <class TNum,class kdefs>
inline 		typename fast_Ring<TNum,kdefs>::sqrtInf_t 	fast_Ring<TNum,kdefs>::sqrt	( const TNum a)  const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		
	#endif
	assert( a.getEps()==0 );
	return 	sqrtTable[ a.getX() ];
}

template <class TNum,class kdefs>
inline 		typename fast_Ring<TNum,kdefs>::sqrtInf_t		fast_Ring<TNum,kdefs>::sqrtRef	( const TNum &a) const
{
	#ifdef SAFE
		assert(Convert(*a)==*a);
	#endif
	assert( a.getEps()==0 );
	return sqrtTable[ a.getX() ];
}


///accMult assumes, that b is not zero, and 
template <class TNum,class kdefs>
inline void fast_Ring<TNum,kdefs>::accMultAddr( TNum* a ,const TNum* b , const TNum* c)  const
{
	#ifdef SAFE
		assert(Convert(a)==a);
		assert(Convert(b)==b);
		assert(Convert(c)==c);
	#endif
	*a = additionTable[getPairIndex(*a, multiplicationTable[ getPairIndex(*b, *c)]) ];
}



/// @todo herausfinden, welchen Datentyp wir nach der %-Operation haben 
///       und ob es Probleme bei unsigned Datentypen gibt.
/// @todo vernuenftige Loesung fuer epsPrecision>0.
/// @todo testen , ob angegebene Charakteristik eine Primzahl ist -
///		 Merkt man das nicht spaetestens wenn es keinen Erzeuger gibt?
template <class TNum,class kdefs>
 TNum fast_Ring<TNum,kdefs>::getGenerator()
{
	//std::cerr << " getGenerator () " << std::endl;

	if (getEpsPrecision()>0)
	{
		std::cerr << " Warning: getGenerator()  is not implemented for epsPrecision>0 !";
		std::cerr << std::endl;

		return TNum::Zero;
	}
	for (int m=1; m<getCharacteristic(); m++)
	{
		TNum erzeuger = m;
		TNum tmp = erzeuger;
		assert (tmp==erzeuger);

		bool erz=true;

		for (int n=0; n<getCharacteristic()-1; n++)
		{

			tmp.setX(     (unsigned int ) 
				        ( (unsigned int )tmp.getX() * (unsigned int )erzeuger.getX() )
					% getCharacteristic()
				);
	
			if ( (tmp.getX()==erzeuger.getX()) && (n<( getCharacteristic() - 2 )) )
			{
				erz=false;
			}
		}

		if (erz)
		{
			tmp = erzeuger;
			return erzeuger;
	
		}
	}
	std::cerr << " Error: kein Erzeuger gefunden !  - da gibt es einen Fehler! " << std::endl;
	exit(0);
	return TNum::Zero;
}

template <class TNum,class kdefs>
 bool fast_Ring<TNum,kdefs>::isGenerator(const TNum & _generator) const
{
        TNum tmp = _generator;
        assert (tmp==_generator);

        bool bIsGenerator=true;
        for (int n=0; n<getCharacteristic()-1; n++)
        {

            tmp.setX(     (unsigned int ) 
                        ( (unsigned int )tmp.getX() * (unsigned int )_generator.getX() )
                    % getCharacteristic()
                );
    
            if ( (tmp.getX()==_generator.getX()) && (n<( getCharacteristic() - 2 )) )
            {
                bIsGenerator = false;
            }
        }
        return bIsGenerator;
}


/** @brief create addition table
*
*  @description Initialisiert die Additionstabelle .
* Es gibt zwei Faelle: epsPrecision==0 und epsPrecision==1
* Bei epsPrecision =0 sind die Eintraege eindimensional (x) und 
* fuer epsPrecision==1 zweidimensional (x,eps).
* 
* Optimierungsmoeglichkeit:
* steht wenig cache zur Verfuegung, so kann der Cache-Bedarf um die Haelfte reduziert werden,
* wenn die Additionstabelle nur Eindimensional angelegt wird und die add-Funktionen 
* bei epsPrecision==1  (x1,eps1)+(x2,eps2) simulieren. 

*/
template <class TNum,class kdefs>
 TNum * fast_Ring<TNum, kdefs>::createAdditionTable()
{

	unsigned short i, j, k , l;  // sollte nicht unb short sein, sondern vom basistyp abhaengen

	TNum * tAdditionTable=0;

	size_t tableSize = getMaxPairIndex() + 1;

	#ifdef DEBUG
	std::cerr << "createAdditionTable::tableSize = " << tableSize << std::endl;
	#endif

	tAdditionTable = new TNum[tableSize];

	for (i=0; i<getCharacteristic(); i++)
		for (j=0; j<((getCharacteristic()-1)*getEpsPrecision())+1; j++)
			for (k=0; k<getCharacteristic(); k++)
				for (l=0; l<((getCharacteristic()-1)*getEpsPrecision())+1; l++)
				{
					TNum 	z1(i, j);
					TNum 	z2(k, l);
					size_t index = getPairIndex(z1, z2);
					assert(index < 	tableSize && index>=0);

					tAdditionTable[index].setX   (
										  ( (int)z1.getX() + (int)z2.getX()  )
										% getCharacteristic() 
									);

					tAdditionTable[index].setEps ( 
										  ((int)z1.getEps() + (int)z2.getEps())
										% getCharacteristic() 
									);
				}
	return tAdditionTable;
}


/// only for epsilon==0
template <class TNum,class kdefs>
 TNum* fast_Ring<TNum,kdefs>::createFastAdditionTable()
{
	TNum * tfastAdditionTable=NULL;

	if (getEpsPrecision()==0)
	{

		tfastAdditionTable = new TNum[ getCharacteristic()*2 ];
	
		for (int m=0; m<getCharacteristic(); m++)
		{
			tfastAdditionTable[m]=m;	
		}
	
		for (int m=0; m<getCharacteristic(); m++)
		{
			tfastAdditionTable[m + getCharacteristic() ]=m;	
		}
	}
	return tfastAdditionTable;
}

/// @todo modulo % getCharacteristic() in initAddInv() unnoetig
template <class TNum,class kdefs>
 TNum* 	fast_Ring<TNum,kdefs>::createAdditiveInverseTable()
{
	int i, j; 

	size_t 	tableSize = getMaxSingleIndex() + 1;

	TNum*	tadditiveInverseTable = new TNum[ tableSize ];

	for (i=0; i<getCharacteristic(); i++)
		for (j=0; j<( ( getCharacteristic()-1 )*getEpsPrecision() )+1; j++)
		{
			TNum 	z1(i, j);

			size_t 	index = getSingleIndex( z1); 
			assert(index<tableSize && index>=0);

			TNum 	z2(	( getCharacteristic() - i) % getCharacteristic(),
					 (getCharacteristic() - j) % getCharacteristic()    );

			tadditiveInverseTable[ index ] = z2;
		}
	return tadditiveInverseTable;
}




/// creates multiplication table
template <class TNum,class kdefs>
TNum * fast_Ring<TNum, kdefs>::createMultiplicationTable()
{	
	int i, j, k, l; 

	TNum * tMultiplicationTable = NULL;
	
	size_t tableSize = getMaxPairIndex() + 1;

	tMultiplicationTable = new TNum[tableSize];

	for (i=0; i< getCharacteristic(); i++)
		for (j=0; j<( (getCharacteristic()-1)*getEpsPrecision() )+1; j++)
			for (k=0; k < getCharacteristic(); k++)
				for (l=0; l<( (getCharacteristic() - 1)* getEpsPrecision() ) + 1; l++)
				{
					TNum 	z1 (i,j);
					TNum 	z2 (k,l);

					size_t 	index = getPairIndex(z1, z2);

					assert(index < 	tableSize && index>=0 );

					TNum result;

					result.setX  (  	((int) z1.getX() * (int)z2.getX() ) 
								% getCharacteristic()
							);
			
					result.setEps( (   (int)z1.getX() * (int)z2.getEps() 
							       + (int)z2.getX() * (int)z1.getEps() 
						         ) % getCharacteristic()	
							);
				
					if (j==0 && l==0 && i==k )
					{
						if (result.getX()==getCharacteristic()-1)
						{
							bContainsImagNum_m=true;
							imagNum_m=z1;
						}
					}

					tMultiplicationTable[ index ] = result;
				}


	/*for (i=0; i< getCharacteristic(); i++)
			for (k=0; k < getCharacteristic(); k++)
				{
					TNum 	z1 (i,0);
					TNum 	z2 (k,0);

					size_t 	index = getPairIndex(z1, z2);

					assert(index < 	tableSize && index>=0 );

					TNum result;

					result.setX  (  	((int) z1.getX() * (int)z2.getX() ) 
								% getCharacteristic()
							);
			
					assert(tMultiplicationTable[ index ] == result);
				}*/

	return tMultiplicationTable;
}


/// creates multiplicative inverse table
/// 
/// @todo in den Initialisierungsfunktionen die Convert-Funktion nutzen!
/// @todo initInv überarbeiten!
template <class TNum,class kdefs>
 TNum* fast_Ring<TNum,kdefs>::createMultiplicativeInverseTable()
{
	int i, j, k, l; 
	
	size_t tableSize = getMaxSingleIndex() + 1;

	TNum * inverses1 = new TNum[ tableSize ];

	for (i=0; i<getCharacteristic(); i++)

		for (j=0; j<((getCharacteristic()-1)*getEpsPrecision())+1; j++)
		{
			TNum 	z1(i, j);
			size_t 	index = getSingleIndex( z1 );

			assert(index<tableSize && index>=0);

			inverses1[ index ]  = TNum::Zero; 

			for (k=0; k<getCharacteristic(); k++)
				for (l=0; l<((getCharacteristic()-1)*getEpsPrecision())+1; l++)
				{
					TNum 	z2(k, l);
					TNum 	z;
					z.setX(  ((int)z1.getX() * (int)z2.getX()) % getCharacteristic()  );

					z.setEps ( (  (int)z1.getX() * (int)z2.getEps() 
						     + (int)z2.getX() * (int)z1.getEps()  
						   ) % getCharacteristic()  );

					if ( z == TNum::One )
					{
						size_t index_2 = getSingleIndex( z1);
						assert( index_2<tableSize && index_2 >=0);
						inverses1 [index_2 ] = z2;
						break;
					}
				}
		}
	return inverses1;
}



	
/// @todo is only correct for epsPrecision==0 and does not really belong in fast_Ring class!
template <class TNum,class kdefs>
 TNum* fast_Ring<TNum,kdefs>::initElementsToExponentsTab(TNum erzeuger)
{
	TNum * tElementsToExponentsTab=NULL;
	
	size_t tableSize = getMaxSingleIndex() + 1;

	tElementsToExponentsTab = new TNum[tableSize];

	tElementsToExponentsTab[0]=TNum::Zero;

	TNum num=erzeuger;

	for (int i=1; i<getCharacteristic(); i++)
	{
		tElementsToExponentsTab[ num.getX() ] = i;
		num.setX( ( (int)num.getX() * (int)erzeuger.getX() ) % getCharacteristic());
	}
	if (num!=erzeuger)
	{
		std::cerr << "num " << num << std::endl;
		std::cerr << "generator " << erzeuger << std::endl;
		num = erzeuger;
		for (int i=1; i<getCharacteristic(); i++)
		{
			std::cerr << "num " << (int)num.getX() * (int)erzeuger.getX() << std::endl;
			num.setX( ( (int)num.getX() * (int)erzeuger.getX() ) % getCharacteristic());
			std::cerr << "num " << num << std::endl << std::endl;
		}

		assert( num==erzeuger );
	}
	return tElementsToExponentsTab;
}


/// @note exponentsToElementTab koennte PerformanceProbleme bereiten, wenn 
///  intern die Indizierungsfunktnion verwendet wuerde
template <class TNum,class kdefs>
 TNum* fast_Ring<TNum,kdefs>::initExponentsToElementTab(TNum erzeuger)
{

	TNum * 	tExponentsToElementTab=NULL;
	
	long 	tableSize=0;

	tableSize = getMaxSingleIndex() + 1;

	tExponentsToElementTab = new TNum[tableSize];

	tExponentsToElementTab[0] = TNum::Zero;

	TNum num = erzeuger;

	for (int i=1; i<getCharacteristic(); i++)
	{
		//exponentsToElementTab[ TNum::getSingleIndex( i ) ] = num;
		tExponentsToElementTab[  i  ] = num;
		num.setX( ( (int)num.getX() * (int)erzeuger.getX() ) % getCharacteristic() );
	}
	assert( num==erzeuger );
	return tExponentsToElementTab;
}


template <class TNum,class kdefs>
typename fast_Ring<TNum,kdefs>::FieldType::ElementType* fast_Ring<TNum,kdefs>::createModuloTable() 
{
	ScalarType * tModuloTable=NULL;
	

	 
	 //moduloTableSize_m = 9*characteristic*characteristic+characteristic;
	// 
	size_t 	_characteristic=characteristic;
	moduloTableSize_m = 15*_characteristic*_characteristic+_characteristic;
	tModuloTable= new ScalarType[ moduloTableSize_m ];

	long currPos=0;
	while( currPos < moduloTableSize_m )
	{
		for (long j=0; j<characteristic;j++  )
		{
			tModuloTable[ currPos ] = j;
			currPos++;
		}
	}
	
	return tModuloTable;
}


template <class TNum,class kdefs>
inline int  fast_Ring<TNum,kdefs>::getLookupModuloTableSize( ) const
{
	return moduloTableSize_m;
}


template <class TNum,class kdefs>
inline typename fast_Ring<TNum,kdefs>::FieldType::ElementType  fast_Ring<TNum,kdefs>::lookupModuloTable(int convertee) const
{
	#ifdef SAFE
		//if (convertee<0 )
		//	std::cerr << "error: convertee = "<<convertee << std::endl;
		//assert(convertee>=0 );
		//assert(convertee<moduloTableSize_m);
		//assert(convertee>=0 && convertee<moduloTableSize_m);
	#endif
	
	if (convertee>=moduloTableSize_m || convertee<0  )
	{
		//std::cerr << "convertee " << convertee << endl;
		return ConvertScalar(convertee);
	}
	//assert(convertee>=0 && convertee<moduloTableSize_m);
	return moduloTable[convertee];
}

 
template <class TNum,class kdefs>
typename fast_Ring<TNum,kdefs>::sqrtInf_t* fast_Ring<TNum,kdefs>::createSqrtTable()
{

	typename fast_Ring<TNum,kdefs>::sqrtInf_t * 	tSqrtTable=NULL;
	
	long 	tableSize=0;

	tableSize = getMaxSingleIndex() + 1;

	tSqrtTable = new typename fast_Ring<TNum,kdefs>::sqrtInf_t[tableSize];

	typename fast_Ring<TNum,kdefs>::sqrtInf_t entry(0,TNum::Zero);

	for (int i=0; i<getCharacteristic(); i++)
	{
		tSqrtTable[ i ] = entry;
	}

	for (int i=0; i<getCharacteristic(); i++)
	{
		int res= (i*i) % getCharacteristic();
		if (tSqrtTable[  res  ].solutions==0)
			tSqrtTable[  res  ].sqrt =TNum(i);
		else
		{
			assert(  getCharacteristic()-i ==tSqrtTable[  res  ].sqrt.getX() );
		}
		tSqrtTable[  res  ].solutions ++;
		
	}

	for (int i=0; i<getCharacteristic(); i++)
	{
		assert(	tSqrtTable[  i  ].solutions==0 ||	
			tSqrtTable[  i  ].solutions==2   ||  tSqrtTable[  i  ].sqrt.getX()==0 || getCharacteristic()==2 );
		
	//	std::cerr<< "tSqrtTable["<< i <<"    ].sqrt.getX()" << (int)tSqrtTable[  i  ].sqrt.getX() << std::endl;
		if (tSqrtTable[  i  ].solutions>0)
			assert(	(tSqrtTable[  i  ].sqrt.getX() *	tSqrtTable[  i  ].sqrt.getX() ) % getCharacteristic() == i	 );
		assert(	tSqrtTable[  i  ].sqrt.getEps()==0);
	}
 
	return tSqrtTable;
}




//----------------------Power---------------------------------------------------

template <class TNum,class kdefs>
inline TNum fast_Ring<TNum,kdefs>::pow(const TNum x,unsigned int exp) const
{
	#ifdef SAFE
		assert(exp>=0);
	#endif

	
	if (exp==1)
		return x;
	if (exp==0)
		return TNum::One;
	
	
	TNum erg = x;
	for (; exp >1; exp--)
	{
		multiplyInPlaceRef( erg, x ) ;			
	}
	return erg;
}



template <class TNum,class kdefs>
inline void fast_Ring<TNum,kdefs>::powInPlace( TNum & x, unsigned int exp) const
{
	#ifdef SAFE
		assert(exp>=0);
	#endif

	
	if (exp==1)
		return ;
	if (exp==0)
		x = TNum::One;
	
	TNum tmp = x;
	for (; exp >1; exp--)
	{
		multiplyInPlaceRef( x, tmp ) ;			
	}
	return;
}


//----------------------Index---------------------------------------------------

template <class TNum,class kdefs>		
inline  size_t 	fast_Ring<TNum,kdefs>::getMaxSingleIndex() const
{ 
	return TNum::getMaxSingleIndex( getCharacteristic());
}



template <class TNum,class kdefs>		
inline  size_t 	fast_Ring<TNum,kdefs>::getMaxPairIndex() const
{ 
	return TNum::getMaxPairIndex( getCharacteristic() );
}



template <class TNum,class kdefs>		
inline  size_t 	fast_Ring<TNum,kdefs>::getSingleIndex(const TNum z1) const
{ 
	// 1. reicht das, wenn ich hier TNum::getPairIndex(z1,z2) einsetze,
	// oder muss ich ueberall TNum::getxxxIndex verwenden damit das Programm schnell laeuft?
	return TNum::getSingleIndex(z1, getCharacteristic());
}

template <class TNum,class kdefs>		
inline  size_t 	fast_Ring<TNum,kdefs>::getSingleIndexByRef(const TNum &z1) const
{ 
	return TNum::getSingleIndexByRef(z1, getCharacteristicRef());
}

template <class TNum,class kdefs>
inline  size_t 	fast_Ring<TNum,kdefs>::getPairIndex(const TNum z1, const TNum z2)  const
{ 
	return TNum::getPairIndex(z1, z2, getCharacteristic() );			
}


template <class TNum,class kdefs>
inline  size_t 	fast_Ring<TNum,kdefs>::getPairIndexByRef(const TNum & z1, const TNum & z2) const
{ 
	return   TNum::getPairIndexByRef(z1, z2, getCharacteristicRef() ) ;
}



