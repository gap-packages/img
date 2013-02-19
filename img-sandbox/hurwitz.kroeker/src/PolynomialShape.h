#pragma once

#include <vector>
#include "FactorPolynomialWrapper.h"

// ShapeRepType expected as  ' vector<INTType> '. Improvement:: could accept a bidirectional iterator.
template<typename ShapeRepType>
inline ShapeRepType    dualPartition(const ShapeRepType dualmulstruct  )
{
    int max=0;
    ShapeRepType    multiplicityStructure(max,0);

    // max =* ( std::max_element( dualmulstruct.begin(), dualmulstruct.end() ) );
    for (int pos= dualmulstruct.size()-1; pos >=0; pos--)
    {
        if (dualmulstruct[pos]>max)
            max=dualmulstruct[pos];
    }


    if (dualmulstruct.size()>0)
    {
        multiplicityStructure.resize(max);

        for (int i=0; i<max; i++)
        {
            
            for (int pos= dualmulstruct.size()-1; pos >=0; pos--)
            {
                if ( dualmulstruct[ pos ] > i )
                    multiplicityStructure[ i ] ++;
            }
            
            
        }
    }   
    return  multiplicityStructure;
}

/*
internalPolynomialMatchesShape (Shape,RingElement,RingElement) := (shape, univarPoly, variable ) -> (
    assert( #(degree univarPoly) == 1); -- otherwise we have a multigraded ring.
    assert( ring univarPoly === ring variable );
    dualMulStructSize := #(shape#"dual");
    gcdF := univarPoly;
    dF := univarPoly;
    pos := 0; -- position in 'shape#"dual"'  array
    --
    -- the while loop is used instead of apply because it is not possible to return from an apply loop
    while (gcdF != 1_(ring variable) and dF != 0_(ring variable)) do (
        prevDegree := degree(variable, gcdF);
        dF = diff( variable, dF );
        gcdF =  gcd( gcdF, dF );
        referenceDegreeDifference := 0;
        if ( pos < dualMulStructSize) then 
            referenceDegreeDifference = (shape#"dual")_pos;
        if  (prevDegree - degree(variable, gcdF) != referenceDegreeDifference)    then 
            return false;
        pos = pos+1;
    );
    -- if (pos >= dualMulStructSize) then the multiplicity structure was fulfilled:
    return (pos >= dualMulStructSize );
)

*/

/// @todo multiplicityStructure as object
template <class TPolynomialRing, typename Shape>
inline  bool        polynomialMatchesShapeStrict(const typename TPolynomialRing::Element     & poly, 
                                    const TPolynomialRing & polynomialRing, const Shape & required_shape)
{
     Shape shape = RationalMapSearch::FLINTFactorPolynomial::computeShape<TPolynomialRing>(poly,polynomialRing);
    
     return ( shape.getShapeRep() == required_shape.getShapeRep() );
}


/// @todo multiplicityStructure as object
template <class TPolynomialRing, typename Shape>
inline  bool        polynomialMatchesShape(const typename TPolynomialRing::Element     & poly, 
                                    const TPolynomialRing & polynomialRing, const Shape & shape)
{

    typedef typename TPolynomialRing::Element       TPolynomial;

    TPolynomial gcdf( poly);
    
    //TPolynomial dF = polynomialRing.diff( gcdf );    
    TPolynomial dF(poly);

    int prevDegree, referenceDegreeDifference;
    size_t pos=0;

    typename Shape::ShapeRepType  dualShapeRep =  shape.getDualShapeRep();
    
    while (! gcdf.isConstant() && !(dF.isZero())  )
    {
        prevDegree = gcdf.getExactDegree();
        //gcdf = polynomialRing.gcd(gcdf,dF);
        dF = polynomialRing.diff( dF );
        gcdf = polynomialRing.fastgcd(gcdf,dF);
        referenceDegreeDifference = 0;
        
        if ( pos < dualShapeRep.size() )   
            referenceDegreeDifference = dualShapeRep[pos];
        if  ( ( prevDegree - gcdf.getExactDegree()) != referenceDegreeDifference )
            return false;
        pos = pos+1;
    }
     return (pos >= dualShapeRep.size() );
}


/// @todo multiplizitätsstruktur wird manchmal immer noch nicht korrekt berechbet ->
/// @todo multiplicityStructure as object
template <class TPolynomialRing, typename ShapeRepType>
inline  ShapeRepType        computeMultiplicityStructureStrict(const typename TPolynomialRing::Element     & poly, 
                                    const TPolynomialRing & polynomialRing)
{

   //factoring is c.a 4 times slower and maybe incorrect
    return RationalMapSearch::FLINTFactorPolynomial::computeShape<TPolynomialRing>(poly,polynomialRing);
}


/// @todo multiplizitätsstruktur wird manchmal immer noch nicht korrekt berechbet ->
/// @todo multiplicityStructure as object
template <class TPolynomialRing, typename ShapeRepType>
inline  ShapeRepType        computeMultiplicityStructure(const typename TPolynomialRing::Element     & poly, 
                                    const TPolynomialRing & polynomialRing)
{

    typedef typename TPolynomialRing::Element       TPolynomial;
     
    ShapeRepType    degrees;

    if ( poly.isConstant() )
        return degrees;

    TPolynomial gcdf( poly);
    
    TPolynomial dF= polynomialRing.diff( gcdf );    

    degrees.push_back( gcdf.getExactDegree() );
    
    while (! gcdf.isConstant() && !(dF.isZero())  )
    {
        //gcdf = polynomialRing.gcd(gcdf,dF);
        gcdf = polynomialRing.fastgcd(gcdf,dF);
        degrees.push_back(gcdf.getExactDegree() );
        dF = polynomialRing.diff( dF );
    }

    ShapeRepType    dualMultiplicityStructure;
    dualMultiplicityStructure.resize( degrees.size()-1  ,0);

    for (int pos= degrees.size()-2; pos>=0; pos--)
    {
        dualMultiplicityStructure[pos]= degrees[pos]- degrees[pos+1]  ;
    }

    return dualPartition( dualMultiplicityStructure );
}
