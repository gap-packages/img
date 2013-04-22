


#include <assert.h>
#include <cstdio>
#include <iostream>
#include <functional>
#include <numeric>
#include <algorithm>
#include <sstream>

/*#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/foreach.hpp>
#include <boost/bind.hpp> 
#include <boost/lambda/if.hpp>
*/

//using namespace boost;
//using namespace boost::lambda;

 


#include "Shape.h"


namespace RationalMapSearch
{


     std::ostream &  operator<<(std::ostream & out, const PolynomialFactorBluePrint& bp)
     {
            out << "( [" << bp.polynomialId_m << "] : " << bp.degree_m << "^" << bp.multiplicity_m << ")" ;
            return out;
     }

     std::ostream &  operator<<(std::ostream & out, const Shape& shape)
     {
            out << "( ";
            for (size_t i=0; i<shape.size(); i++)
                out << shape[i] << " " ;
            out << " ) ";
            return out;
     }

    // todo: test, ob es mit liste und vector ok ist.
     Shape::Shape( const    Shape& shape ) :  shape_m ( shape.shape_m ) , 
                                                    dual_m  ( shape.dual_m  ) , 
                                                    degree_m( shape.degree_m),
                                                    multiplicityDegreeMap_m( shape.multiplicityDegreeMap_m )
    {
    }


    Shape::Shape( const    ShapeRepType& vec ) :  shape_m( initShapeRep(vec) ),
                                                    dual_m( conjugate(vec) ),
                                                    degree_m ( std::accumulate( vec.begin(), vec.end(), 0 )),
                                                    multiplicityDegreeMap_m( createMultiplicityDegreeRep(vec) )
    {
             
    }

    
    

    // todo: hast vergessen, dass wenn ein shape nur unendlich hat, nix mehr übrig bleibt.
    Shape::ShapeRepType Shape::conjugate(const    ShapeRepType& vec)
    {
       
        ShapeRepType sortedVec = vec; 
        std::sort( sortedVec.begin(), sortedVec.end(), std::greater<int>() );

        // Macaulay2: dualPartition = L -> apply(max L, i-> #(select(L,j->j>i)))


        ScalarType max= 0;
        if ( sortedVec.size()>0 )
            max = sortedVec[0];

        ShapeRepType conjugateVec( max );
        for (int pos=0; pos<max ; pos++)
        {   
            //boost::function<bool(ScalarType)> f = ( boost::lambda::_1 > pos );
            //// sortedVec[pos] = (ScalarType) count_if ( sortedVec.begin(), sortedVec.end(), [](ScalarType x) { return x > pos; });
            //// conjugateVec[pos] = (ScalarType)  std::count_if ( sortedVec.begin(), sortedVec.end(), f );
            //conjugateVec[pos] = (ScalarType)  std::count_if ( sortedVec.begin(), sortedVec.end(),  boost::lambda::_1 > pos  );
            int counter=0;
            for  (ShapeRepType::const_iterator it= sortedVec.begin(); it!= sortedVec.end(); it++)
            {
                if  ((*it)>pos )
                    counter++;
            }
            conjugateVec[pos] = counter;
        }
        return conjugateVec;
    }

    bool  Shape::hasExponent(ScalarType   exp) const
    {
        for (size_t pos =0;pos!=shape_m.size();pos++)
        {
            if (  shape_m[pos] == exp  )
                return true;
        }
        return false;
    }

    //removes given exponent from the Shape, if possible and returns the result as new shape.
    Shape   Shape::removeExponent(ScalarType   exp) const
    {
        if ( ! hasExponent(exp))
            return Shape(*this);
        
        bool first=true;
       //  boost::function<bool(ScalarType)> expMatches = (boost::lambda::if_then_else_return ((_1==exp && first), !first=false,  false ));
       // boost::remove_copy_if (shape_m.begin(),shape_m.end(),result,expMatches);
        
        ShapeRepType result;
        for (size_t pos =0; pos!=shape_m.size(); pos++)
        {
            if (  shape_m[pos] == exp && first )
            {
                first=false;
                continue;
            }
            result.push_back( shape_m[pos] );
        }
        return Shape(result);
    }
        
    bool  Shape::operator==(const Shape &shape) const
    {
        return (shape.shape_m==this->shape_m);
    }

    Shape::ShapeRepType Shape::initShapeRep(const    ShapeRepType& vec)
    {
        ShapeRepType retval =vec; 
        std::sort( retval.begin(),retval.end(), std::greater<int>() );
        ShapeRepType::iterator it;
        for (it= retval.begin(); it != retval.end();it++)
        {
            assert((*it)>0);
        }
        return retval;
    }

    Shape::MultiplicityDegreeHashType  Shape::createMultiplicityDegreeRep(const Shape::ShapeRepType &ShapeRepType) 
    {
        //std::vector<std::pair> result;
    
        Shape::MultiplicityDegreeHashType multiplicityDegreeMap;

        for ( size_t pos =0; pos!=ShapeRepType.size(); pos++)
        {

            if ( multiplicityDegreeMap.find( ShapeRepType[pos] ) == multiplicityDegreeMap.end() )
                multiplicityDegreeMap.insert(std::pair< MultiplicityDegreeRepKey, MultiplicityDegreeRepValue > (ShapeRepType[pos],1) );
            else 
                multiplicityDegreeMap[ ShapeRepType[pos] ]++;
        }

        return multiplicityDegreeMap;
    }

    bool    Shape::hasNaturalNormalizableFactor() const
    {
        Shape::MultiplicityDegreeHashType::const_iterator it;
        for ( it = multiplicityDegreeMap_m.begin(); it != multiplicityDegreeMap_m.end(); it++ )
        {   
                if ((*it).second==1)
                    return true;
        }
        return false;
    }

    // Todo: move computation to constructor and a member variable? 
    Shape::ScalarType   Shape::getMaxFactorDegree() const
    {
        Shape::MultiplicityDegreeHashType::const_iterator it;
        Shape::ScalarType maxFactorDegree = 1;
        for ( it = multiplicityDegreeMap_m.begin(); it != multiplicityDegreeMap_m.end(); it++ )
        {   
                if ((*it).second>maxFactorDegree)
                    maxFactorDegree=(*it).second;
        }
        return maxFactorDegree;
    }

    std::string     Shape::toString()   const
    {
        std::stringstream strstr;
        
        strstr << "{ " ;
        for (size_t pos=0; pos <shape_m.size(); pos++)
        {
            if (pos>0 )
                strstr << ", " ;
            strstr << shape_m[pos];
        }
        strstr << "} " ;
        return strstr.str();
    }

    void Shape::test()
    {
   
        // removed initializer list usage for backward compatibility
        //        std::vector< Shape::ScalarType >    partition = { 4,3,2,2,2 };

        // std::vector< Shape::ScalarType >    dualPartition =  { 5,5,2,1 };

        int ar[]={ 4,3,2,2,2 };
        const int TotalItems = sizeof(ar)/sizeof(ar[0]);
        std::vector< Shape::ScalarType >    partition(ar, ar+TotalItems);

        int dar[]= { 5,5,2,1 };
        const int TotalItemsDar = sizeof(dar)/sizeof(dar[0]);
        std::vector< Shape::ScalarType >    dualPartition(dar, dar+TotalItemsDar);

        
        Shape shape(partition);
        Shape::ShapeRepType dual = shape.getDualShapeRep();
        Shape::ShapeRepType storedPartition = shape.getShapeRep();

        assert(shape.hasNaturalNormalizableFactor() );
        assert(shape.conjugate(shape.conjugate(storedPartition))==storedPartition);

        // copy(dual.begin(), dual.end(), std::ostream_iterator<Shape::ScalarType>(std::cout, "\n"));
        assert( dualPartition==dual );
        assert( partition == storedPartition );
        assert( shape.getDegree()==13 );
        assert( shape.getMaxFactorDegree()==3 );

        Shape reducedShape = shape.removeExponent(2);

        //Shape refReducedShape = { 4,3,2,2 };

        int  redar[]={ 4,3,2,2 };
        const int TotalItemsRed = sizeof(redar)/sizeof(redar[0]);
        std::vector< Shape::ScalarType >    preRefReducedShape(redar, redar+TotalItemsRed);
        Shape refReducedShape(preRefReducedShape);

        assert( refReducedShape==refReducedShape );      

        //Shape shapeN = { 4,4,3,3,2,2,2 };

        int  Nar[]={ 4,3,2,2 };
        const int TotalItemsNar = sizeof(Nar)/sizeof(Nar[0]);
        std::vector< Shape::ScalarType >    preShapeN(Nar, Nar+TotalItemsNar);
        Shape shapeN(preShapeN);

        assert(! shapeN.hasNaturalNormalizableFactor() );

        MultiplicityDegreeHashType::const_iterator it;

        assert(shapeN.multiplicityDegreeMap_m[4]==2);
        assert(shapeN.multiplicityDegreeMap_m[3]==2);
        assert(shapeN.multiplicityDegreeMap_m[2]==3);

        for (it=shapeN.multiplicityDegreeMap_m.begin();it!=shapeN.multiplicityDegreeMap_m.end();it++)
        {
            std::cerr << (*it).first << "^" << (*it).second << std::endl;
        }
        std::cerr << "Shape test succeeded!\n";
    }


    ShapeList::ShapeList( const  std::vector<Shape> shapeList )    : shapeList_m(shapeList), reordered_m(false)
    {
        checkConsistency();   
    }
    ShapeList::ShapeList( const ShapeList& shapeList): shapeList_m(shapeList.shapeList_m),
                                                                reordered_m(shapeList.reordered_m)
    {   
       checkConsistency();
    } 

    ShapeList   ShapeList::removeExponent( unsigned int polynomId, ScalarType   exp) const
    {
            assert( size()>polynomId ) ;
            ShapeList slcopy = *this;
            slcopy.shapeList_m[polynomId] = slcopy.shapeList_m[polynomId].removeExponent(exp);
            return slcopy;
    }

    void ShapeList::checkConsistency() const
    {
        assert( shapeList_m.size()>1 );
        if (reordered_m) 
        {
            assert(reorderMap_m.size()==shapeList_m.size() );
            for (size_t pos=0;pos<reorderMap_m.size(); pos++)
            {   
                assert(1==(int) count ( reorderMap_m.begin() , reorderMap_m.end(), pos ));
            }
        }

        size_t polDegree = shapeList_m[0].getDegree();

        for (ShapeListRepType::const_iterator it=shapeList_m.begin(); it!=shapeList_m.end(); it++)
        {
            assert( (*it).getDegree() == polDegree );
        }
        /*BOOST_FOREACH( Shape shape, shapeList_m )
        {
            assert( shape.getDegree() == polDegree );
        }*/
        //boost::function<void(Shape)>  checkDegree = ( assert( (boost::bind( &Shape::getDegree,boost::lambda::_1) ==polDegree ))() ) ;
        //std::for_each(shapeList_m.begin(),shapeList_m.end(),checkDegree  );
        //std::for_each(shapeList_m.begin(),shapeList_m.end(), boost::bind( &Shape::getDegree,boost::lambda::_1) ==polDegree  );
    
    }

     //       ShapeList( const ShapeList&,  std::vector<size_t> reorderMap );

    ShapeList::ScalarType  ShapeList::getMaxFactorDegree( ) const
    {
        ShapeList::ScalarType maxFactorDegree = 1;
       
      /* BOOST_FOREACH( Shape shape, shapeList_m )
        {
            ShapeList::ScalarType tmpFactor = shape.getMaxFactorDegree();
            if ( tmpFactor >maxFactorDegree)
                maxFactorDegree = tmpFactor;
        }*/

         for (ShapeListRepType::const_iterator it=shapeList_m.begin(); it!=shapeList_m.end(); it++)
        {
              ShapeList::ScalarType tmpFactor = (*it).getMaxFactorDegree();
                if ( tmpFactor >maxFactorDegree)
                     maxFactorDegree = tmpFactor;
        }

        return maxFactorDegree;
    }

    
    ShapeList::ScalarType  ShapeList::getConstructionMaxFactorDegree( ) const
    {
        
        ShapeList::ScalarType maxFactorDegree = 1;
                
        if (  shapeList_m[0].getMaxFactorDegree()>maxFactorDegree)
                maxFactorDegree = shapeList_m[0].getMaxFactorDegree();
        if (  shapeList_m[1].getMaxFactorDegree()>maxFactorDegree)
                maxFactorDegree = shapeList_m[1].getMaxFactorDegree();
        
        return maxFactorDegree;
    }

    /// @TODO BUG! - what does mean BUG, and why the hell maxDegree was initially = 2 ?
    ShapeList::ScalarType  ShapeList::getDegree( ) const
    {
        int maxDegree = 1;
       
        /* BOOST_FOREACH( Shape shape, shapeList_m )
        {
            if (shape.getDegree()>maxDegree) 
                maxDegree = shape.getDegree()  ;
        }
        */

        for (ShapeListRepType::const_iterator it=shapeList_m.begin(); it!=shapeList_m.end(); it++)
        {
            if ((*it).getDegree()>maxDegree) 
                maxDegree = (*it).getDegree()  ;
        }

        return maxDegree;
    }

    /// lower characteristic bound for 'polynomialMatchesShape'-check
    /// @todo: need to prove correctness - not correct (bound too low) because all critival values obviously also have to be distinct, 
    /// and  also  factors where there are no other factors in the same polynomial with same multiplicity ( they also have to be distinct ).
    ///  minCharacteristic is  Max(computeLowerCharacBound, #cv, #factors where there are no other factors in the same polynomial with same multiplicity)
    /// of course, then lowerCharBound will still not be fully correct, but much more accurate. 
    /// 
    /// @note first two shapes are leaved out, because they are constructed correct a priori by the search algorithm.
    /// @note the function really depends on search algorithm and multiplicity check ('polynomialMatchesShape') implementations.
    
    ShapeList::ScalarType  ShapeList::computeLowerCharacBound( ) const
    {   
        ScalarType lowerCharBound = 2;
        
        // for each requires gcc 4.6
        //for (Shape shape : shapeList_m)
    
        //    BOOST_FOREACH( Shape shape, shapeList_m )
        for (size_t pos=2;pos< shapeList_m.size(); pos++)
        {
            if (shapeList_m[pos].getMaxExponent()>lowerCharBound-1) 
                lowerCharBound = shapeList_m[pos].getMaxExponent() + 1;
        }
        return lowerCharBound;
    }
    
    ShapeList::ScalarType  ShapeList::getMaxExponent() const
    {
        ScalarType maxExponent = 2;
        /*BOOST_FOREACH( Shape shape, shapeList_m )
        {
            if (shape.getMaxExponent()>maxExponent) 
                maxExponent = shape.getMaxExponent()  ;
        }*/

        for (ShapeListRepType::const_iterator it=shapeList_m.begin(); it!=shapeList_m.end(); it++)
        {
            if ((*it).getMaxExponent()>maxExponent) 
                maxExponent = (*it).getMaxExponent()  ;
        }

        return maxExponent;
    }
    

    const Shape&   ShapeList::operator[](size_t pos) const
    {
        if (reordered_m) 
        {           
            return  shapeList_m[reorderMap_m[pos]];
        }
        return shapeList_m[pos];
    }

    void ShapeList::test()
    {   
        //Shape shape = { 4,3,2,2,2 };

        int ar[]={ 4,3,2,2,2 };
        const int TotalItems = sizeof(ar)/sizeof(ar[0]);
        std::vector< Shape::ScalarType >    partition(ar, ar+TotalItems);
    
        Shape shape(partition);
        
        std::vector< Shape >  preShapeList ( 3, shape );
        
        ShapeList   shapeList(preShapeList);
        assert(shapeList.computeLowerCharacBound()>4);

        std::cerr << "shapeList.getMaxFactorDegree()" << shapeList.getMaxFactorDegree()<<std::endl;
        assert(shapeList.getMaxFactorDegree()==3);

         Shape reducedShape= shapeList[0].removeExponent(2) ;

        Shape reducedShape2 (reducedShape);

        reducedShape2=reducedShape;

        preShapeList.push_back(reducedShape );
        //ShapeList   illegalShapeList(preShapeList);

    }

    
    

}















