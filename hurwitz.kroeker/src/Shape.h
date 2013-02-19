
#pragma once

#include <vector>
#include <assert.h>
#include <cstdlib>
#include <map>
#include <ext/hash_map> //TODO: replace with <unordered_map>!
#include <hash_map>
#include <iostream>
#include <functional>
#include <numeric>

namespace RationalMapSearch
{


    struct PolynomialFactorBluePrint
    {
        uint multiplicity_m;
        uint degree_m;   
        uint polynomialId_m;   
        PolynomialFactorBluePrint  (uint multiplicity, uint degree,uint polynomialId) : multiplicity_m(multiplicity), 
                                                                                degree_m(degree),
                                                                                polynomialId_m(polynomialId)
        {};

        PolynomialFactorBluePrint()
        {
            //assert(false);
        }

        static int  minDegreeMaxMultiplicityLower(const PolynomialFactorBluePrint &a, const PolynomialFactorBluePrint &b)
        {
            if ( a.degree_m != b.degree_m )
                return ( a.multiplicity_m < b.multiplicity_m );
            return ( a.degree_m < b.degree_m );
        }
    };

    std::ostream &  operator<<(std::ostream & out, const PolynomialFactorBluePrint& bp);

  

    struct MultiplicityDegreePair
    {
        uint multiplicity_m;
        uint degree_m;   

        MultiplicityDegreePair  (uint multiplicity, uint degree) : multiplicity_m(multiplicity), degree_m(degree)
        {};
        
    };

    ///@note though it is possible to parametrize Shape with 

    class Shape
    {
    
        public:
            typedef int ScalarType;

      
        typedef int32_t     MultiplicityDegreeRepKey;
        typedef int32_t     MultiplicityDegreeRepValue;

        typedef std::map< MultiplicityDegreeRepKey, MultiplicityDegreeRepValue >  MultiplicityDegreeHashType;
        
            typedef  std::vector< ScalarType > ShapeRepType;
        
        private:

            // problem: can't be const when using assignment...
            //const ShapeRepType    shape_m;
            //const ShapeRepType    dual_m;

          // problem: can't be const when using assignment...
            ShapeRepType    shape_m;
            ShapeRepType    dual_m;

        
            ScalarType  degree_m;
        
            MultiplicityDegreeHashType multiplicityDegreeMap_m;

    
        public:
           
          /*  outcommented initializer lists for backward compatibility
             template <class E>
            inline Shape(std::initializer_list<E> s);
          */
           

           Shape ( const    Shape& shape );

            Shape(const     ShapeRepType & vec);


            // to use std::vector:  
            //This means I must define a operator= with unclear semantics just because it is never used. Oh my noodles! 
            // see for reference http://blog.copton.net/archives/2007/10/13/stdvector/index.html
          /*  Shape& operator=(const Shape& rhs) { 
                std::cerr<< "abort in shape";
                abort();
                return *this;
            }*/

    
            // probably not this way, but additional class ShapeParser, which takes a string and returns a vector
            // template TVec   Shape(const std::string     &   str);
        
            //todo: degree 
            ScalarType  getDegree() const       {   return degree_m;    }        

            
            ShapeRepType    getShapeRep() const     {   return shape_m;    }      

            ShapeRepType    getDualShapeRep() const {   return dual_m;    }

            ScalarType  getMaxExponent() const  { assert(shape_m.size()>0);  return shape_m[0]; }

    
            /// todo: gibt es eigentlich eine Begründung für den DatentypWahl ScalarType für getMaxFactorDegree?
            ScalarType  getMaxFactorDegree() const;
            

            const ScalarType&  operator[](size_t pos) const { assert(pos<shape_m.size() ); return shape_m[pos] ; }
          
            size_t      size()   const          { return shape_m.size();  };

            
           const  ShapeRepType &   getShapeRepRef() const     {   return shape_m;    }      

            const ShapeRepType &   getDualShapeRepRef() const {   return dual_m;    }            

            
            Shape   removeExponent(ScalarType   exp) const;

            bool   hasExponent(ScalarType   exp) const;
            
            static MultiplicityDegreeHashType  createMultiplicityDegreeRep(const     ShapeRepType & vec);

            MultiplicityDegreeHashType  getMultiplicityDegreeRep() const  {return this->multiplicityDegreeMap_m;};

            bool    hasNaturalNormalizableFactor() const;

    
            static ShapeRepType         conjugate(const    ShapeRepType& vec);

            static ShapeRepType     initShapeRep(const    ShapeRepType& vec)  ; 

            std::string         toString()  const;

            static void test();

            // sortExponentsByDegree      
            
            bool operator==(const Shape &shape) const;
    };

   /*  outcommented initializer lists for backward compatibility
    template <class E>
    inline Shape::Shape(std::initializer_list<E> s):                  shape_m (s.begin(),s.end()),
                                    dual_m( conjugate(shape_m) ),
                                                    degree_m ( std::accumulate( shape_m.begin(), shape_m.end(), 0 )),
                                                    multiplicityDegreeMap_m( createMultiplicityDegreeRep(shape_m) )
        {
                
                
            }
*/

  std::ostream &  operator<<(std::ostream & out, const Shape& shape);

    class ShapeList
    {
        public:
           typedef  std::vector<Shape> ShapeListRepType;

        private:
         
            ShapeListRepType  shapeList_m;
         

            std::vector<size_t>    reorderMap_m;
            bool reordered_m;


            void    checkConsistency() const;

        public:
            typedef Shape::ScalarType ScalarType;

         
            ShapeList   removeExponent( unsigned int polynomId, ScalarType   exp) const;
    
    
            ShapeList( const  std::vector<Shape> shapeList );
            ShapeList( const ShapeList& );
            ShapeList( const ShapeList&,  std::vector<size_t> reorderMap );

        /*     outcommented initializer lists for backward compatibility
            template <class E>
            inline ShapeList(std::initializer_list<E> s);

           template <class E>
            inline ShapeList(std::initializer_list<E> s1,std::initializer_list<E> s2,std::initializer_list<E> s3);


            template <class E>
            inline ShapeList(std::initializer_list< std::initializer_list<E> > s);

            template <class E>
            inline  ShapeList&  operator=(std::initializer_list<E> s);

*/

            ScalarType  getDegree( ) const;

            ScalarType  getMaxFactorDegree( ) const;
            ScalarType  getMaxExponent() const;

            ScalarType  getCharLowerBound() const;
            ScalarType  getConstructionMaxFactorDegree()    const;
            ScalarType  computeLowerCharacBound( ) const;
        
            const Shape&   operator[](size_t pos) const;

            size_t  size() const {  return shapeList_m.size();  }
            
            // reorderShapeList ()
            // dropReordering()
            // getReorderMap();
            // std::vector<size_t> computeOptimalReorderingMap() const; // todo: introduce parameters ...
            // findMinPenaltyShapePos
    
            static void test();
    };

/*  outcommented initializer lists for backward compatibility
        template <class E>
        inline  ShapeList::ShapeList(   std::initializer_list<E> s1,
                                        std::initializer_list<E> s2,
                                        std::initializer_list<E> s3)     :  shapeList_m ( {s1,s2,s3} ),  
                                                                            reordered_m(false)

        {
            checkConsistency();
        }
        
        template <class E>
        inline ShapeList::ShapeList( std::initializer_list<E> s ):        shapeList_m (s.begin(),s.end()),  
                                                                        reordered_m(false)
        {
            checkConsistency();
        }

        template <class E>
        inline ShapeList::ShapeList( std::initializer_list< std::initializer_list<E> > s ):   shapeList_m (s.begin(),s.end()), 
                                                                                            reordered_m(false)
        {
            checkConsistency();
        }
     */
        
}

