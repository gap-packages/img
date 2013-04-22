#pragma once

#include <algorithm>

#include <iostream>

#include <ostream>

#include <vector>


extern "C" {
#undef __cplusplus
    # include "nmod_poly.h"
#define __cplusplus
}



namespace RationalMapSearch
{

    template <typename TPolRingType>
    class Power
    {

        public:
            typedef typename TPolRingType::Element BaseType;
            

        private:
            


            BaseType base_m;
            int exponent_m;

        public:

             Power(const typename  TPolRingType::Element &el  ) :base_m(el),exponent_m(1)            
            {
            }

            Power(const typename  TPolRingType::Element &el , int exponent) :base_m(el),exponent_m(exponent)            
            {
            }

            Power(const Power& pow) :base_m(pow.base_m),exponent_m(pow.exponent_m)            
            {
            }

            /*Power& operator = (const Power& pow)
            {
                if (&pow==this)
                    return *this;
                 BaseType base=pow.base_m;
                int exponent=pow.exponent_m;
                base_m=base;
                exponent_m=exponent;
                 return *this;
            }*/

            BaseType getBase() const 
            {
                return base_m;
            }

            BaseType first() const
            {   return getBase(); }

          
            int getExponent() const 
            {
                return exponent_m;
            }

            int second() const
            {   return getExponent(); }
    };

    template <typename TPower>
    class Product
    {
        public:
            typedef TPower FactorType;

        private:
            std::vector< TPower >  factors_m;

        public:
            //Product() {};

            void push_back(const TPower & pow )
            {
                factors_m.push_back(pow);
            }

            std::vector< TPower > getFactors()  const
            {
                return factors_m;
            }
        
            const TPower operator[](size_t pos)   const
            {   
                return factors_m.at(pos);
            }
 
            TPower operator[](size_t pos)  
            {   
                return factors_m.at(pos);
            }

            size_t size() const {return factors_m.size(); } 
        
    };

   template <typename PolRingType>
    class CIFactorPolynomial
    {
        public:
           typedef std::vector<typename  PolRingType::Element > (* IFactorPolynomial)(const typename PolRingType::Element & , const PolRingType &  );

            typedef Product<  Power<  PolRingType>  > (* IExtendedFactorPolynomial)(const typename PolRingType::Element & , const PolRingType &  );  
    };

  

    class FLINTFactorPolynomial
    {

        public:

            static const bool threadsafe ;

            //typedef typename IrreduciblePolTableType::PolRingType    TPolRing;
            //typedef typename TPolRing::ElementType  VecElemType;
           
            // ok, it seems that factorPolynomial will fail for non-square-free polynomials 
            // a check is also : multiply all factors and then it should be the original polynomial.
            template <typename PolRingType>
            static std::vector<typename  PolRingType::Element > factorPolynomial(const typename  PolRingType::Element & pol, const PolRingType & polRing)  
            {
 
                    typedef typename  PolRingType::ElementType  VecElemType;

                    int degree = pol.getExactDegree();

                      // pol.getCoeffRing()
                    // check if is Zero  pol(a= for a in 0...char coeffring
                    const typename PolRingType::CoeffRingType & coeffRing = polRing.getCoeffRing();

                    nmod_poly_t x ;
                    nmod_poly_init(x ,coeffRing.getCharacteristic() );
    
    
                    // riscy: what if pol.getCoeff is not nonnegative?
                    for (int currDegree= 0; currDegree<=degree; currDegree++)
                    {
                         assert( pol.getCoeff( currDegree).getX()>=0 );
                         nmod_poly_set_coeff_ui(x,   currDegree,  pol.getCoeff( currDegree).getX() );
                    }
                    nmod_poly_factor_t factors;
                    nmod_poly_factor_init( factors );
                    nmod_poly_factor ( factors, x );
                    
                    std::vector<typename  PolRingType::Element >    factorList;
                    for (long currFactorPos = 0; currFactorPos< factors[0].num_factors; currFactorPos++)
                    {   
                        //std::cerr << " currFactorPos " << currFactorPos << std::endl;
                        const nmod_poly_struct * currFactor = factors->factors[currFactorPos];
                        long currDegree = nmod_poly_degree (currFactor);
                        assert( currDegree>=0 );
                        typename PolRingType::ElementType    pol=   typename PolRingType::ElementType(currDegree);
                        for (int currExp=0;currExp<=currDegree;currExp++)
                        {
                          //  std::cerr << " currExp " << currExp << std::endl;
                            
                            ulong currCoeff= nmod_poly_get_coeff_ui(currFactor,currExp);

                            //std::cerr << " currCoeff " << currCoeff << std::endl;
                            pol.setCoeff(currExp, polRing.getCoeffRing().Convert(currCoeff));
                        }
                        factorList.push_back(pol);
                    }

                    nmod_poly_clear (x);
                    nmod_poly_factor_clear (factors);

                    return factorList ;
           };

            template <typename PolRingType>
            static std::vector<int > computeShape(const typename  PolRingType::Element & pol, const PolRingType & polRing) 
            {

                std::vector<int > shape;
                
                 std::vector<typename  PolRingType::Element > elements = factorPolynomial(pol,polRing);

                for (size_t pos=0;pos<elements.size(); pos++)
                {
                    int exponent=0;
                    std::pair<typename  PolRingType::Element,typename  PolRingType::Element> divisionResult;
                    divisionResult.first = pol;
                    do 
                    
                    {    
                         divisionResult = polRing.divide( divisionResult.first, elements[pos] );
                        if (not divisionResult.second.isZero() )
                            break;
                        exponent++;
                        
                    }
                    while ( divisionResult.second.isZero() );
                    assert( exponent>0 );
                    for (int deg=1; deg<=elements[pos].getExactDegree(); deg++)
                         shape.push_back( exponent );
                    
                }

                 return shape;
                    
            }

            template <typename PolRingType>
            static Product<  Power< PolRingType > >  extendedFactorPolynomial(const typename  PolRingType::Element & pol, const PolRingType & polRing) 
            {
                 Product<  Power< PolRingType >  >        prod;
                
                 std::vector<typename  PolRingType::Element > elements = factorPolynomial(pol,polRing);

                for (size_t pos=0;pos<elements.size(); pos++)
                {
                    int exponent=0;
                    std::pair<typename  PolRingType::Element,typename  PolRingType::Element> divisionResult;
                    divisionResult.first = pol;
                    do 
                    
                    {    
                         divisionResult = polRing.divide( divisionResult.first, elements[pos] );
                        if (not divisionResult.second.isZero() )
                            break;
                        exponent++;
                        
                    }
                    while ( divisionResult.second.isZero() );
                    assert( exponent>0 );
                
                    prod.push_back( Power< PolRingType >( elements[pos],exponent) );
                    
                }
                 return prod;
            }

          



            static void test()
            {
                std::cerr << " test factor polynomial with FLINT) " << std::endl;
                typedef TPolRingType::CoeffRingType CoeffRingType;
        
                int characteristic = 7;
                int epsPrecision = 0;

                CoeffRingType *    coeffRing  = new CoeffRingType(characteristic, epsPrecision );
                    TPolRingType     polRing (*coeffRing);

               TPolRingType::ElementType   pol(4);
                pol.setCoeff(0,1);
                pol.setCoeff(2,5);
                pol.setCoeff(4,1);
                std::cerr << "pol: " << pol << std::endl;
                std::vector< TPolRingType::ElementType> factors = FLINTFactorPolynomial::factorPolynomial( pol, polRing );

                std::vector< TPolRingType::ElementType>::iterator it;
        
                std::cerr << " factors " << std::endl;
                for (  it=factors.begin(); it!=factors.end(); it++)
                {
                    cout << *(it) << ",  " ; 
                }
                
                std::cerr  << std::endl;
                //  std::for_each( factors.begin(), factors.end(), []( TPolRingType::ElementType  n) { cout << n << " "; } );

                Shape shape = computeShape(pol,polRing);
                std::cerr << "shape " << shape;
                
            }

    };

   
}

