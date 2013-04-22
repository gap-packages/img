#pragma once

 
#include <vector>
#include <algorithm>
#include <list>
#include <assert.h>
#include <cstdlib>
#include <map>
#include <ext/hash_map>
#include <hash_map>
#include <iostream>

extern "C" {
#undef __cplusplus
    # include "nmod_poly.h"
#define __cplusplus
}


#include "DebugLogger.h"


namespace RationalMapSearch
{
        inline void  printMpz(mpz_t num)       
            {                 
                   
                    char*  str= mpz_get_str( NULL, 10, num);
                
                    std::string counterStr  (str);
                    
                    std::cout << counterStr << ";" << std::endl;   
                    delete[] str;
            };
 
            inline void    getIrredCount(int degree, int cardinality, mpz_t result)   
            {
                //std::cerr << "getIrredCount (" << degree << ")" << std::endl;           
                mpz_set_ui(result,0);
                mpz_t tmp ;
                mpz_init(tmp);  

                for (int d =1; d<= degree; d++)
                {
                    if (degree % d ==0 )
                    {
                        mpz_set_ui(tmp,1);
                        for (int dexp=d; dexp>0; dexp--)
                        {
                            mpz_mul_ui( tmp, tmp, cardinality );                         
                        }
                        int degd= degree/d;                                        
                        int moebius = n_moebius_mu(degd);
                        //mpz_set_si( tmp2, moebius );
                        mpz_mul_si(tmp, tmp, moebius );

                        mpz_add(result,result,tmp);
                    }
                }
             
                mpz_div_ui(result,result,degree);         
                return  ;
            }

           
            inline size_t    getIrredCount(int degree, int cardinality) 
            {
                mpz_t irredListSizeMpz;
                mpz_init( irredListSizeMpz );
                getIrredCount( degree, cardinality, irredListSizeMpz );
                size_t result = mpz_get_ui ( irredListSizeMpz );
                mpz_t tmp;
                mpz_init( tmp );
                mpz_set_ui( tmp,result );
                
                if ( mpz_cmp(tmp,irredListSizeMpz)==0 )
                    return result;
                else
                {
                    std::cerr << " too many irreducibles... ";
                    assert(false);
                }
            }

    /// @implementation only for the rational map search needs. Doesn't have to be efficient, because only factorisations of very small integers are required.
    /// alternatively use FLINT or http://www.boo.net/~jasonp/qs.html or interface to an CAS.
    class IntFactorTable
    {

        public:
            typedef unsigned int IntType ;
            typedef std::vector<IntFactorTable::IntType> IntFactorVecType ;        
        
            typedef std::map< IntType, IntFactorVecType >  IntFactorTableType;

            static const IntType MaxInt  ;
            
             bool *  isPrime_m  ;

            IntFactorTable()  ;

        private:
            IntFactorTableType  intFactorTable_m;

            
            IntFactorVecType computeIntFactors(IntType integer);
          
        
            void addTableEntry(IntType integer)
            {
                  intFactorTable_m.insert(std::pair< IntFactorTable::IntType, IntFactorVecType > ( integer, computeIntFactors(integer) ) );
            }
            
        public:

            // todo: templatize return type?
            IntFactorVecType getFactors(IntType integer);

    
            static void test();
            

    };

  
    //http://www.parashift.com/c++-faq-lite/pointers-to-members.html
    // Zitat: use a Typedef for member functions !!!

  ///deprecated
  template <typename IrreduciblePolTableType>
    class PseudoIrreducibleExpensiveTest
    {
        public:

            static const bool threadsafe;

            //typedef typename IrreduciblePolTableType::PolRingType    TPolRing;

            //typedef typename TPolRing::ElementType  VecElemType;
          
              typedef bool (* IrredTesterFktType)(typename IrreduciblePolTableType::PolRingType::Element & , const typename IrreduciblePolTableType::PolRingType & ,   IrreduciblePolTableType & );

          
            inline static bool isIrreducible(typename IrreduciblePolTableType::PolRingType::Element & pol, const typename IrreduciblePolTableType::PolRingType & polRing,   IrreduciblePolTableType & irredTable)  
            {
                    typedef typename IrreduciblePolTableType::PolRingType    TPolRing;
                    typedef typename  IrreduciblePolTableType::PolRingType::ElementType  VecElemType;

                    int degree = pol.getExactDegree();
                    
                    // check if is Zero  pol(a= for a in 0...char coeffring
                    const typename TPolRing::CoeffRingType & coeffRing = polRing.getCoeffRing();
    
                    for (int coeff = 0; coeff< coeffRing.getCharacteristic(); coeff++ )
                    {
                        if (  polRing.evalAt(pol, coeffRing.Convert(coeff)  )== TPolRing::RingType::ElementType::Zero) 
                            return false;
                    }               
                    for (int currDegree= 2; currDegree<= degree /2 ;currDegree++)
                    {
                           typename IrreduciblePolTableType::IrredTableType & irredPolList = irredTable.getIrredPolList(currDegree);

                           bool finished = false;
                           #pragma omp parallel for shared(finished)
                           for (int pos=0; pos< irredPolList.size();pos++)
                           {
                                if (!finished)
                                {
                                
                                    typename TPolRing::Element gcdf = polRing.gcd( pol, *(irredPolList[pos])  );
                                    if ( not gcdf.isConstant() )
                                    {
                                            finished=true;
                                    }
                                }
                               
                            }
                            if (finished)
                                    return false;
                    
                    }
                    return true;
            };

    };

    template <typename IrreduciblePolTableType>
    const bool PseudoIrreducibleExpensiveTest<IrreduciblePolTableType>::threadsafe = true;

    template <typename IrreduciblePolTableType>
    class FLINTIrreducibleTest
    {

        public:

            static const bool threadsafe ;

            //typedef typename IrreduciblePolTableType::PolRingType    TPolRing;
            //typedef typename TPolRing::ElementType  VecElemType;

              typedef bool (* IrredTesterFktType)(typename IrreduciblePolTableType::PolRingType::Element & , const typename IrreduciblePolTableType::PolRingType & ,   IrreduciblePolTableType & );
          
            inline static bool isIrreducible(typename IrreduciblePolTableType::PolRingType::Element & pol, const typename IrreduciblePolTableType::PolRingType & polRing,   IrreduciblePolTableType & irredTable)  
            {
                    //std::cerr << "FLINT::isIrreducible" << std::endl;
                    typedef typename IrreduciblePolTableType::PolRingType    TPolRing;
                    typedef typename  IrreduciblePolTableType::PolRingType::ElementType  VecElemType;

                    int degree = pol.getExactDegree();

                      // pol.getCoeffRing()
                    // check if is Zero  pol(a= for a in 0...char coeffring
                    const typename TPolRing::CoeffRingType & coeffRing = polRing.getCoeffRing();

                    nmod_poly_t x ;
                    nmod_poly_init(x ,coeffRing.getCharacteristic() );
    
                    //typename TPolRing::Element::CoefficientType coeff = TPolRing::Element::CoefficientType::Zero;
    
                    for (int coeff = 0; coeff< coeffRing.getCharacteristic(); coeff++ )
                    {
                        if (  polRing.evalAt(pol, coeffRing.Convert(coeff)  )== TPolRing::RingType::ElementType::Zero) 
                            return false;
                    }               
    
                    // riscy: what if pol.getCoeff is not nonnegative?
                    for (int currDegree= 0; currDegree<=degree; currDegree++)
                    {
                         assert( pol.getCoeff( currDegree).getX()>=0 );
                         nmod_poly_set_coeff_ui(x,   currDegree,  pol.getCoeff( currDegree).getX() );
                    }
                    bool result = nmod_poly_is_irreducible(x);
                    nmod_poly_clear (x);

                       return result ;

                    
            };

    };

   template <typename IrreduciblePolTableType>
    const bool FLINTIrreducibleTest<IrreduciblePolTableType>::threadsafe = true;



   ///deprecated, checks only for linear and quadratic factors.
   template <typename IrreduciblePolTableType>
    class PseudoIrreducibleTest
    {
        public:

            static const bool threadsafe;

            typedef typename IrreduciblePolTableType::PolRingType    TPolRing;

            typedef bool (* IrredTesterFktType)(typename IrreduciblePolTableType::PolRingType::Element & , const  typename IrreduciblePolTableType::PolRingType &,  IrreduciblePolTableType &);

            
            inline static bool isIrreducible(typename IrreduciblePolTableType::PolRingType::Element & pol, const typename IrreduciblePolTableType::PolRingType & polRing,  IrreduciblePolTableType & irredTable)  
            {
                    // check if is Zero  pol(a= for a in 0...char coeffring
                    const typename  IrreduciblePolTableType::PolRingType::CoeffRingType & coeffRing = polRing.getCoeffRing();
        
                    for (int coeff = 0; coeff< coeffRing.getCharacteristic(); coeff++ )
                    {
                        if (  polRing.evalAt(pol, coeffRing.Convert(coeff)  )== TPolRing::RingType::ElementType::Zero) 
                            return false;
                    }               
                    return true;
            }
    };
    template <typename IrreduciblePolTableType>
    const bool PseudoIrreducibleTest<IrreduciblePolTableType>::threadsafe = true;

    // todo: factory class which creates IrreduciblePolTable. Eventually it is sufficient to derive classes instead template parametrizing.

    // todo: in construct polynomial Set : somehow get rid of recursion if possible,
    template <class TPolRing>
    class IrreduciblePolTable
    {

        public:

            typedef TPolRing    PolRingType;
            typedef bool (* IrredTesterFktType)(typename IrreduciblePolTable<TPolRing>::PolRingType::Element & , const TPolRing & , IrreduciblePolTable<TPolRing> & );

            typedef std::vector<const typename TPolRing::Element * > IrredVecListType;

            typedef std::map< unsigned int, IrredVecListType >  IrredTableType;

        private: 
   
        
         TPolRing polRing_m;

        typename TPolRing::CoeffRingType  coeffRing_m;
        // TPolIrredTester irredTester_m;

        //bool (* irredTester_m)(typename TPolRing::Element & , const TPolRing &  );

            IrredTesterFktType  irredTester_m;
 
            IrredTableType      irredPolTable_m;

        bool parallelize_m;


        public:
                IrreduciblePolTable( const TPolRing & polRing, 
                                       IrredTesterFktType irredTester ) :  polRing_m(polRing),
                                                                        coeffRing_m( polRing.getCoeffRing() ),
                                                                      irredTester_m(irredTester),
                                                                        parallelize_m(false)
                {};
                
                IrreduciblePolTable( const TPolRing & polRing, 
                                       IrredTesterFktType irredTester,
                                        bool parallelize ) :  polRing_m(polRing),
                                                                        coeffRing_m( polRing.getCoeffRing() ),
                                                                      irredTester_m(irredTester),
                                                                        parallelize_m(parallelize)
                {};

                IrreduciblePolTable( const TPolRing & polRing, 
                                       IrredTesterFktType irredTester,
                                        unsigned int maxDegree ) :  polRing_m(polRing),
                                                                        coeffRing_m( polRing.getCoeffRing() ),
                                                                        irredTester_m(irredTester),
                                                                        parallelize_m(false)
                {
                    for (int deg = 1; deg<=maxDegree;deg++)
                    {   
		        size_t irredCount = getIrredCount( deg, polRing.getCardinality() );
			DebugLogger::logStream() << "irredCount to compute for degree " << deg << " is " << irredCount << std::endl;
                      
                        computeIrredPolList(deg);
                    }
                };


                IrreduciblePolTable( const TPolRing & polRing, 
                                       IrredTesterFktType irredTester,
                                        unsigned int maxDegree,
                                        bool parallelize ) :  polRing_m(polRing),
                                                                        coeffRing_m( polRing.getCoeffRing() ),
                                                                        irredTester_m(irredTester),
                                                                        parallelize_m(parallelize)
                {
                    for (int deg = 1; deg<=maxDegree;deg++)
                    {   
                        // just to make sure there are not too many irreducibles.
                        size_t irredCount = getIrredCount( deg, polRing.getCardinality() );
			DebugLogger::logStream() << "irredCount to compute for degree " << deg << " is " << irredCount << std::endl;
                        computeIrredPolList(deg);
                    }
                };
                
                void computeIrredPolList(unsigned int degree)  
                {
                       
                       #ifdef SAFE
                        assert(degree>0);
                       #endif
                      if ( irredPolTable_m.find(degree ) == irredPolTable_m.end() )
                      {          

                            
                            std::list< const typename TPolRing::Element *>    irredList;
			    
			    // TODO: it is possible to precompute the size of irredList, but then 
			    // how it should be parallelized?
                            if (degree==1) 
                            {
                                typename  TPolRing::Element irredPol(typename  TPolRing::Element(1) ) ;
                                for (int root = 0; root<  coeffRing_m.getCharacteristic(); root++ )
                                {                                
                                    irredPol.setCoeff( 0, coeffRing_m.addInv( coeffRing_m.Convert(root) ) );
                                    irredPol.setCoeff( 1, TPolRing::CoeffRingType::ElementType::One );
                                    irredList.push_back( new  const typename TPolRing::Element ( irredPol)  ); 
                                }           
                            }
                            else
                            {
                                // loop through all polynomials and check for Irreducibility.
                                // 
                                // two possibilities: either construct all polynomials first
                                // or use nextPolynomial (will be slower.)
                                // should decide which one to use by looking at the total count and memory usage
                                typename TPolRing::Element   pol = typename TPolRing::Element(degree) ;
                                // pol.SetCoeff( 0,TPolRing::Element::ElementType::One );
                              
                                
                                pol.setCoeff( degree, TPolRing::CoeffRingType::ElementType::One );


                                //parallelize_m=true;

                                if (parallelize_m)
                                {
                                
                                    size_t tmpVecSize=64;
                                    std::vector< typename TPolRing::Element >    tmpPolVec(tmpVecSize,pol);
                                    bool done = false;
                                    while ( !done )
                                    {
                                        for (size_t i=0; i<tmpVecSize; i++)
                                        {
                                            tmpPolVec[i]=pol;
                                            done = ! pol.nextInPlace( coeffRing_m, true);
                                            if (done)
                                            {
                                                tmpPolVec.resize(i+1);
                                                break;
                                            }
                                        }
                                        #pragma omp parallel for shared(irredList)
                                        for (size_t j = 0 ; j< tmpPolVec.size();j++) 
                                        {
                                            if ( irredTester_m( tmpPolVec[j], polRing_m, *this ) ) 
                                            {
                                                #pragma omp critical
                                                irredList.push_back( new  typename TPolRing::Element ( tmpPolVec[j] ) ); 
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    do  
                                    {
        
                                            //std::cerr << "pol " << pol << std::endl;
                                            //getchar();
                                            // speed improvement: estimate number of irreducibles and use a vector instead of a list and then shrink the vector.
                                            if ( irredTester_m( pol, polRing_m, *this ) ||  degree==1) 
                                                irredList.push_back(  new const typename TPolRing::Element ( pol) ); 
                                    }
                                    while ( pol.nextInPlace( coeffRing_m, true) );        
                                }                 
                            }
                            // now, copy the list to a vector. 
                            //#pragma omp critical irredTable
                            #ifdef VERBOSE
                            std::cerr << "computed " << degree << std::endl;
                            #endif
			    DebugLogger::logStream() << "computed " << degree << std::endl;
			    
                            irredPolTable_m[degree]  = IrredVecListType( irredList.begin(), irredList.end() );  
                            
                      }
                }

                void updateIrredPolList(unsigned int degree)  
                {
                       #ifdef SAFE
                        assert(degree>0);
                       #endif
 
                        for (unsigned int deg=1; deg<=degree;deg++)
                        {
                            if ( irredPolTable_m.find( deg ) == irredPolTable_m.end() )
                            {       
                                computeIrredPolList(deg);
                                #ifdef VERBOSE
                                std::cerr << "irredPolTable_m[degree].size()" << irredPolTable_m[deg].size() <<  std::endl;
                                #endif
				DebugLogger::logStream() <<"irredPolTable_m[degree].size()" << irredPolTable_m[deg].size() <<  std::endl;
                            }
                        }
                        
                }

             
                const IrredVecListType * getIrredPolListPtr(unsigned int degree)  
                {
                       #ifdef SAFE
                        assert(degree>0);
                       #endif
                     typename IrredTableType::const_iterator it=irredPolTable_m.find(degree);  
                     assert( it!= irredPolTable_m.end() );
                     return  ( &((*it).second) );
                }
                
                const IrredVecListType & getIrredPolList(unsigned int degree)  const
                {
                    //std::cerr << "getIrredPolList (" << degree << ")" << std::endl;
                    typename IrredTableType::const_iterator it=irredPolTable_m.find(degree);  
                    #ifdef VERBOSE
                    if (it == irredPolTable_m.end() )                     
                        std::cerr << "getIrredPolList (" << degree << ")" << std::endl;
                    #endif
                    assert( it != irredPolTable_m.end() );
                      return  ( (*it).second );
                }
            
                static void test()
                {
                    int characteristic = 11;
                    int epsPrecision = 0;

                    ///@todo hier gibt es dublicate deletes, wenn coeffRing nicht als Pointer definiert wird !!!
                    typename TPolRing::CoeffRingType *    coeffRing  = new typename TPolRing::CoeffRingType(characteristic, epsPrecision );
                    TPolRing     polRing (*coeffRing);
                    
                   // IrreduciblePolTable poltable( polRing, PseudoIrreducibleTest<TPolRing>::isIrreducible ) ;
 
                    typedef IrreduciblePolTable<TPolRing> IrreduciblePolTableType;

                   // typedef PseudoIrreducibleTest< IrreduciblePolTableType  > IrreducibleTestType;
                   //typedef PseudoIrreducibleExpensiveTest< IrreduciblePolTableType  > IrreducibleTestType;

 
                    typedef FLINTIrreducibleTest< IrreduciblePolTableType  > IrreducibleTestType;
                    //typedef GAPIrreducibleTest< IrreduciblePolTableType  > IrreducibleTestType;


                   //  PseudoIrreducibleExpensiveTest< IrreduciblePolTable<TPolRing> >::IrredTesterFktType fkt= PseudoIrreducibleExpensiveTest< IrreduciblePolTable<TPolRing> >::isIrreducible< IrreduciblePolTable<TPolRing> >;
                 

                    IrreduciblePolTable<TPolRing> poltable( polRing, IrreducibleTestType::isIrreducible ) ;

                    int maxDegree = 7;
                    poltable.updateIrredPolList( maxDegree );

                    for (int degree=1; degree < (maxDegree/ 2) ; degree++)
                         poltable.getIrredPolList(degree);
                   

                
                    assert( poltable.getIrredPolList(1).size()==11);
                    assert( poltable.getIrredPolList(2).size()==55);
                    assert( poltable.getIrredPolList(3).size()==440);
                    
                    
                }
        
            virtual ~IrreduciblePolTable() 
            {
                // todo: free memory.
            };

    };


}