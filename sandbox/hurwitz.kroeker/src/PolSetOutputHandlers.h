
#pragma once

  
#include "hmfTypedefs.h"
#include "OutputMode.h"

#include <iostream>
#include <sstream>

namespace RationalMapSearch
{


    template<class TPolRingTypePar>
    class FiniteFieldSearch;


    
    // preconditions:   
    //-   TPolRingTypePar::ElementType
    //-   TPolRingTypePar::CoeffRingType
    //-   TPolRingTypePar::CoeffRingType::operator[ exponent ]()
    //-   TPolRingTypePar::CoeffRingType::ElementType
    //-   TPolRingTypePar::CoeffRingType::getExactDegree()
    //-   TPolRingTypePar::getCoeffRing()
    //-   TPolRingTypePar::CoeffRingType::elemToGeneratorExponent()

    template<class TPolSet>
    class IOutputHandler
    {
        public:
            virtual void startOutput( )=0;
            virtual void finishOutput( )=0;
            virtual void print(const TPolSet& polSet)=0;

            virtual ~IOutputHandler() {};
    };

    template<class TPolSet>
    class EmptyOutputHandler  : public IOutputHandler<TPolSet>
    { 
            public:
        
            EmptyOutputHandler( ) {};
            virtual void startOutput( ) {};
            virtual void finishOutput( ) {};
            virtual void print(const TPolSet& polSet) 
            {};
            virtual ~EmptyOutputHandler() {};
    };

    template<class TPolSet>
    class GAPOutputHandler  : public IOutputHandler<TPolSet>
    { 
        private:
            bool printFirst_m;
        public:

 

            GAPOutputHandler():printFirst_m(true)
            {}
        
      
            std::ostream & start( std::ostream & os)
            {   
                  os << " [ " ;
                  return os;
            }

            std::ostream & finish( std::ostream & os)
            {   
                  os << " ] ;" << std::endl;
                  return os;
            }

            virtual void startOutput( )
            {   
                start(std::cout);
            }
            virtual void finishOutput( )
            {   
                finish(std::cout);
            }


            std::ostream &  printPolynomial(std::ostream &os, const typename TPolSet::PolynomialRingType &  polRing, 
                                                              const typename TPolSet::ElementType   &  pol, 
                                                              int maxDegree)
            {

                typedef typename TPolSet::PolynomialRingType::CoeffRingType::ElementType    CoeffType;
                  int degree = pol.getExactDegree();
                 assert(degree>=0);
                 bool first=true;
                os << " [ " ;
                int currDegree=0;
                 for (  currDegree=0 ; currDegree<=degree; currDegree++)
                {
                    if (! first)
                                        os << ", " ;
                    first = false;
                    {
                        CoeffType coeff = pol[ currDegree ] ;
                         if (coeff.isZero() )
                              os <<   " 0*Z(" <<  polRing.getCoeffRing().getCardinality() << ") " ;
                        else
                        {
                            std::stringstream oss;
                        
                            int exponent =    (polRing.getCoeffRing() ).elemToGeneratorExponent(coeff);
                             oss << exponent    ;
                        
                              os <<  "  Z(" <<  (polRing.getCoeffRing().getCardinality()) << ")^" << oss.str() << " ";   
                        }
                    }
                }
                for (    ; currDegree<=maxDegree; currDegree++)
                { 
                    if (! first)
                        os << ", " ;
                    first = false;
                    os <<   " 0*Z(" <<  polRing.getCoeffRing().getCardinality() << ") " ;
                }
                os << " ] " << std::endl;
                return os;
            }
            
    protected:
           std::ostream & iprint( std::ostream &os, const TPolSet & polSet )
           {
               
                if (! printFirst_m)
                {
                    os << " , " ;
                }
                os << " [ " ;
            
                bool firstPolynomial=true;

                TPolSet  polSetCopy = polSet;

                typename TPolSet::ElementType pol = polSet[0];
                
 
                for (size_t pos=0; pos< polSetCopy.size(); pos++)
                {
                    if (! firstPolynomial)    os << " , " ;

                    firstPolynomial=false;
                    // print polynomial in some style, dependent on the polynomial.
                    printPolynomial(os, polSet.getRing() , polSetCopy[pos], polSet.getDegree() );
                }
                os << " ] " ;

                printFirst_m = false;               
                return os;
           }

        public :
           virtual void print( const   TPolSet&   polSet  )
           {
                {
                  iprint(std::cout,  polSet );
                }
           }
            virtual ~GAPOutputHandler() {};
    };
  template<class TPolSet>
    class M2OutputHandler  : public IOutputHandler<TPolSet>
    { 
        private:
            bool printFirst_m;
        public:

 

            M2OutputHandler():printFirst_m(true)
            {}
        
      
            std::ostream & start( std::ostream & os)
            {   
                  os << " { " ;
                  return os;
            }

            std::ostream & finish( std::ostream & os)
            {   
                  os << " } ;" << std::endl;
                  return os;
            }

            virtual void startOutput( )
            {   
                start(std::cout);
            }
            virtual void finishOutput( )
            {   
                finish(std::cout);
            }

        
            std::ostream &  printPolynomial(std::ostream &os, const typename TPolSet::PolynomialRingType &  polRing, 
                                                              const typename TPolSet::ElementType   &  pol, 
                                                              int maxDegree)
            {

                typedef typename TPolSet::PolynomialRingType::CoeffRingType::ElementType    CoeffType;
                int i;
            
                CoeffType    z;
            
                bool first = true;
            
                for (i=0; i<= pol.getExactDegree() ; i++)
                {
                    
                    z = pol.getCoeff(i); 
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
                return os;
            }
            
    protected:
           std::ostream & iprint( std::ostream &os, const TPolSet & polSet )
           {
               
                if (! printFirst_m)
                {
                    os << " , " ;
                }
                os << " { " ;
            
                bool firstPolynomial=true;

                TPolSet  polSetCopy = polSet;

                typename TPolSet::ElementType pol = polSet[0];
                
 
                for (size_t pos=0; pos< polSetCopy.size(); pos++)
                {
                    if (! firstPolynomial)    os << " , " ;

                    firstPolynomial=false;
                    // print polynomial in some style, dependent on the polynomial.
                    printPolynomial(os, polSet.getRing() , polSetCopy[pos], polSet.getDegree() );
                }
                os << " } " ;

                printFirst_m = false;               
                return os;
           }

        public :
           virtual void print( const   TPolSet&   polSet  )
           {
                {
                  iprint(std::cout,  polSet );
                }
           }
            virtual ~M2OutputHandler() {};
    };


}
