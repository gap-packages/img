#pragma once

#include "Shape.h"
#include <iostream>

namespace RationalMapSearch
{
    /*enum class NormalizationValue
        {   
            infinity=-1,
            zero = 0,
            one = 1
        };
    */

    /// is an enum emulation; not completely error-prone but acceptable to provide backward compatibility instead of enums.
   class NormalizationValue
        {   
           

            private:
                int enumVal_m;

           
             NormalizationValue(int val)
            {
                enumVal_m=val;
            };

            protected:
                static NormalizationValue createInfinityValue()
                {
                    return NormalizationValue(-1);
                }
                static   NormalizationValue createZeroValue()
                {
                    return NormalizationValue(0);
                }
                static NormalizationValue createOneValue()
                {
                    return NormalizationValue(1);
                }

            public:
                bool operator==(const NormalizationValue& value) const
                {   
                    return value.enumVal_m==enumVal_m;
                }
              bool operator!=(const NormalizationValue& value) const
                {   
                    return value.enumVal_m!=enumVal_m;
                }

                static const NormalizationValue infinity;
                static const NormalizationValue zero;
                static const NormalizationValue one;
        };



    /// potenzielle Probleme: Konvertierung von und zu Shape::ScalarType
    class NormalizationRule
    {
        private:

              int                   polynomialId_m ;
              int     exponent_m;
              NormalizationValue    value_m;
            
        
            
        public:
    
  
            static const int dontcare ;

            NormalizationValue getValue() const { return value_m; };

            int     getPolynomialId()   const { return polynomialId_m; };
            int     getExponent()       const    { return exponent_m; };

            

            /// which polynomial to normalize  can be chosen arbitrarily if polynomialId==-1, same holds for exponent.
            NormalizationRule(  int polynomialId, 
                                int exponent,
                                NormalizationValue  value);

        
    
            void checkNormalizationRule(const ShapeList& shapelist) const;

            bool matches(const ShapeList& shapelist,int polynomialId ) const;

            bool matches(const Shape& shape ) const;

            bool matches(int polynomialId , int exponent ) const;


            //void setPolynomialId(int polynomialId)  ;

            void clearPolynomialId();
            void clearExponent();

        // to use std::vector:  
            //This means I must define a operator= with unclear semantics just because it is never used. Oh my noodles! 
            // see for reference http://blog.copton.net/archives/2007/10/13/stdvector/index.html
          /*  NormalizationRule& operator=(const NormalizationRule& rhs) { 
                std::cerr<< "abort in NormalizationRule" << std::endl;
                abort();
                return *this;
            }*/

    };

   


    class NormalizationRuleList
    {
        private:
            std::vector<NormalizationRule> normRuleList_m;
            bool strictNormalization_m;

        public:

            static int countValue(const std::vector<NormalizationRule>& nrl, NormalizationValue val) ;
            ///@optional jede mögliche Normalisierung probieren, so dass möglichst viele angewandt werden können.

            NormalizationRuleList(const std::vector<NormalizationRule> &normRuleList, bool strictNormalization=false ) ;

            //containsMatchingRules

            static NormalizationRuleList constructDefault(const ShapeList & shapelist);

            std::vector<NormalizationRule>  getNormalizationRuleListAsVector()  const {   return normRuleList_m;  }

                bool strictNormalization()     const     {  return strictNormalization_m;   };

            NormalizationRule operator [](size_t pos)
            {
                //assert( pos<normRuleList_m.size() );
                return normRuleList_m.at(pos);
            }

            NormalizationRule getRuleByPolynomialId(int polynomialId)
            {
                for (size_t pos=0; pos< normRuleList_m.size(); pos++)
                {   
                    if (normRuleList_m[pos].getPolynomialId()==polynomialId)            
                        return normRuleList_m[pos];
                }
                assert (false);
            }

            NormalizationRule getRuleByNormValue(NormalizationValue val)
            {   
                for (size_t pos=0; pos< normRuleList_m.size(); pos++)
                {   
                    if (normRuleList_m[pos].getValue()==val)            
                        return normRuleList_m[pos];
                }
                assert (false);
            }

            size_t size() const { return normRuleList_m.size(); }
    };

};
    