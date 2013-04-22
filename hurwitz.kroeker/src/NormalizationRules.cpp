


#include "NormalizationRules.h"
/*#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/foreach.hpp>
#include <boost/bind.hpp> 
#include <boost/lambda/if.hpp>
*/

namespace RationalMapSearch
{    

            const NormalizationValue NormalizationValue::infinity=NormalizationValue::createInfinityValue();
            const NormalizationValue NormalizationValue::zero=NormalizationValue::createZeroValue();
            const NormalizationValue NormalizationValue::one=NormalizationValue::createOneValue();

            const int NormalizationRule::dontcare = -1 ; 

            /// which polynomial to normalize  can be chosen arbitrarily if polynomialId==-1, same holds for exponent.
            NormalizationRule::NormalizationRule(  int polynomialId, 
                                int exponent,
                                NormalizationValue  value):     polynomialId_m(polynomialId),  
                                                                exponent_m(exponent),
                                                                value_m(value)
                                                                
            {
                assert(exponent_m>0 || exponent_m == NormalizationRule::dontcare );
                // only first three polynomials considered as candidates for normalization.
                assert( (polynomialId_m>=0 && polynomialId_m<3) ||  polynomialId_m == NormalizationRule::dontcare );
            }
    
            void NormalizationRule::checkNormalizationRule(const ShapeList& shapelist) const
            {
                if (polynomialId_m>=0)
                    assert((size_t)(polynomialId_m)< shapelist.size() );
                if (polynomialId_m>=0 && exponent_m>0)     //dann muss der Exponent im entsprechenden Polynom auch vorkommen?
                    assert(shapelist[polynomialId_m].hasExponent(exponent_m));            
            }        

            bool NormalizationRule::matches(const ShapeList& shapelist,int polynomialId ) const
            {
                    assert( polynomialId>=0  && (size_t)(polynomialId)< shapelist.size()  );               
                    if (polynomialId==polynomialId_m || polynomialId_m== NormalizationRule::dontcare )
                    {
                        if (exponent_m== NormalizationRule::dontcare )
                            return true;
                        else
                            return  shapelist[polynomialId].hasExponent(exponent_m);
                    }
                return false;
            }

            bool NormalizationRule::matches(const Shape& shape ) const
            {
                        if (exponent_m== NormalizationRule::dontcare )
                            return true;
                        else
                            return  shape.hasExponent(exponent_m);
            }

            bool NormalizationRule::matches(int polynomialId , int exponent ) const
            {   
                    assert( polynomialId>=0 && exponent>0);
                    if (polynomialId==polynomialId_m || polynomialId_m== NormalizationRule::dontcare )
                    {
                        if (exponent_m== NormalizationRule::dontcare || exponent==exponent_m)
                            return true;
                    }
                return false;
            }
        
        /*void    NormalizationRule::setPolynomialId(int polynomialId)  
        {
            assert( (polynomialId>=0 && polynomialId<3) ||  polynomialId == NormalizationRule::dontcare );

            polynomialId_m=polynomialId;
        }*/
 
        void    NormalizationRule::clearPolynomialId()
        {     
            polynomialId_m= NormalizationRule::dontcare ;
        }

        void    NormalizationRule::clearExponent()
        {     
            exponent_m= NormalizationRule::dontcare ;
        }


      int NormalizationRuleList::countValue(const std::vector<NormalizationRule>& nrl, NormalizationValue val)  
    {
        int count=0;
        std::vector<NormalizationRule>::const_iterator it;
        for (it= nrl.begin(); it!=nrl.end(); it++ )
            if ((*it).getValue()==val)
                count++;
        return count;
    }
   

    NormalizationRuleList::NormalizationRuleList(const std::vector<NormalizationRule> &normRuleList,bool strictNormalization ) : strictNormalization_m( strictNormalization )
    {

        assert (NormalizationValue::infinity!=NormalizationValue::one && 
                NormalizationValue::infinity!=NormalizationValue::zero &&
                NormalizationValue::zero!=NormalizationValue::one);

       // assert( 1 >= std::count_if ( normRuleList.begin(), normRuleList.end(),    boost::bind(&NormalizationRule::getValue(),_1)(boost::lambda::_1) == NormalizationValue::infinity  ));
       // assert( 1 >= std::count_if ( normRuleList.begin(), normRuleList.end(),    boost::bind(&NormalizationRule::getValue(),_1)(boost::lambda::_1) == NormalizationValue::zero  ));
       // assert( 1 >= std::count_if ( normRuleList.begin(), normRuleList.end(),    boost::bind(&NormalizationRule::getValue(),_1)(boost::lambda::_1) == NormalizationValue::one  ));

        assert( NormalizationRuleList::countValue(normRuleList,NormalizationValue::infinity )==1);
        assert( NormalizationRuleList::countValue(normRuleList,NormalizationValue::zero )==1 );
        assert( NormalizationRuleList::countValue(normRuleList,NormalizationValue::one )==1 );

        std::vector<NormalizationRule> normRuleListCopy= normRuleList;

        std::vector<NormalizationRule> normRuleListTmp;
        assert (normRuleList_m.size() <=3 );

        for ( size_t pos = 0;pos < normRuleListCopy.size() ;pos++)
        {
            if (normRuleListCopy[pos].getPolynomialId()>=0 &&normRuleListCopy[pos].getExponent() >0)
                normRuleList_m.push_back(normRuleListCopy[pos]);
            else
                normRuleListTmp.push_back(normRuleListCopy[pos]);
        }
        normRuleListCopy = normRuleListTmp;
        //normRuleListTmp.resize(0);
        normRuleListTmp.clear();
         
        for ( size_t pos = 0;pos < normRuleListCopy.size() ;pos++)
        {
            if (normRuleListCopy[pos].getPolynomialId()>=0 )
                normRuleList_m.push_back(normRuleListCopy[pos]);
            else
                normRuleListTmp.push_back(normRuleListCopy[pos]);
        }
        normRuleListCopy = normRuleListTmp;
        //normRuleListTmp.resize(0);
        normRuleListTmp.clear();
        
        for ( size_t pos = 0;pos < normRuleListCopy.size() ; pos++ )
        {
            if (normRuleListCopy[pos].getExponent() >=0 )
                normRuleList_m.push_back(normRuleListCopy[pos]);
            else
                normRuleListTmp.push_back(normRuleListCopy[pos]);
        }

            normRuleListCopy = normRuleListTmp;
         //normRuleListTmp.resize(0);
        normRuleListTmp.clear();
        
        for ( size_t pos = 0;pos < normRuleListCopy.size() ;pos++)
        {
                normRuleList_m.push_back( normRuleListCopy[pos] );
        }
        // check NormalizationValue
    }

    /// todo: es gibt zwei sorten von normalisierungsregeln:
    /// 1. bei der Polynomsuche und zweitens bei der Anwendung!
    NormalizationRuleList NormalizationRuleList::constructDefault(const ShapeList & shapelist)
    {
        int polynomialId = NormalizationRule::dontcare;  
        int exponent = NormalizationRule::dontcare;  
        //NormalizationValue value;
        //NormalizationRule   first( polynomialId  = NormalizationRule::dontcare, exponent = NormalizationRule::dontcare,   value = NormalizationValue::infinity);
        NormalizationRule   first( polynomialId  = 0, exponent = NormalizationRule::dontcare,    NormalizationValue::infinity);
        NormalizationRule   second( polynomialId = 1, exponent = NormalizationRule::dontcare,    NormalizationValue::zero);
        NormalizationRule   third( polynomialId  = NormalizationRule::dontcare, exponent = NormalizationRule::dontcare,   NormalizationValue::one);
        std::vector<NormalizationRule> normRuleList;
        normRuleList.push_back(first);
        normRuleList.push_back(second);
        normRuleList.push_back(third);
        assert(normRuleList.size()==3);
        return NormalizationRuleList(normRuleList);
    }

          
                
};
    
