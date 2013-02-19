#pragma once

#include "PolynomialShape.h"
#include "FactorPolynomialWrapper.h"
#include "DebugLogger.h"

namespace RationalMapSearch 
{

    template <class TPolRingTypePar>
    std::list< PolynomialFactorBluePrint>     FiniteFieldSearch<TPolRingTypePar>::createPolFactorBlueprintList(  )
            {   

                 const ShapeList & shapelist = hurwitzMapSearchProblem_m.getConstShapeListRef();

                 std::list< PolynomialFactorBluePrint> polfactorBlueprintList;

                 for (short shapepos=0; shapepos<=1; shapepos++ )
                 {
                    Shape::MultiplicityDegreeHashType  muldegrep = shapelist[shapepos].getMultiplicityDegreeRep();
                    Shape::MultiplicityDegreeHashType::const_iterator it;
                    for ( it = muldegrep.begin(); it != muldegrep.end(); it++ )
                    {   
                            //std::cerr <<  (*it).first<< "^" <<(*it).second << std::endl;
                            const PolynomialFactorBluePrint bp( (*it).first, (*it).second, shapepos);
                            //std::cerr << bp<< std::endl;
                            polfactorBlueprintList.push_back(bp);
                    }
                }   
                return polfactorBlueprintList;
            }

             
        template <class TPolRingTypePar>
            FiniteFieldSearch<TPolRingTypePar>::FiniteFieldSearch(const HMSProblem & hms_, const SearchOptions& so, 
                                                                                     const TPolRingTypePar &  polynomialRing): hurwitzMapSearchProblem_m( hms_),
                                                                                                            searchOptions_m( so ),
                                                                                                            outputHandler_m(NULL),
                                                                                                            characteristic_m( polynomialRing.getCoeffRing().getCharacteristic() ),
                                                                                                          
                                                                                                            //polynomialRing_m( getPolRingRef(characteristic ) ),
                                                                                                            polynomialRing_m(polynomialRing),
                                                                                                            field_m( polynomialRing.getCoeffRing() )
                                                                                                        
                                                                                                         
            {
                mpz_init(counter_m);
                mpz_init(counterMod_m);
                mpz_set_ui(counter_m,0);
                mpz_set_ui(counterMod_m,1);
   
                setPolRingRef( & polynomialRing );

                assert(hms_.getNormalizationRuleList().getNormalizationRuleListAsVector().size()>0);
                assert(hurwitzMapSearchProblem_m.getNormalizationRuleList().getNormalizationRuleListAsVector().size()>0);

                IrreduciblePolTableType * irredPolTablePtr;
               
                irredPolTablePtr = getIrredPolTablePtr( field_m.getCharacteristic() );           


		DebugLogger::logStream() << "# constructionMaxFactorDegree " << hurwitzMapSearchProblem_m.getShapeList().getConstructionMaxFactorDegree() << std::endl;

		// todo: print for HMSProblem.
                //std::cerr << "constructionMaxFactorDegree " << hurwitzMapSearchProblem_m.getShapeList().getConstructionMaxFactorDegree() << std::endl;

                if ( not   searchOptions_m.dryRun() )
                    irredPolTablePtr->updateIrredPolList ( hurwitzMapSearchProblem_m.getShapeList().getConstructionMaxFactorDegree() );
                else
                    irredPolTablePtr->updateIrredPolList ( 1 );


                for (size_t pos=0; pos< hurwitzMapSearchProblem_m.getMinimalPolynomials().size() ;pos++)
                {
                   typename TPolRingTypePar::ElementType pol = convertPolRepToRingElem( hurwitzMapSearchProblem_m.getMinimalPolynomials()[pos] );
                    minimalPolynomials_m.push_back(pol);
                }
                
                
                
                if (searchOptions_m.outputMode()==OutputMode::GAPOutput)  //requires nested enumerators
                     outputHandler_m = new GAPOutputHandler<PolynomialSet<TPolRingTypePar> >();
                else if (searchOptions_m.outputMode()==OutputMode::M2Output)  //requires nested enumerators
                      outputHandler_m = new M2OutputHandler<PolynomialSet<TPolRingTypePar> >();
                //else
                //    outputHandler_m= new EmptyOutputHandler<PolynomialSet<TPolRingTypePar> >();
                assert( outputHandler_m != NULL);
		DebugLogger::logStream() << " FiniteFieldSearch initialized!" << std::endl;
              }

            /// convention: polRep is a vector with polynomial coefficients and polrep[i] is the coefficient of monomial x^i .
            template <class TPolRingTypePar>
            typename TPolRingTypePar::ElementType   FiniteFieldSearch<TPolRingTypePar>::convertPolRepToRingElem(const typename  HMSProblem::PolynomRepType & polRep) const
            {
                    typename TPolRingTypePar::ElementType pol( polRep.size() -1 );
                    for (size_t pos=0; pos< polRep.size(); pos ++)
                    {
                        pol.setCoeff( pos, field_m.Convert( polRep[pos] ) );
                    }
                    return pol;
            }
 
            template <class TPolRingTypePar>
            std::list< PolynomialFactorBluePrint>   FiniteFieldSearch<TPolRingTypePar>::createPolFactorConstructionRules(std::vector<int > partition, uint exponent, uint destpolynomial)
            {
                std::list< PolynomialFactorBluePrint> polFactorConstructionRuleList;
        
                for (size_t pos = 0; pos<partition.size() ; pos++)
                {
                   polFactorConstructionRuleList.push_back(PolynomialFactorBluePrint(exponent, partition[pos], destpolynomial) );
                }
                return polFactorConstructionRuleList;
            }

            // extract a PolynomialFactorBluePrint to which is is possible to apply the Normalization rule.
            template <class TPolRingTypePar>
            PolynomialFactorBluePrint*  FiniteFieldSearch<TPolRingTypePar>::extractMatchingRule(std::list< PolynomialFactorBluePrint> &polFactorConstructRules, NormalizationRule rule)
            {
                    //std::cerr << "extractMatchingRule" << std::endl;

                    PolynomialFactorBluePrint * res=NULL;

                    polFactorConstructRules.sort( PolynomialFactorBluePrint::minDegreeMaxMultiplicityLower );
                    std::list< PolynomialFactorBluePrint>::iterator it;

                    // wenn man erase benutzt, kann man  keine for-schleife verwenden.
                    it = polFactorConstructRules.begin();
                    while ( it != polFactorConstructRules.end()  )
                    {
                        /*if ( rule.matches( (*it).polynomialId_m, (*it).multiplicity_m ) )
                        {
                            std::cerr << "(*it).polynomialId_m" << (*it).polynomialId_m << std::endl;
                            std::cerr << "(*it).multiplicity_m" << (*it).multiplicity_m << std::endl;
                            std::cerr << "matches";
                        }*/
                           
                        if ( rule.matches( (*it).polynomialId_m, (*it).multiplicity_m ) &&   (*it).degree_m==1 )
                        
                        {   
                            res = new PolynomialFactorBluePrint( (*it) );                     
                            //std::cerr << "polFactorConstructRules size" << polFactorConstructRules.size() << std::endl;
                            it = polFactorConstructRules.erase(it);
                            //std::cerr << "polFactorConstructRules size" << polFactorConstructRules.size()<< std::endl;
                            return res;
                        }
                        else
                            it++;
                    }                 
                  
                    return res;
            }

            template <class TPolRingTypePar>
            bool     FiniteFieldSearch<TPolRingTypePar>::computeScalingFactors(const PolSetType &  polSet,   std::vector<  typename TPolRingTypePar::CoeffRingType::ElementType > & scalingFactors)
            {
                
                assert( polSet.size()>2 );
                scalingFactors.resize( polSet.size()-2 );

                bool computed=false;

                for ( size_t pos = 2; pos<polSet.size(); pos++ )
                {
                    computed=false;
                    for (int coeff = 1; coeff< field_m.getCharacteristic(); coeff++)
                    {
                        typename TPolRingTypePar::ElementType pol= polynomialRing_m.add(polSet[1], polynomialRing_m.scalarMultiply( field_m.Convert(coeff), polSet[0] ));
                        if (pol==polSet[pos])
                        {
                            scalingFactors[pos-2] = field_m.addInv(field_m.Convert(coeff));
                            computed=true;
                            break;
                        }
                    }
                    if (! computed)
                     break;
                }
                return computed;
            }


            /*
              for alphaFactorPos in 1..#alphaFactors-1 do
            (
                currScalingFactor:= alphaFactors#alphaFactorPos;
                if (normalizedByFirstScalar) then
                    currScalingFactor=currScalingFactor*(alphaFactors#0)^-1;
                if (sub (minPolyList#(alphaFactorPos-1),matrix {{ currScalingFactor}} )!=0 ) then
                (
            */

            

            /// @todo: polSetBlueprint wird nur gebraucht, wenn man auch eine Liste der Faktoren mitführen will
            template <class TPolRingTypePar>
            bool    FiniteFieldSearch<TPolRingTypePar>::processNormalizationRules( PolSetBlueprintType & polsetblueprint, 
                                                                                    std::list< PolynomialFactorBluePrint> &polFactorConstructRules,
                                                                                    std::vector<int> & degreeOneRootList,
                                                                                    ShapeList & shapeList,
                                                                                    FiniteFieldSearch<TPolRingTypePar>::PolSetType & polSet )
            {
                NormalizationRuleList   nrl = hurwitzMapSearchProblem_m.getNormalizationRuleList();
               
                std::vector<NormalizationRule>  nrlvec = nrl.getNormalizationRuleListAsVector();

                assert(nrlvec.size()>0);

                IrreduciblePolTableType & irredPolTableRef = getIrredPolTableRef( field_m.getCharacteristic() );
                

                std::vector<NormalizationRule>::iterator it =nrlvec.begin();

   
            
                while (  it != nrlvec.end() )
                {
                    //overwrite polynomialId if polynomialId>=2)
                    if ( (*it).getPolynomialId()>=2 || ! hurwitzMapSearchProblem_m.strictNormalization() )
                    {
                        (*it).clearPolynomialId(  );
                        (*it).clearExponent(  );
                    }
                    it++;
                }
                it = nrlvec.begin();
            
                while (  it != nrlvec.end() )
                {
                     PolynomialFactorBluePrint * polFactorBlueprint_p = extractMatchingRule( polFactorConstructRules, (*it) );
                     if ( polFactorBlueprint_p !=NULL)
                     {
                            const size_t polynomialId = (* polFactorBlueprint_p).polynomialId_m;
                            if ((*it).getValue() != NormalizationValue::infinity) 
                            {    
                                const typename IrreduciblePolTableType::IrredVecListType &  irredVecList =  irredPolTableRef.getIrredPolList((* polFactorBlueprint_p).degree_m);
                                std::vector<int>::iterator degOneVecIterator = degreeOneRootList.end();
                                if ((*it).getValue()== NormalizationValue::one )
                                {
                                    degOneVecIterator = std::find(degreeOneRootList.begin(),degreeOneRootList.end(), 1 );
                                    //std::cerr << "irredVecList[1] = " << irredVecList[1] << std::endl;
                                    assert( polynomialRing_m.evalAt( *(irredVecList[1]), field_m.Convert(1))== TPolRingTypePar::RingType::ElementType::Zero);
                                }
                                if ((*it).getValue()== NormalizationValue::zero )
                                {
                                    degOneVecIterator = std::find(degreeOneRootList.begin(),degreeOneRootList.end(), 0 );
                                    assert( polynomialRing_m.evalAt( *(irredVecList[0]), field_m.Convert(0))== TPolRingTypePar::RingType::ElementType::Zero);
                                }
                                assert(degOneVecIterator != degreeOneRootList.end() );
                                degreeOneRootList.erase( degOneVecIterator );
        
                               typename  TPolRingTypePar::Element ringElem(typename  TPolRingTypePar::Element(1) ) ;
                                assert(ringElem.getDegree()>0);
                                // todo: problem: degree von polynomial
                                assert((*it).getValue()== NormalizationValue::one || (*it).getValue()==  NormalizationValue::zero );
                               
                                ringElem.setCoeff(0,   TPolRingTypePar::RingType::ElementType::Zero );
                                if ((*it).getValue()== NormalizationValue::one )
                                    ringElem.setCoeff(0, field_m.addInv( TPolRingTypePar::RingType::ElementType::One ));
                                ringElem.setCoeff(1, TPolRingTypePar::RingType::ElementType::One );
                                if ((*it).getValue() == NormalizationValue::one )
                                    assert(ringElem == *(irredVecList[1]) );
                                if ((*it).getValue() == NormalizationValue::zero )
                                    assert(ringElem == *(irredVecList[0]) );

                                TPolFactorPowerType powfactor( ringElem, (* polFactorBlueprint_p).multiplicity_m );
                                
                                assert(polynomialId>=0 && polynomialId<polsetblueprint.size());
                                polsetblueprint[ polynomialId].push_back( powfactor );
                          
                                polSet[ polynomialId ] = polynomialRing_m.multiply( polSet[ polynomialId ], polynomialRing_m.pow(ringElem, (* polFactorBlueprint_p).multiplicity_m )  );
                            }
                            else
                            {
                                //correct shapeList: remove infinity exponent since it is invisible for univariate polynomials .
                                shapeList = shapeList.removeExponent( (* polFactorBlueprint_p).polynomialId_m, (* polFactorBlueprint_p).multiplicity_m );
                            }
                            it=nrlvec.erase(it);    
                            
                            delete polFactorBlueprint_p;
                     }
                     else
                        it++;
                }
                if ( ! hurwitzMapSearchProblem_m.strictNormalization() || nrlvec.size()==0 )
                    return true;
                if ( hurwitzMapSearchProblem_m.strictNormalization() && nrlvec.size()==1 )
                {   
                    /// todo: setzt voraus, dass bei den Normalisierungsregeln als polynomialid nur 0,1,2 zugelassen ist!
                    /// Wenn eine Normalisierungsregel übrig geblieben ist, ist das in Ordnung.
                    if ( nrlvec[0].getPolynomialId()==2 ||  nrlvec[0].getPolynomialId()==NormalizationRule::dontcare ) 
                        return true;
                }
                return false;
            }


            template <class TPolRingTypePar>
            typename TPolRingTypePar::ElementType   FiniteFieldSearch<TPolRingTypePar>::getScalarFromInt(int val)   const
            {
                typename TPolRingTypePar::ElementType   scalar(0);
                //a1.setCoeff(1, TPolynomialRing::RingType::ElementType::One );
                scalar.setCoeff(0, val );
                return scalar;
            }

            template <class TPolRingTypePar>
            void        FiniteFieldSearch<TPolRingTypePar>::removeConstantFactorsInPlace(typename FiniteFieldSearch<TPolRingTypePar>::PolSetType &  polSet)   const
            {
                for (size_t pos =0;pos<polSet.size();pos++)
                {
                typename TPolRingTypePar::ElementType & pol=polSet[pos];
                    polynomialRing_m.normalizeInPlace( pol );
                }
            }

            template <class TPolRingTypePar>
            template <typename TProduct>
            bool    FiniteFieldSearch<TPolRingTypePar>::normalizationRuleMatches(const std::vector<TProduct> & prodVec , const NormalizationRule & nr,  int degree ) const
            {             

                int baseDegree=0;
                int currDegree=0;

                for (size_t prodVecPos = 0; prodVecPos< prodVec.size(); prodVecPos++ )
                {
                        if (  nr.getPolynomialId() != NormalizationRule::dontcare &&  nr.getPolynomialId()!=prodVecPos)
                            continue;

                    for (size_t factorPos = 0; factorPos< prodVec[prodVecPos].size(); factorPos++ )
                    {

                        typename TPolRingTypePar::Element el= prodVec[prodVecPos][ factorPos ].first();
                        baseDegree = el.getDegree();
                        
                        if ( nr.getValue()== NormalizationValue::one || nr.getValue()== NormalizationValue::zero )
                        {
                            int normVal=0;
                            if ( nr.getValue()== NormalizationValue::one )
                                normVal=1;

                            if ( nr.getPolynomialId()==NormalizationRule::dontcare || nr.getPolynomialId()==prodVecPos   )
                            {
                                //std::cerr << "pos matches" << std::endl;
                                if ( polynomialRing_m.evalAt( prodVec[prodVecPos][factorPos].first() , field_m.Convert(normVal) ) == TPolRingTypePar::RingType::ElementType::Zero)
                                {
                                    //std::cerr << "value matches" << std::endl;                                    
                                    //return true;
                                    if ( nr.getExponent()==NormalizationRule::dontcare || nr.getExponent()==prodVec[ prodVecPos ][ factorPos ].second() )
                                    {
                                        //std::cerr << "exponent matches" << std::endl;

                                        return true;
                                    }
                                }
                                
                                    
                            }
                        }
    
                        currDegree += baseDegree*prodVec[prodVecPos][factorPos].second();
                    }

                    if ( nr.getValue()==NormalizationValue::infinity && currDegree< degree &&   
                        ( nr.getPolynomialId()==NormalizationRule::dontcare || nr.getPolynomialId()==prodVecPos  ) )
                    {
                        //std::cerr << "infinity: value matches" << std::endl;
                        return true;
                    }

                
                }

                return false;
            }
            
            template <class TPolRingTypePar>
            template <typename TProduct>
            bool    FiniteFieldSearch<TPolRingTypePar>::normalizationRulesMatches(const std::vector<TProduct> & prodVec  ) const
            {
                NormalizationRuleList nrl = hurwitzMapSearchProblem_m.getNormalizationRules();
                
                assert( nrl.size()<=prodVec.size() );

                bool allMatches=true;
                for (size_t pos = 0; pos < nrl.size(); pos++ )
                {
                    bool  matches  = normalizationRuleMatches( prodVec, nrl[pos], hurwitzMapSearchProblem_m.getMapDegree() )  ;
                    
                    if ( not matches )
                        return false;
                  //     std::cerr << "rule "<< pos << " matches" << std::endl;
                }
                //std::cerr << "allMatches  "  << allMatches <<  std::endl;
                //std::cerr << "--" <<  std::endl;
                
                return allMatches;
            }

            template <class TPolRingTypePar>
            void  FiniteFieldSearch<TPolRingTypePar>::renormalize(  const ShapeList & shapeList,  FiniteFieldSearch<TPolRingTypePar>::PolSetType & polSetCopy )
            {
                    #ifdef VERBOSE
                        std::cerr << "# renormalize " << std::endl;
                        #endif

                bool                  polSetMultiplicityStructureIsOk=true;

                        typename FiniteFieldSearch<TPolRingTypePar>::PolSetType secondCopy = polSetCopy;
                        
                        typename TPolRingTypePar::ElementType substpol= typename TPolRingTypePar::ElementType(1)  ;  
                        substpol.setCoeff(0,field_m.Convert(0));
                        substpol.setCoeff(1,field_m.Convert(1));

                        bool first=true;
                        do 
                        {
  
                            secondCopy[0] = polynomialRing_m.subst(polSetCopy[0],  substpol);
                            secondCopy[1] = polynomialRing_m.subst(polSetCopy[1],  substpol);
                            secondCopy[2] = polynomialRing_m.subst(polSetCopy[2],  substpol);
                            if (first)
                            {
                                assert(secondCopy[0]==polSetCopy[0]);
                                assert(secondCopy[1]==polSetCopy[1]);
                                assert(secondCopy[2]==polSetCopy[2]);
                                first=false;
                            }
                       
                            // todo: define check as filter!
                            if (  polynomialMatchesShapeStrict( secondCopy[0], polynomialRing_m, shapeList[0] ) )
                            if (  polynomialMatchesShapeStrict( secondCopy[1], polynomialRing_m, shapeList[1] ) )
                            if (  polynomialMatchesShapeStrict( secondCopy[2], polynomialRing_m, shapeList[2] ) )

                            /// todo: gcd check is only correct for three critival values.
                            if ( polynomialRing_m.fastgcd( secondCopy[0],secondCopy[1]).isConstant() )
                            if ( polynomialRing_m.fastgcd( secondCopy[1],secondCopy[2]).isConstant() )
                            if ( polynomialRing_m.fastgcd( secondCopy[0],secondCopy[2]).isConstant() )
                        
                            {
                               typedef  Product< Power<TPolRingTypePar> > TProduct;
                                std::vector <TProduct> prodVec;
                                for (size_t prodVecPos =0; prodVecPos < secondCopy.size(); prodVecPos++)
                                {
                                       TProduct prod = FactorizerType::extendedFactorPolynomial<TPolRingTypePar>( secondCopy[prodVecPos], polynomialRing_m );
                                       prodVec.push_back(prod);
                                }
                              
                                // now check if prod_i mathes the normalization Rule for the   polynomial i.
                             
                                if ( ! normalizationRulesMatches(  prodVec ) )
                                    continue;
    
                               

                                #ifdef VERBOSE
                                std::cerr << "substpol " << substpol<< std::endl;
                                #endif

                                for ( size_t polpos=0; polpos<polSetCopy.size(); polpos++)
                                    secondCopy[polpos]= polynomialRing_m.subst(  polSetCopy[polpos], substpol );

                                std::vector< typename CoeffRingType::ElementType >    scalingRelations ;

                                if ( ! computeScalingFactors(secondCopy, scalingRelations))
                                {
                                    #ifdef VERBOSE
                                    std::cerr << "failed to compute scaling factors" << std::endl;
                                    #endif
                                    continue;
                                }
                            
                                // check B-A==C: (lambda =1 !);
                                assert( secondCopy[2] == polynomialRing_m.add( secondCopy[1], polynomialRing_m.addInv( secondCopy[0]) ) ) ;
 
                                #ifdef VERBOSE
                                for  (size_t spos=0;spos<scalingRelations.size(); spos++)
                                    std::cerr << "scaling factor [ " << spos << "] = " << scalingRelations[spos] << std::endl;
                                #endif
                                
                                if (     minimalPolynomials_m.size()>0)
                                {   
                                   
                                    // optional:
                                    // removeConstantFactorsInPlace( polSetCopy );
                                    // normalizeInPlace( polSetCopy );
                                    // normalize
                                    // then compute and apply scaling factors.
                                    //
                                    // removeConstantFactorsInPlace( polSetCopy );

                                    for (size_t minPolPos=0; minPolPos<minimalPolynomials_m.size(); minPolPos++)
                                    {   
                                        /// todo: eigentlich weiss nur ein Ring, ob ein Element=0 ist oder nicht ...bisher falsch umgesetzt...
                                        if (not polynomialRing_m.evalAt(minimalPolynomials_m[minPolPos], 
                                                                        field_m.multiply(scalingRelations[minPolPos+1] ,
                                                                                                    field_m.multInv(scalingRelations[0] )
                                                                                                )
                                                                        ).isZero() 
                                            )
                                        {
                                            polSetMultiplicityStructureIsOk = false;
                                            break;
                                        }
                                    }
                                }
                            
                                if (! polSetMultiplicityStructureIsOk)   continue;       
                    
                        
                                #pragma omp critical    
                                {

                                    //localCounter++ ;
                                    #ifdef VERBOSE
                                    std::cerr << "found solution : " << std::endl;
                                    std::cerr << secondCopy[0] << std::endl;
                                    std::cerr << secondCopy[1] << std::endl;
                                    std::cerr << secondCopy[2] << std::endl;
                                    #endif
                                
                                    /// ensure that all polynomial coefficient lists have the same size:
                                    /// @todo: this should be a part of polSetObject 
                                    /// @todo : use everywhere field_m.getZero() and field_m.getOne() for 0 and 1 or similar. 
                                    /*for (size_t polId=0;    polId< secondCopy.size();polId++)
                                    {
                                        secondCopy[polId].setDegree( shapeList.getDegree() );
                                        //for  (size_t missingcoeffPos = secondCopy[polId].size();  missingcoeffPos<  shapeList.getDegree();      missingcoeffPos++)
                                        //       secondCopy[polId].push_back( field_m.getZero() );                               
                                    }*/
        
                                    ////
                                      PolynomialSet<TPolRingTypePar>  polSetObject=PolynomialSet<TPolRingTypePar>( hurwitzMapSearchProblem_m,
                                                                                 searchOptions_m,
                                                                                 polynomialRing_m,
                                                                                 secondCopy  
                                                                                );

                                    outputHandler_m->print( polSetObject );
                                    solutionCandidates_m.push_back( secondCopy );
                                    
                                }
                                if ( not searchOptions_m.allNormalizations() )
                                    break;
                            }

                        }
                        while ( substpol.nextInPlace(field_m) );
            }

            ///@note es ist zwar prinzipiell egal, ob das Zweite Polynom als Nullstelle inf hat, aber um konstante Faktoren be den Polynomen zu vermeiden,
            /// sollte das schon so gehandhabt werden.            
            template <class TPolRingTypePar>
            void  FiniteFieldSearch<TPolRingTypePar>::last_search_level(  const ShapeList & shapeList,  FiniteFieldSearch<TPolRingTypePar>::PolSetType & polSet )
            {
                //std::cerr << "last_Search_Level" << std::endl;  

                size_t mid = shapeList.size()-2; 

         
                size_t nonzeroesNum =  field_m.getCharacteristic()-1 ;

                if (mid>nonzeroesNum) 
                //    return 0;
                return;

                if (   searchOptions_m.dryRun() )
                    assert( false ); //should not happen
                //    return 1;

                //CounterType localCounter=0;
                #ifdef SAFE
                    assert( polynomialMatchesShape( polSet[0], polynomialRing_m, shapeList[0] ) );
                    assert( polynomialMatchesShape( polSet[1], polynomialRing_m, shapeList[1] ) );
                #endif

             ScalarCombinationType    nonzeroes;
                /// todo: the nonzeroes should be provided by the field_m (e.g. a galois field) 
                nonzeroes = ScalarCombinationType( field_m.getCharacteristic()-1 );

                for (int num=1; num< field_m.getCharacteristic(); num++)
                    nonzeroes[num-1]=num;
            
                ScalarCombinationType   combination(mid);
              

                for (size_t i=0; i< mid ; i++)
                    combination[i] = i;
 
                 int tmpCombSize = 32;

                std::vector<ScalarCombinationType > tmpCombinations;
                tmpCombinations =  std::vector<ScalarCombinationType >(tmpCombSize);
                

                /// checking all examples means walk through all nonzero combinations of size mid. 
                /// To parallelize the computations a portion of tmpCombSize=64 combinations is computed and then distributed to all local CPU's by OpenMP
                /// The first idea was to compute a list of all combinations and then parallelize, but that could result in a problem since for greater characteristic and 
                /// a lot of branch points the size of the combinations explodes combinatorically
                /// parallelizing the computation across a computer cluster needs more serious thought.
                bool done = false;
                while (!done)
                {
                    for (int counter=0; counter<tmpCombSize; counter++)
                    {
                        tmpCombinations[counter] = combination;
                        done = ! naive_next_combination_vec( nonzeroesNum, mid , combination);
                        if (done)
                        {
                            tmpCombinations.resize(counter+1);
                            break;
                        }
                    }
                    tmpCombSize = tmpCombinations.size();                

                    // es is nicht sehr effizient, die innerste Schleife zu parallelisieren 
                    // wenn die Anzahl der Verzweigungen nicht zu gross ist, sollte die Größe von nonzeroCombinationVector_m nicht explodieren.
                    //#ifdef OPENMP
                    //    #pragma omp parallel for
                    //#endif

                    for (int combinationPos = 0; combinationPos< tmpCombSize ; combinationPos++)      
                    {

                        typename FiniteFieldSearch<TPolRingTypePar>::PolSetType polSetCopy;
                       
                        polSetCopy = polSet;
                        
                        bool polSetMultiplicityStructureIsOk=true;
                        
                        combination = tmpCombinations[ combinationPos ] ;

                        for (size_t pos=0; pos<mid; pos++)
                        {
                            typename TPolRingTypePar::Element secondSummand =  polynomialRing_m.multiply( getScalarFromInt( nonzeroes[combination[ pos ] ] ),
                                                                                                        polSetCopy[0] 
                                                                                                    ) ;
            
                            // wieso eigentlich additiveInverse?
                            //polSetCopy[pos+2] = polynomialRing_m.addInv( polynomialRing_m.add (polSetCopy[0],  secondSummand )      );
                            polSetCopy[pos+2] =  polynomialRing_m.add (polSetCopy[1],  secondSummand )  ;

                            if (searchOptions_m.logStructure() )
                            {
                                Shape::ShapeRepType tmpShapeRep =  computeMultiplicityStructure<TPolRingTypePar,Shape::ShapeRepType>( polSetCopy[2], polynomialRing_m ) ;
                                assert(tmpShapeRep.size()>0);
                                Shape logShape = Shape(tmpShapeRep);
                                std::string key = logShape.toString();
    
                                #pragma omp critical    
                                {
                                    if ( logShape.getDegree() != shapeList.getDegree() )
                                    {
                                        //std::cerr << "logged shape degree " << logShape.getDegree()<< std::endl;
                                        //std::cerr << "shapeList degree " << shapeList.getDegree()<< std::endl;
                                        //std::cerr << " unextected shape degree  "  << logShape << std::endl << polSetCopy[2] << std::endl;
                                        //getchar();
                                    }
                                    if ( singleStructureHashmap_m.find(key) == singleStructureHashmap_m.end() )
                                            singleStructureHashmap_m.insert(std::pair<std::string, uint64_t> (key,1));
                                        else 
                                            singleStructureHashmap_m[key]++;
                                    if ( fullStructureHashmap_m.find(key) == fullStructureHashmap_m.end() )
                                            fullStructureHashmap_m.insert(std::pair<std::string, uint64_t> (key,1));
                                        else 
                                            fullStructureHashmap_m[key]++;
                                }
                            } 


    
                            // check multiplicityStructure of polSetCopy[pos+2]
                            if ( ! polynomialMatchesShape( polSetCopy[pos+2], polynomialRing_m, shapeList[pos+2] ) )
                            {
                                polSetMultiplicityStructureIsOk=false;
                                break;
                            }
                            
                        }

                        
        
                    
                        if (!polSetMultiplicityStructureIsOk) continue;
                        //#pragma omp critical
                        {
                            //localCounter++;
                        }

                         // passe W_inf an.
                         polSetCopy[0]=   polynomialRing_m.addInv(polynomialRing_m.multiply( getScalarFromInt( nonzeroes[combination[ 0 ] ] ),
                                                                                            polSetCopy[0] 
                                                                                        ) 
                                                                  );


                        #ifdef VERBOSE
                        std::cerr << "found solution candidate" << std::endl;
                        #endif
                        // check multiplicityStructure of polSetCopy[0] and polSetCopy[1] and gcd (polSet[0], polSet[0]) 
                        // because in case preudoIrreduciblePolynomial lists are used, the multiplicityStructure may be wrong 
                        // or gcd (polSet[0], polSet[0]) not a constant! Both conditions are mandatory!
                        
                        // BUT concider following: 
                        // when checking polynomials fast, the characteristic has to be of sufficient size depending on the destination polynomial shape
                        // otherwise a slower factorization is required to get the correct result
                        // since polSetCopy[0] and polSetCopy[1] may have a shape that is not testable in a fast way
                        // a factorization is required here for a correct result. 
                        
                        if ( ! polynomialMatchesShapeStrict( polSetCopy[0], polynomialRing_m, shapeList[0] ) )
                         continue;

                      

                        if ( ! polynomialMatchesShapeStrict( polSetCopy[1], polynomialRing_m, shapeList[1] ) )
                              continue;

                      
                        
                            typename TPolRingTypePar::ElementType gcdf = polynomialRing_m.fastgcd(polSetCopy[0],polSetCopy[1]);
                            if ( !gcdf.isConstant() )
                                 continue;
                        
 
                                    
                                        
                        #ifdef VERBOSE
                        std::cerr << " pre solution (not normalized) " << std::endl;
                        std::cerr << polSetCopy[0] << std::endl;
                        std::cerr << polSetCopy[1] << std::endl;
                        std::cerr << polSetCopy[2] << std::endl;
                        #endif


                        renormalize(shapeList,polSetCopy);
                      
                    }// end for 
               } //end while not done
                //return localCounter ;
                return;
            }


            //weitere idee: die sortedPolFactorConstructRules so sortieren, dass oben die Liste mit dem größten grad ist und nur das erste mal parallelisieren.
            template <class TPolRingTypePar>
            void    FiniteFieldSearch<TPolRingTypePar>::third_search_level( const std::vector< std::pair<int,const BPVecTYPE * > > &  sortedPolFactorConstructRules,
                                                                            const std::vector<int> &    degreeOneRootList,
                                                                            const ShapeList & shapeList,
                                                                            std::vector< std::pair<int, const BPVecTYPE * > >::const_reverse_iterator    mapIterator,
                                                                            FiniteFieldSearch<TPolRingTypePar>::PolSetType & polSet,
                                                                            mpz_t tmpCounterPar, mpz_t localCounter)
            {   
                    //std::cerr << "third_Search_Level" << std::endl;
		   mpz_t tmpCounter;
            
    
                    if ( mapIterator==sortedPolFactorConstructRules.rend() )
                    {
                        //#pragma omp critical
                       

                        if (   searchOptions_m.dryRun() )
                        {
			     
			      mpz_init( tmpCounter );
			      mpz_set(tmpCounter, tmpCounterPar);

			    {
			      mpz_init( localCounter );
			    }
                            mpz_set( localCounter, tmpCounter );
                            size_t mid = shapeList.size()-2; 
                            mpz_t tmpMpz;
                            mpz_init( tmpMpz );
                        
                            for (size_t tmp =field_m.getCardinality()-1; tmp> field_m.getCardinality()-1 - mid ;tmp--)
                            {
                                 mpz_set_ui( tmpMpz,tmp );
                                 mpz_mul(localCounter, localCounter, tmpMpz);
                            }
                            return;
                            
                        }
                        last_search_level(shapeList, polSet);
                        return  ;
                    }
                    else
                    {
                     
                       int degree = (*mapIterator).first;       
                         
                       const BPVecTYPE  & polFactorConstRuleList = *((*mapIterator).second);

                        mapIterator++;

                       size_t mid =   polFactorConstRuleList.size()  ;

                       size_t   irredListSize;

                        if (   searchOptions_m.dryRun() )
                        {
			      
			      mpz_init( tmpCounter );
			      mpz_set(tmpCounter, tmpCounterPar);
			      
                                mpz_t irredListSizeMpz;
                                mpz_init( irredListSizeMpz );
                                getIrredCount( degree, field_m.getCardinality(), irredListSizeMpz );
                                
                                 if (degree==1)
                                    {
                                        mpz_set_ui(irredListSizeMpz, degreeOneRootList.size() );
                                    }

                                irredListSize = mpz_get_ui(irredListSizeMpz);
                                getIrredCount( degree, field_m.getCardinality() ); // will assert if irredListSize to big

                                //int irredListSizeDebug = mpz_get_ui(irredListSizeMpz);
                                //assert(irredListSizeDebug==irredListSize);


                                #ifdef VERBOSE
                                std::cerr << "degree " << degree << std::endl;
                                std::cerr << "irredListSize " << irredListSize << std::endl;
                                #endif

                                mpz_t midMpz;
                                mpz_init( midMpz );
                                mpz_set_ui(midMpz,mid);

                                if ( mpz_cmp(midMpz,irredListSizeMpz)>0 ) 
                                {
                                    mpz_init( localCounter );
                                    return;
                                }
                                
                                mpz_t tmpMpz;
                                mpz_init( tmpMpz );
                                mpz_t tmpCmpMpz;
                                mpz_init( tmpCmpMpz );
                                mpz_set(tmpMpz,irredListSizeMpz);
                               
                                mpz_add(tmpCmpMpz,tmpMpz,midMpz);

                                mpz_t minusOneMpz;
                                mpz_init( minusOneMpz );
                                mpz_set_si(minusOneMpz,-1);
                            
                                while ( mpz_cmp(tmpCmpMpz, irredListSizeMpz)>0 )
                                {
                                    mpz_mul( tmpCounter, tmpCounter,tmpMpz ); 
                                    mpz_add(tmpMpz,tmpMpz,minusOneMpz);
                                    mpz_add( tmpCmpMpz,tmpMpz,midMpz); // tmpCmpMpz=tmp+mid
                                  
                                    
                                };
                                /*for (size_t  tmp = irredListSize; tmp>irredListSize-mid; tmp-- )
                                {
                                    mpz_set_ui(tmpMpz,tmp);
                                    mpz_mul(tmpCounter, tmpCounter,tmpMpz);   
                                }*/

                                mpz_t tmpCounterCopy;
                                mpz_init(tmpCounterCopy);
                                mpz_set(tmpCounterCopy,tmpCounter);

                                third_search_level(sortedPolFactorConstructRules, degreeOneRootList, shapeList, mapIterator, polSet, tmpCounterCopy, localCounter);        
                                return;
                        }
 
    
                      const IrreduciblePolTableType * irredPolTablePtr= getIrredPolTableConstPtr( field_m.getCharacteristic() );    
                        const typename IrreduciblePolTableType::IrredVecListType &  irredVecList =  irredPolTablePtr->getIrredPolList(degree);
                       

                        /// note: irredListSize, mid , combination must have the same type for working naive_next_combination!
                       
                        irredListSize = irredVecList.size();
                  

                        if (degree==1)
                            irredListSize= degreeOneRootList.size();

           
                        if ( mid>irredListSize ) 
                        {
                            #ifdef VERBOSE
                            std::cerr << "degree " << degree << std::endl;
                            std::cerr << "irredListSize " << irredListSize << std::endl;
                            std::cerr << "polFactorConstRuleListSize " << mid << std::endl;
                            #endif
                            return ;//return 0;
                        }
 
                        ScalarCombinationType combination( mid );
                            for (size_t i=0; i< mid ; i++)
                                combination[i] = i;

                        ///todo: bei dry run: zähle einfach Kombinationen und permutationen anstatt alle durchzulaufen. DONE
 
                        int tmpCombSize = 64;
                        std::vector<ScalarCombinationType > tmpCombinations(tmpCombSize);

                        bool done = false;
                        
                        while (!done)
                        {
                            for (int counter=0; counter<tmpCombSize; counter++)
                            {
                                tmpCombinations[counter] = combination;
                                done = ! naive_next_combination_vec( irredListSize, mid , combination);
                                if (done)
                                {
                                    tmpCombinations.resize(counter+1);
                                    break;
                                }
                            }
                        
                            tmpCombSize = tmpCombinations.size();

                        
                            //#ifdef OPENMP
                            //    #pragma omp parallel for  reduction(+: localCounter)
                            //#endif
                            #ifdef OPENMP
                                #pragma omp parallel for 
                            #endif

                            // idee: laufe alle combinations durch; integer in der combination sind indices für die irreducible liste...
                            // falls es sich um deg1 irred handelt, doppelte indirektion:  irredlist[ degreeOneRootList[pos] ]
                            for (int combinationPos = 0; combinationPos< tmpCombSize ; combinationPos++)      
                            {
                                // for each permutation...
                                size_t permutation[ mid ];
                                for (size_t i=0; i< mid ; i++)  permutation[i] = tmpCombinations[combinationPos] [i];
                                
                                do 
                                {
                                    if ( not searchOptions_m.dryRun() )
                                    {
                                           std::vector<typename TPolRingTypePar::ElementType> polSetCopy( polSet.begin(), polSet.end() ) ;
                                        assert( polSet.size()>0 );
                                        assert( polSetCopy.size()>0 );
                                        // here: construct polynomials: apply all rules from polFactorConstRuleList.
                                        for (size_t rulePos = 0 ;  rulePos < polFactorConstRuleList .size(); rulePos++ )  
                                        {
                                            #ifdef SAFE
                                                assert( polFactorConstRuleList[rulePos].polynomialId_m < polSetCopy.size() );
                                                assert(  polFactorConstRuleList[rulePos].polynomialId_m >= 0 );
                                            #endif
                                            typename TPolRingTypePar::ElementType & polRef = polSetCopy[ polFactorConstRuleList[rulePos].polynomialId_m ];
                                            if (degree==1)
                                                polRef =  polynomialRing_m.multiply(  polRef ,polynomialRing_m.pow( *(irredVecList[ degreeOneRootList[permutation[rulePos ]] ]), polFactorConstRuleList[rulePos].multiplicity_m)  )  ;
                                            else
                                                polRef =  polynomialRing_m.multiply(  polRef ,polynomialRing_m.pow( *(irredVecList[permutation[rulePos ] ]), polFactorConstRuleList[rulePos].multiplicity_m)  )  ;
                                            
                                        }
                                          third_search_level(sortedPolFactorConstructRules, degreeOneRootList, shapeList, mapIterator, polSetCopy, tmpCounter,localCounter );    
                                    }
                                    else
                                    {
                                        third_search_level(sortedPolFactorConstructRules, degreeOneRootList, shapeList, mapIterator, polSet, tmpCounter,localCounter);    
                                    }
    
                                } while (next_permutation (permutation, permutation + mid ) );
                                
                            
                            }// end for   
                        }// end while !done
                     }// end not all construction rule processed.
                    return;
            }

            
            template <class TPolRingTypePar>
            void    FiniteFieldSearch<TPolRingTypePar>::second_search_level(std::list< PolynomialFactorBluePrint> polFactorConstructRules)
            {

                #ifdef VERBOSE
                std::cerr << "second_Search_Level" << std::endl;
                #endif
                std::list< PolynomialFactorBluePrint>::iterator it;
                
                /*for (it=polFactorConstructRules.begin(); it != polFactorConstructRules.end(); it++)
                {
                    std::cout << (*it) << std::endl;
                }*/
                
                PolSetBlueprintType polsetblueprint;
                polsetblueprint.resize( hurwitzMapSearchProblem_m.getPolSetSize() );

                #ifdef VERBOSE
                std::cerr << " polFactorConstructRules.size() = " << polFactorConstructRules.size() << std::endl;
                #endif
                std::vector<int>  degreeOneRootList;

                for (int num= 0; num< field_m.getCharacteristic();num++)
                     degreeOneRootList.push_back(num);

                ShapeList modifiedShapeList = hurwitzMapSearchProblem_m.getShapeList();

                        
                PolSetType  polSet(modifiedShapeList.size(), TPolRingTypePar::ElementType::getOne()  );
                //assert( TPolRingTypePar::ElementType::One.isOne() );
                assert( TPolRingTypePar::ElementType::getOne().isOne() );
                assert( polSet.size()>0 ); 
                if ( processNormalizationRules( polsetblueprint, polFactorConstructRules, degreeOneRootList, modifiedShapeList, polSet) )
                {
                    // sort 
                    for (it=polFactorConstructRules.begin(); it != polFactorConstructRules.end(); it++)
                    {
                        #ifdef VERBOSE
                        std::cerr << (*it) << std::endl;
                        #endif
			DebugLogger::logStream() << (*it) << std::endl;
                    }
                
                    std::list< PolynomialFactorBluePrint>::iterator it;
                    DegSortedConstructionRuleTableType   sortedPolFactorConstructRules;

                    for (it=   polFactorConstructRules.begin();it!=polFactorConstructRules.end();it++)
                    {
               
                        if (sortedPolFactorConstructRules.find( (*it). degree_m )==sortedPolFactorConstructRules.end())
                        {
                            sortedPolFactorConstructRules[(*it). degree_m]=   std::vector< PolynomialFactorBluePrint>();
                        }
                        sortedPolFactorConstructRules[(*it). degree_m].push_back( *it );
                    }
 
          
                    std::vector< std::pair<int, const BPVecTYPE *> >  sortedPolFactorConstructRulesVec;

                DegSortedConstructionRuleTableType::iterator dit = sortedPolFactorConstructRules.begin();
                while (dit!=sortedPolFactorConstructRules.end())
                {

                    int degree = (*dit).first;
                    const BPVecTYPE * vecbp = &((*dit).second);
                    std::pair<int, const BPVecTYPE* > pp (degree, vecbp);
                     //sortedPolFactorConstructRulesVec.push_back( std::pair<int,BPVecTYPE > (degree, vecbp ) );
                    sortedPolFactorConstructRulesVec.push_back( pp );
                    dit++;
                }

                    std::vector< std::pair<int, const BPVecTYPE * > >::const_reverse_iterator mapIterator = sortedPolFactorConstructRulesVec.rbegin();
                
                    mpz_t tmpCounter;
                    mpz_init(tmpCounter); 
                    mpz_set_ui(tmpCounter,1);
                    mpz_t counterMpz;
                    third_search_level( sortedPolFactorConstructRulesVec , degreeOneRootList, modifiedShapeList, mapIterator, polSet ,tmpCounter, counterMpz );

                    
                    #pragma omp critical
                    {
                        if (searchOptions_m.dryRun())
                        {
                            mpz_add( counter_m, counter_m, counterMpz );
                            if ( mpz_cmp(counter_m , counterMod_m )   >0)
                            {
                                #ifdef VERBOSE
                                //std::cerr << " counter_m = " << counter_m << " ! " << std::endl;
                                //std::cerr << " counter_m = " << counter_m << " ! " << std::endl;
                                
                                #endif
                                mpz_mul_ui(counterMod_m,counterMod_m,10);
                                //counterMod_m *= 10;
                            }
                        }
                    }
                }
                else
                    assert(false);
            }

            template <class TPolRingTypePar>
            void    FiniteFieldSearch<TPolRingTypePar>::first_search_level( std::list< PolynomialFactorBluePrint>     polfactorBlueprintList,  
                                      std::list< PolynomialFactorBluePrint>   polFactorConstructRules)
            {    
               
                if (polfactorBlueprintList.size()==0)
                {
                    second_search_level(polFactorConstructRules);
                    printStructureMap(singleStructureHashmap_m);
                    singleStructureHashmap_m.clear();
                }
                else
                {
                    PolynomialFactorBluePrint bp = polfactorBlueprintList.front();
    
                    polfactorBlueprintList.pop_front();

                    std::vector<int> partition ( bp.degree_m ,1);
                    do 
                    {
                        std::list< PolynomialFactorBluePrint>  polFactorConstructRulesCopy=polFactorConstructRules,
                     
                        currentPolFactorConstructRules  = createPolFactorConstructionRules( partition, bp.multiplicity_m ,bp.polynomialId_m);

                        polFactorConstructRulesCopy.splice( polFactorConstructRulesCopy.end(), currentPolFactorConstructRules);
                        first_search_level(polfactorBlueprintList, polFactorConstructRulesCopy );
                    }
                    while (next_partition_desc( & partition )) ;
                }
           
            }
           
            template <class TPolRingTypePar>
            void    FiniteFieldSearch<TPolRingTypePar>::run()   
            {
                //counter_m = 0;  
                mpz_set_ui(counter_m,0);
                #ifdef VERBOSE
                std::cerr << "run()" << std::endl;
                #endif
                std::list< PolynomialFactorBluePrint>    polfactorBlueprintList = createPolFactorBlueprintList(  ) ;

                std::list< PolynomialFactorBluePrint> polFactorConstructRules;

                if (! searchOptions_m.dryRun() )
                    outputHandler_m->startOutput();
                
                first_search_level(polfactorBlueprintList, polFactorConstructRules);

                if (! searchOptions_m.dryRun() )
                    outputHandler_m->finishOutput();      
            }

              template <class TPolRingTypePar>
              FiniteFieldSearch<TPolRingTypePar>::~FiniteFieldSearch()
            {
                delete outputHandler_m;

            }

}


