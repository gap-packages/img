
#pragma once


#include "Shape.h"
#include "NormalizationRules.h"

#include "hmfTypedefs.h"
#include "combinatorics/partition.h"
#include "IrreduciblePolTable.h"
#include <map>
#include <ext/hash_map>
#include <list>
#include "combinatorics/combinations.h"

#include "PolSetOutputHandlers.h"



#include "FactorPolynomialWrapper.h"


// todo: wenn degree polSet[2]<12, dann polynom zusammen mit der Multiplizitätsstruktur ausgeben.


namespace RationalMapSearch
{

    

    class HMSProblem
    {

        public:
            typedef std::vector<int> PolynomRepType;


    private:
          ShapeList               shapeList_m;
          NormalizationRuleList   normRuleList_m;
          bool strictNormalization_m;
        
          std::vector<PolynomRepType>     minimalPolynomials_m;
        
    public:

        HMSProblem( const ShapeList & sl ,
                        std::vector<PolynomRepType> minimalPolynomials =  std::vector<PolynomRepType>() ) : shapeList_m(sl),
                                                                                normRuleList_m( NormalizationRuleList::constructDefault(shapeList_m) ),
                                                                                minimalPolynomials_m(minimalPolynomials)
        {
            checkConsistency();
        };

        HMSProblem( const ShapeList &sl, const NormalizationRuleList& nrl,
                            std::vector<PolynomRepType> minimalPolynomials =  std::vector<PolynomRepType>() ): shapeList_m(sl), 
                                                                                    normRuleList_m(nrl) ,
                                                                                     minimalPolynomials_m(minimalPolynomials)
        {
            checkConsistency();
        };

      

        void checkConsistency() const
        {
            assert( minimalPolynomials_m.size()==0 || minimalPolynomials_m.size()== shapeList_m.size()-3 );
        }
        
        const ShapeList&  getConstShapeListRef() const {  return shapeList_m; };

        ShapeList  getShapeList() const         {  return shapeList_m; };

        size_t      getShapeListSize() const         {  return shapeList_m.size(); };

        const ShapeList&   getShapeListConstRef() const         {  return shapeList_m; };

        NormalizationRuleList  getNormalizationRuleList() const {  return normRuleList_m; };

        NormalizationRuleList  getNormalizationRules() const {  return normRuleList_m; };

        int     getPolSetSize() const   {   return shapeList_m.size();  };

        bool    strictNormalization()     const     {  return normRuleList_m.strictNormalization();   };

        int getMapDegree() const
        {
            return shapeList_m.getDegree() ;
        }
        int getLowerCharacBound()   const
            {
                return shapeList_m.computeLowerCharacBound();
            }

        const std::vector<PolynomRepType> &    getMinimalPolynomials() const {  return minimalPolynomials_m;    };

    };

   

    class SearchOptions
    {

      private:
            bool dryRun_m;
            bool logStructure_m;
    
            bool strictNormalization_m;

            OutputMode  outputMode_m;

            bool bAllNormalizations_m;

        public:
            SearchOptions():    dryRun_m(false),
                                logStructure_m(false),
                                outputMode_m( OutputMode::defaultOutput ), //   requires c++0x
                                bAllNormalizations_m(false)
            {};

             SearchOptions(bool dryRun, 
                            bool logStructure, 
                            bool strictNormalization,
                            OutputMode outputMode=OutputMode::defaultOutput ): //requires c++0x
                                            dryRun_m(dryRun),
                                            logStructure_m(logStructure),
                                            strictNormalization_m(strictNormalization),
                                            outputMode_m(outputMode),
                                            bAllNormalizations_m(false)
            {};

      
            inline OutputMode outputMode()  const {   return outputMode_m;    }
            inline bool logStructure() const {   return logStructure_m;  }
            inline  bool dryRun()       const {   return dryRun_m;  }
    
            inline  bool strictNormalization()       const {   return strictNormalization_m;  }

            inline  bool allNormalizations()       const {   return bAllNormalizations_m;  }
    
            void    print(std::ostream & os) const
            {
                os << "SearchOptions( dryRun: " <<  dryRun() << ", logStructure " << logStructure() << ")"<< std::endl;
            }
        
    };

    template <class TPolRingTypePar>
    class PolynomialSet
    {
        public:
            typedef TPolRingTypePar PolynomialRingType;

            typedef typename TPolRingTypePar::ElementType   PolynomialType;

            typedef typename TPolRingTypePar::ElementType   ElementType;
    
            typedef  std::vector< PolynomialType >     PolSetType;
        
            

        private:
            HMSProblem      hmsProblem_m;
            SearchOptions   searchOptions_m;

            PolSetType  polSet_m;
        
            const PolynomialRingType & polynomialRing_m;
            
        public:
            PolynomialSet  ( HMSProblem hms, 
                            SearchOptions so, 
                            const PolynomialRingType & polRing,
                            PolSetType  & polSet  ): hmsProblem_m(hms),    
                                                     searchOptions_m(so),
                                                    polynomialRing_m(polRing),
                                                    polSet_m(polSet)
            {}

        typename TPolRingTypePar::ElementType operator[](size_t pos)
        {
            if ( pos<polSet_m.size() )
                return polSet_m[pos];
            assert(false);
        }

        const typename TPolRingTypePar::ElementType &  operator[](size_t pos) const
        {
            if ( pos<polSet_m.size() )
                return polSet_m[pos];
            assert(false);
        }
        
        int getDegree() const 
        {
            return (hmsProblem_m.getShapeList().getDegree() );
        }

        const PolynomialRingType & getRing()    const
        {
            return polynomialRing_m;
        }

        inline size_t size()   const
        {       
            return polSet_m.size();
        }
        
    };

    // todo: beim normalisieren wird nicht der höchste Exponent ausgewählt...
    
    // note : idea to test the FiniteFieldSearch: templatize with level-1, level2- and last-stage-functions!
    // preconditions:  the first or the second Shape-List should contain a degree 1-Factor to apply the infinityNormalization.  
    // currently applying the infinityNormalization to the first or second polynomial is mandatory for correct work.
    /// todo: es gibt eine Abhängigkeit zwischen Shape und minCharacteristic->schaue im Macaulay2-code nach, was es war.
    
    template <class TPolRingTypePar>
    class FiniteFieldSearch
    {
      public:
         typedef std::map<std::string, int64_t>                         StructureHashMapType;
         typedef   std::map< unsigned int , const  TPolRingTypePar* >  PolRingHashTableType;


        typedef     IrreduciblePolTable< TPolRingTypePar > IrreduciblePolTableType;

        typedef   std::map< unsigned int ,  IrreduciblePolTableType* >  IrredPolHashTablesType;  


        typedef  FLINTFactorPolynomial  FactorizerType;


        //typedef    typename PseudoIrreducibleTest< IrreduciblePolTableType >::IrredTesterFktType          IrredTesterFktType;
        typedef    typename PseudoIrreducibleExpensiveTest< IrreduciblePolTableType >::IrredTesterFktType          IrredTesterFktType;

    
        typedef std::vector<PolynomialFactorBluePrint > BPVecTYPE;

        typedef    std::map< int, BPVecTYPE  >      DegSortedConstructionRuleTableType;

        typedef std::vector<size_t>    ScalarCombinationType;

   
        std::vector<typename TPolRingTypePar::ElementType>   minimalPolynomials_m;

        protected:
 
            static PolRingHashTableType    polRingHashTable_m;
    
            static IrredPolHashTablesType    irredPolHashTable_m;
    
            StructureHashMapType        singleStructureHashmap_m;
            StructureHashMapType        fullStructureHashmap_m;

               template <typename TProduct>
            bool normalizationRulesMatches(const std::vector<TProduct> & prodVec  ) const;

   template <typename TProduct>
            bool normalizationRuleMatches(const std::vector<TProduct> & prodVec , const NormalizationRule & nr,  int degree ) const;


        public:

            typedef     TPolRingTypePar TPolRingType;
            typedef     typename TPolRingTypePar::CoeffRingType CoeffRingType;
    
            typedef   std::pair<typename TPolRingType::Element , uint>         TPolFactorPowerType;

            // will contain a list and eachirredTable_m list will contain a list of polynomial roots over fp.
            typedef std::vector< std::vector< TPolFactorPowerType > >     PolSetBlueprintType;

            typedef std::vector<typename TPolRingTypePar::ElementType >     PolSetType;

        private:
 
            std::list<PolSetType>   solutionCandidates_m;

           const    HMSProblem &     hurwitzMapSearchProblem_m;
            const   SearchOptions   searchOptions_m;
            
             IOutputHandler< PolynomialSet<TPolRingTypePar> >* outputHandler_m;

            int     characteristic_m;

           const TPolRingType &  polynomialRing_m;
    
          const typename TPolRingType::RingType & field_m;

            typedef int64_t CounterType;
            mpz_t   counter_m;        ///todo: use mpz integers.        
            mpz_t   counterMod_m;        ///todo: use mpz integers.      


           typename TPolRingTypePar::ElementType    getScalarFromInt(int val)   const;
        
            std::list< PolynomialFactorBluePrint>     createPolFactorBlueprintList(  );

            static const TPolRingType & getPolRingRef(uint cardinality)
            {
                  if ( polRingHashTable_m.find(cardinality ) == polRingHashTable_m.end() )
                  {
                        const typename TPolRingType::RingType* field = new  typename TPolRingType::RingType(cardinality,0);
                        const   TPolRingType * ring = new TPolRingType(*field);
                        FiniteFieldSearch::polRingHashTable_m[cardinality] = ring;
                  }
                  return *(polRingHashTable_m[cardinality]);
            }
            static const TPolRingType & getPolRingRef(uint cardinality, uint generator)
            {
                  if ( polRingHashTable_m.find(cardinality ) == polRingHashTable_m.end() )
                  {
                        const typename TPolRingType::RingType* field = new  typename TPolRingType::RingType(cardinality,0,generator);
                        const   TPolRingType * ring = new TPolRingType(*field);
                        FiniteFieldSearch::polRingHashTable_m[cardinality] = ring;
                  }
                  assert( polRingHashTable_m[cardinality]->getGenerator()==generator);
                  return *(polRingHashTable_m[cardinality]);
            }

            static const void  setPolRingRef(const TPolRingType * polRing)
            {

                 polRingHashTable_m[polRing->getCoeffRing().getCardinality() ] =  polRing;
                 return;
            }


           typename FiniteFieldSearch::IrreduciblePolTableType * getIrredPolTablePtr(uint cardinality)
            {
                  if ( irredPolHashTable_m.find(cardinality ) == irredPolHashTable_m.end() )
                  {
                         IrreduciblePolTable<TPolRingTypePar>* irredPolTable = new   IrreduciblePolTable<TPolRingTypePar>( getPolRingRef(cardinality)  ,  getIrredTester().first, getIrredTester().second );
                        FiniteFieldSearch::irredPolHashTable_m[cardinality] = irredPolTable;
                  }
                  return ( irredPolHashTable_m[cardinality] );
            }

            const typename FiniteFieldSearch::IrreduciblePolTableType * getIrredPolTableConstPtr(uint cardinality)
            {
                  if ( irredPolHashTable_m.find(cardinality ) == irredPolHashTable_m.end() )
                  {
                         IrreduciblePolTable<TPolRingTypePar>* irredPolTable = new   IrreduciblePolTable<TPolRingTypePar>( getPolRingRef(cardinality)  , getIrredTester().first, getIrredTester().second );
                        FiniteFieldSearch::irredPolHashTable_m[cardinality] = irredPolTable;
                  }
                  return ( irredPolHashTable_m[cardinality] );
            }

            typename FiniteFieldSearch::IrreduciblePolTableType & getIrredPolTableRef(uint cardinality)
            {
                  return *(getIrredPolTablePtr(cardinality));
            }

         
            
            /// todo: getIrredTester() kann parametrisiert werden werden!
            /// returns a pair: irreducibleTest itself and a value if the irreducibleTest is threadsafe.
            std::pair<typename FiniteFieldSearch::IrredTesterFktType, bool >    getIrredTester()  
            {
                
                //return PseudoIrreducibleExpensiveTest< IrreduciblePolTableType  >::isIrreducible;
                //return GAPIrreducibleTest< IrreduciblePolTableType  >::isIrreducible;
                return std::pair<typename FiniteFieldSearch::IrredTesterFktType,bool>(FLINTIrreducibleTest< IrreduciblePolTableType  >::isIrreducible,
                                                                                        FLINTIrreducibleTest< IrreduciblePolTableType  >::threadsafe)  ;

            }

            typename TPolRingTypePar::ElementType   convertPolRepToRingElem( const typename  HMSProblem::PolynomRepType & polRep) const;

        public:

            void normalizeInPlace( FiniteFieldSearch<TPolRingTypePar>::PolSetType &  polSet)    const;

            void removeConstantFactorsInPlace( FiniteFieldSearch<TPolRingTypePar>::PolSetType &  polSet)    const;

            ///@note formula \$ \frac{1}{degree}\sum_{d|degree}{\nue(degree/d)*q^d}  \$ where q is the cardinality of a finite field, is given by Gauss, see
            /// "Untersuchungen über höhere Arithmetik", second edition, reprinted, Chelsea publishing company, New York 1981 
            //CounterType  getIrredCount(int degree)   const;
            
            FiniteFieldSearch(const HMSProblem & hms, const SearchOptions& so,  const TPolRingTypePar & );


            std::list< PolynomialFactorBluePrint>  createPolFactorConstructionRules(std::vector<int > partition, uint exponent, uint destpolynomial);

            // extract a PolynomialFactorBluePrint to which is is possible to apply the Normalization rule.
            PolynomialFactorBluePrint* extractMatchingRule(std::list< PolynomialFactorBluePrint> &polFactorConstructRules, NormalizationRule rule);

            
            bool processNormalizationRules( PolSetBlueprintType & polsetblueprint, 
                                            std::list< PolynomialFactorBluePrint> &polFactorConstructRules, 
                                            std::vector<int>  & degreeOneRootList,
                                            ShapeList & shapeList,
                                            FiniteFieldSearch<TPolRingTypePar>::PolSetType &  polSet);

            void last_search_level( const ShapeList & shapeList, PolSetType & polSet );

            void renormalize( const ShapeList & shapeList, PolSetType & polSet );
    
            void first_third_search_level(DegSortedConstructionRuleTableType & polFactorConstructRules,  
                                    int maxDegree,
                                    const std::vector<int> & degreeOneRootList, 
                                    const ShapeList & shapeList,
                                    PolSetType & polSet);

            void  third_search_level(const std::vector< std::pair<int, const BPVecTYPE * > > & polFactorConstructRules,  
                                    const std::vector<int> & degreeOneRootList, 
                                    const ShapeList & shapeList,
                                    std::vector< std::pair<int, const BPVecTYPE * > >::const_reverse_iterator    vecIterator,
                                    PolSetType & polSet,
                                    mpz_t tmpCounter ,   mpz_t retCounter );

            void second_search_level(std::list< PolynomialFactorBluePrint> polFactorConstructRules);

            void first_search_level( std::list< PolynomialFactorBluePrint>     polfactorBlueprintList,  
                                      std::list< PolynomialFactorBluePrint>   polFactorConstructRules);
           
            void printStructureMap()    const
            {
                printStructureMap(fullStructureHashmap_m);  
            }

            bool computeScalingFactors(const PolSetType &  polSet,   std::vector< typename TPolRingTypePar::CoeffRingType::ElementType > & scalingRel);

            void printStructureMap(StructureHashMapType map)    const
            {

                size_t exampleBunchCount=0;
                StructureHashMapType::const_iterator it;
                for ( it=map.begin(); it!=map.end(); it++)
                {
                    exampleBunchCount += (*it).second;
                    std::cout << (*it).first << "=>" << (*it).second << ", "  << std::endl;
                    
                }
                #ifdef VERBOSE
                std::cerr << "example count " << exampleBunchCount << std::endl;
                #endif
            }
            void run()   ;

            void  printCount()   const   
            {                 
                    #ifdef VERBOSE
                    //    std::cerr << "# count:= " << counter_m << ";" << std::endl;   
                    #endif
                    char*  str= mpz_get_str( NULL, 10, counter_m);
                
                    std::string counterStr  (str);
                    
                    std::cout << counterStr << ";" << std::endl;   
                    delete[] str;
            };

            void  printMpz(mpz_t num)   const   
            {                 
                    #ifdef VERBOSE
                    //    std::cerr << "# count:= " << counter_m << ";" << std::endl;   
                    #endif
                    char*  str= mpz_get_str( NULL, 10, num);
                
                    std::string counterStr  (str);
                    
                    std::cout << counterStr << ";" << std::endl;   
                    delete[] str;
            };

            int64_t getCounter()    const   
            {
                int64_t result = mpz_get_ui(counter_m);   
                mpz_t tmp;
                mpz_init( tmp );
                mpz_set_ui( tmp,result );
                
                if ( mpz_cmp(tmp,counter_m)==0 )
                    return result;
                else
                {
                    std::cerr << " converting mpz_t failed ... ";
                    assert(false);
                }
            }

            virtual   ~FiniteFieldSearch();


    };

    template <typename TPolRingTypePar>
    typename FiniteFieldSearch<TPolRingTypePar>::PolRingHashTableType   FiniteFieldSearch<TPolRingTypePar>::polRingHashTable_m ;

    template <typename TPolRingTypePar>
    typename FiniteFieldSearch<TPolRingTypePar>::IrredPolHashTablesType     FiniteFieldSearch<TPolRingTypePar>::irredPolHashTable_m ;

 
    /// todo: idea: - first compute all possible factorDegree combinations, then sort them by max factorDegree and  
    /// then start checking with the smallest factorDegree ; also compute irreducible polynomial lists just before they are required
    /// 
    class HurwitzMapFinder
    {

       
    
        public:
            HurwitzMapFinder() {} ;
            
            template <class TPolRingTypePar>
            void finiteFieldSearch(const HMSProblem & hurwitzMapSearchProblem, const SearchOptions& searchOptions, const TPolRingTypePar & polRing)
            {
                assert( hurwitzMapSearchProblem.getLowerCharacBound() <= polRing.getCoeffRing().getCardinality() );

                assert( hurwitzMapSearchProblem.getNormalizationRuleList().getNormalizationRuleListAsVector().size()>0 );

                FiniteFieldSearch< TPolRingTypePar > ffs(hurwitzMapSearchProblem, searchOptions, polRing);
                assert( hurwitzMapSearchProblem.getNormalizationRuleList().getNormalizationRuleListAsVector().size()>0 );
                ffs.run();
                
                if (searchOptions.logStructure() )
                  ffs.printStructureMap();

                if (searchOptions.dryRun() )
                    ffs.printCount();
            }

            static void test()
            {    
               // //std::vector< Shape::ScalarType >    partition = { 4,3,2,2,2 };                      
               // Shape shape = { 4,3,2,2,2 }; //requires initializer lists

                //std::vector< Shape::ScalarType >    partition = { 4,3,2,2,2 };

                int ar[]={ 4,3,2,2,2 };
                const int TotalItems = sizeof(ar)/sizeof(ar[0]);
                std::vector< Shape::ScalarType >    partition(ar, ar+TotalItems);
                Shape    shape(partition);
                
                std::vector< Shape >  preShapeList(3,shape);
                
                ShapeList   shapeList(preShapeList);

                /* initializing shapelist using initializer lists
                //ok:
                 shapeList = std::vector<Shape>{ { 4,3,2,2,2 }, { 4,3,2,2,2 }, { 4,3,2,2,2 } } ;              
                //ok:
                 shapeList = ShapeList( { { 4,3,2,2,2 }, { 4,3,2,2,2 }, { 4,3,2,2,2 } } );

               shapeList = ShapeList(  { 4,3,2,2,2 }, { 4,3,2,2,2 }, { 4,3,2,2,2 }  );

                // ok:
                shapeList  = { { 4,3,2,2,2 }, { 4,3,2,2,2 }, { 4,3,2,2,2 } } ;

                */
                //ShapeList   shapeList  = { shape, shape, shape } ;

                HMSProblem  hurwitzMapSearchProblem(shapeList);

                assert(hurwitzMapSearchProblem.getNormalizationRuleList().getNormalizationRuleListAsVector().size()>0);
    
                bool dryRun,logStructure,strictNormalization;

                const SearchOptions   searchOptions=SearchOptions( dryRun=false, logStructure=true, strictNormalization=false);
                //const SearchOptions   searchOptions=SearchOptions( dryRun=false, logStructure=false, strictNormalization=false);
                HurwitzMapFinder    hmf;
                int characteristic = 7;
                
               const   TPolRingType::RingType* field = new    TPolRingType::RingType(characteristic,0);
                const   TPolRingType * ring = new TPolRingType(*field);

                FiniteFieldSearch<TPolRingType> ffs= FiniteFieldSearch<TPolRingType>(hurwitzMapSearchProblem, searchOptions, *ring);

                assert( getIrredCount(1, field->getCardinality() ) == 7 );
                std::cerr  << "ffs.getIrredCount(2)" <<  getIrredCount(2, field->getCardinality()) << std::endl;
                assert( getIrredCount(2, field->getCardinality()) == 21 );
                assert( getIrredCount(3, field->getCardinality()) == 112 );

                ffs.run();

                assert( ffs.getCounter() == 1280160 );

             
                
                hmf.finiteFieldSearch(hurwitzMapSearchProblem, searchOptions, *ring);
                //hmf.finiteFieldSearch(hurwitzMapSearchProblem, searchOptions, characteristic=7 );
                     return;
            }

            static void solve43222(int charac)
            {
                //std::vector< Shape::ScalarType >    partition = { 4,3,2,2,2 }; requires initializer lists   c++0x
                // std::vector< Shape::ScalarType >    partition = { 3,3,3,2,2 }; 


                int ar[]={ 4,3,2,2,2 };
                const int TotalItems = sizeof(ar)/sizeof(ar[0]);
                std::vector< Shape::ScalarType >    partition(ar, ar+TotalItems);
           
                Shape shape(partition);
                
                std::vector< Shape >  preShapeList(3,partition);
                
                ShapeList   shapeList(preShapeList);
                HMSProblem  hurwitzMapSearchProblem(shapeList);

                assert(hurwitzMapSearchProblem.getNormalizationRuleList().getNormalizationRuleListAsVector().size()>0);
    
                bool dryRun,logStructure,strictNormalization;

                const SearchOptions   searchOptions=SearchOptions( dryRun=false, 
                                                                    logStructure=true, 
                                                                    strictNormalization=false, 
                                                                     RationalMapSearch::OutputMode::M2Output   
                                                              
                                                                );
                //const SearchOptions   searchOptions=SearchOptions( dryRun=false, logStructure=false, strictNormalization=false);
                HurwitzMapFinder    hmf;
                int characteristic = charac;
                  
              const  TPolRingType::RingType* field = new   TPolRingType::RingType(characteristic,0);
                const   TPolRingType * ring = new TPolRingType(*field);

                hmf.finiteFieldSearch(hurwitzMapSearchProblem, searchOptions, * ring );
                //hmf.finiteFieldSearch(hurwitzMapSearchProblem, searchOptions, characteristic=7 );
                return;
            }
 

            //------------------------ 
            //  1 11  1     9 1       13  3  4 3 2 2 2 0 4 3 2 2 2 0 4 3 2 2 2 0 

            //------------------------  { (11 7 5 4 3 2^3),  (3^12), (2^18 )}

            //  3  5  1    1 3 1    0 2    36   3  11 7 5 4 3 2 2 2  0   3 3 3 3 3 3 3 3 3 3 3 3 0  2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0  # 488569597440   # noch machbar
            //  3  7  1    1 4 1    0 3    36   3  11 7 5 4 3 2 2 2  0   3 3 3 3 3 3 3 3 3 3 3 3 0  2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0  # 848877404860800 # schon zu schwer: faktor 100 fehlt         
            

            //------------------------ { (1^7  29), (3^12), (2^18 )}
            //  3  3  1    1 1 1    0 2    36   3  1 1 1 1 1 1 1 29 0   3 3 3 3 3 3 3 3 3 3 3 3 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0   # 492568168
            //  3  5  1    1 3 1    0 2    36   3  1 1 1 1 1 1 1 29 0   3 3 3 3 3 3 3 3 3 3 3 3 0  2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0  # 20140838436064 # gerade noch machbar
            //  3  5  1    1 3 1    0 2    36   3  1 1 1 1 1 1 1 29 0   2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0  3 3 3 3 3 3 3 3 3 3 3 3 0  # 319533755920577952 
            //  3  7  1    1 4 1    0 3    36   3  1 1 1 1 1 1 1 29 0   3 3 3 3 3 3 3 3 3 3 3 3 0  2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0  # 18925417228717488      
            //  3  11  1    1 9 1    0 2    36   3  1 1 1 1 1 1 1 29 0   3 3 3 3 3 3 3 3 3 3 3 3 0  2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 # 170504511106086922360

            //--------------------- { (1^5 21),  (3^8 2), (2^13) }
            // 3  3  1    1 1 1    0 2  26   3    1 1 1 1 1 21 0 3 3 3 3 3 3 3 3 2 0 2 2 2 2 2 2 2 2 2 2 2 2 2 0      # 566256
            //  3  5  1    1 3 1    0 2   26   3    1 1 1 1 1 21 0 3 3 3 3 3 3 3 3 2 0 2 2 2 2 2 2 2 2 2 2 2 2 2 0    # 1243688928
            //  3  7  1      1 4 1    0 3   26   3    1 1 1 1 1 21 0 3 3 3 3 3 3 3 3 2 0 2 2 2 2 2 2 2 2 2 2 2 2 2 0  # 168884296848
            //   3  11  1    1 9 1    0 2  26   3    1 1 1 1 1 21 0 3 3 3 3 3 3 3 3 2 0 2 2 2 2 2 2 2 2 2 2 2 2 2 0   # 111017632038720 # geschwindigkeitsfaktor 10 fehlt.
  

            //---------------------  { (7 5 4 3^2 2) ,( 3^8), (2^12) }

            //  3  3  1    1 1 1    0 2    24   3  7 5 4 3 3 2  0   3 3 3 3 3 3 3 3  0 2 2 2 2 2 2 2 2 2 2 2 2 0  # count : 22560
            //  3  5  1    1 3 1    0 2    24   3  7 5 4 3 3 2  0   3 3 3 3 3 3 3 3  0 2 2 2 2 2 2 2 2 2 2 2 2 0      #  93201120
            //  3  5  1    1 3 1    0 2    24   3  3 3 3 3 3 3 3 3  0    2 2 2 2 2 2 2 2 2 2 2 2 0  7 5 4 3 3 2 0     # 95931929420224 and also characteristic too small, needs slow multiplicity check
            //  3  7  1    1 4 1    0 3    24   3  7 5 4 3 3 2  0   3 3 3 3 3 3 3 3  0 2 2 2 2 2 2 2 2 2 2 2 2 0   # 16034945760
            //  3  7  1    1 4 1    0 3    24   3  3 3 3 3 3 3 3 3  0 2 2 2 2 2 2 2 2 2 2 2 2 0  7 5 4 3 3 2 0    #  125280700756918104, and also characteristic too small, needs slow multiplicity check
            //  3  11  1    1 9 1    0 2   24   3  7 5 4 3 3 2  0   3 3 3 3 3 3 3 3  0 2 2 2 2 2 2 2 2 2 2 2 2 0   # 11004081580800
            //  3  11  1    1 9 1    0 2   24   3  3 3 3 3 3 3 3 3  0 2 2 2 2 2 2 2 2 2 2 2 2 0  7 5 4 3 3 2 0     # 1806941616406887141260  

            //---------------------- { (7 5 4 3 2^2 1) , (3^8), (2^12) }
            //  3  3  1    1 1 1    0 2    24   3  7 5 4 3 2  2 1  0   3 3 3 3 3 3 3 3  0 2 2 2 2 2 2 2 2 2 2 2 2 0 # 0
            //  1  5  1     3 1     24   3  7 5 4 3 2  2 1  0   3 3 3 3 3 3 3 3  0 2 2 2 2 2 2 2 2 2 2 2 2 0  # 147415680
            //  1  7  1     4 1     24   3  7 5 4 3 2  2 1  0   3 3 3 3 3 3 3 3  0 2 2 2 2 2 2 2 2 2 2 2 2 0  # 50705740800
            //  1  11  1    9 1     24   3  7 5 4 3 2  2 1  0   3 3 3 3 3 3 3 3  0 2 2 2 2 2 2 2 2 2 2 2 2 0 # 75157785014400

            //---------------
            //  1  7  1    1 4 1    0 3  3 4     2 1 0    2 1 0    2 1 0    2 1 0    1

            static void testSchuett()
            {
          
              /*  Shape shape1 ={ 7,5,4,3,3,2 };
                Shape shape2 = { 3,3,3,3,3,3,3,3 };
                Shape shape3 = { 2,2,2,2,2,2,2,2,2,2,2,2 } ;
            */

                int ar1[]={ 7,5,4,3,3,2 };
                const int TotalItems1 = sizeof(ar1)/sizeof(ar1[0]);
                std::vector< Shape::ScalarType >    shape1data(ar1, ar1+TotalItems1);

                int ar2[]={ 3,3,3,3,3,3,3,3 };
                const int TotalItems2 = sizeof(ar2)/sizeof(ar2[0]);
                std::vector< Shape::ScalarType >    shape2data(ar2, ar2+TotalItems2);

                int ar3[]={ 2,2,2,2,2,2,2,2,2,2,2,2 } ;
                const int TotalItems3 = sizeof(ar3)/sizeof(ar3[0]);
                std::vector< Shape::ScalarType >    shape3data(ar3, ar3+TotalItems3);

                 std::vector< Shape > preshapeList ;
                preshapeList.push_back(Shape(shape1data));
                preshapeList.push_back(Shape(shape2data));
                preshapeList.push_back(Shape(shape3data));

                
                ShapeList   shapeList (preshapeList) ;

                /* using initializer lists
                ShapeList   shapeList(  { 7,5,4,3,3,2 },  
                                        { 3,3,3,3,3,3,3,3 } , 
                                        { 2,2,2,2,2,2,2,2,2,2,2,2 } 
                                     );
                */


                HMSProblem  hurwitzMapSearchProblem(shapeList);

                assert(hurwitzMapSearchProblem.getNormalizationRuleList().getNormalizationRuleListAsVector().size()>0);
    
                bool dryRun,logStructure,strictNormalization;

                 const SearchOptions   searchOptions=SearchOptions( dryRun=true,
                                                                     logStructure=false, 
                                                                     strictNormalization=false, 
                                                                     RationalMapSearch::OutputMode::M2Output  // requires  nested enumerators 
                                                                   );
                // const SearchOptions   searchOptions=SearchOptions( dryRun=false, logStructure=false, strictNormalization=false);
                HurwitzMapFinder    hmf;          // 4,3,2,2,2 char =23 :  103347928992
                                                                 
                //int characteristic = 7;                         //   scharf: 130 Std., keine Beispiele...
                // int characteristic = 5;                             //      20616540, keine Beispiele
                //int characteristic = 11;                       //     8014067740800 (mit genauen irreduziblen polynomen, 17321min)  
                                            // unsers mith char =23      103347928992
                // char 13                                             87611161223520
                                                        // char 7         20931320340
                                              
                int characteristic = 13;
          
                const  TPolRingType::RingType* field = new   TPolRingType::RingType(characteristic,0);
                const   TPolRingType * ring = new TPolRingType(*field);
                hmf.finiteFieldSearch(hurwitzMapSearchProblem, searchOptions, *ring );
               
            }
        
    
    };


};

#include "HurwitzMapFinder.hpp"


/*
//deg 36 example char 7:
18925417228717488
char 5:
20140838436064
char 3:
492568168
char 11:
4483814442700957816

165781036554829200 
11004081580800 //7 5 4 3 3 2 0
111017632038720 // deg26 char 11
75157785014400 // 7543221 0
7562712967133163520 //deg36 char 13
3822308396606137632 //deg36 char 13 erstes beispiel.

3695170689490405488 // deg26 char 23
111017632038720 // deg26 char 11


*/

/* rfs liefert:
{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1} =>825330
{ 1, 1, 1, 1, 1, 1} =>11
{ 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1} =>58468
{ 2, 1, 1, 1, 1} =>1
{ 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1} =>7383
{ 2, 2, 2, 1, 1, 1, 1, 1, 1, 1} =>1254
{ 2, 2, 2, 2, 1, 1, 1, 1, 1} =>166
{ 2, 2, 2, 2, 2, 1, 1, 1} =>23
{ 2, 2, 2, 2, 2, 2, 1} =>1
{ 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1} =>7609
{ 3, 2, 1, 1, 1, 1, 1, 1, 1, 1} =>368
{ 3, 2, 2, 1, 1, 1, 1, 1, 1} =>58
{ 3, 2, 2, 2, 1, 1, 1, 1} =>10
{ 3, 2, 2, 2, 2, 1, 1} =>2
{ 3, 3, 1, 1, 1, 1, 1, 1, 1} =>127
{ 3, 3, 2, 1, 1, 1, 1, 1} =>13
{ 3, 3, 2, 2, 1, 1, 1} =>1
{ 3, 3, 3, 1, 1, 1, 1} =>7
{ 3, 3, 3, 3, 1} =>1
{ 4, 1, 1, 1, 1, 1, 1, 1, 1, 1} =>1080
{ 4, 2, 1, 1, 1, 1, 1, 1, 1} =>42
{ 4, 2, 2, 1, 1, 1, 1, 1} =>4
{ 4, 2, 2, 2, 1, 1, 1} =>2
{ 4, 3, 1, 1, 1, 1, 1, 1} =>4
{ 4, 4, 1, 1, 1, 1, 1} =>3
{ 5, 1, 1, 1, 1, 1, 1, 1, 1} =>164
{ 5, 2, 1, 1, 1, 1, 1, 1} =>10
{ 5, 2, 2, 1, 1, 1, 1} =>2
{ 6, 1, 1, 1, 1, 1, 1, 1} =>14
{ 6, 2, 1, 1, 1, 1, 1} =>2

*/
