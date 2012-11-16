

#include "hmfTypedefs.h"
#include "HurwitzMapFinder.h"

//#include "FactorPolynomialWrapper.h"

#include "DebugLogger.h"


template<class TPolynomial>
std::ostream& printUnivarPolynomial(std::ostream& os, const TPolynomial &pol, std::string varName = std::string("a") )
{
    bool first=true;
    for (size_t pos=0; pos<pol.size(); pos++)
    {

        if (not first) 
            if ( pol[ pol.size()-1-pos ]>0)
                os << "+" ;

        if ( pos==  pol.size()-1 )
        {
            if (pol[pol.size()-1-pos] !=0 )
            {
                os << pol[ pol.size()-1-pos];        
                first=false;
            }
        }
        else
        {
            if ( pol[ pol.size()-1-pos ] !=0 )
            {
                first=false;

                if (pol[ pol.size()-1-pos ] != 1 )
                    os << pol[ pol.size()-1-pos ] ;
                os << varName;
                if ( pol.size()-1-pos >1)
                    os  <<"^" << pol.size()-1-pos;
            }
        }
    }
    os << std::endl;
    return os;
}


/// @todo: improve: indroduce a parameter object, which is initialized by reading an input stream and detects ill-formed input.
void simpleCommandLineInterface(int argc, char* argv[])
{


    if (argc>1)
    {        
        std::string strFirstCmd(argv[1]);
        if ( strFirstCmd.compare( std::string("--help") )==0 || strFirstCmd.compare( std::string("-h") )==0  )
        {
            std::cerr << "# input format:    "  << std::endl;
            std::cerr << "#  flags {finite field prime}  {prime field extension degree}  {finite field mod polynomial}  {shapes: polynomial degree}  {branch value count}  {first shape} .. {last shape}   {prime field-reduced branch value approx v_4} .. v_{branch value count} "  << std::endl;
            std::cerr << "# a shape is a list of factor exponents  separated by spaces  and terminated by '0'  "  << std::endl;
            std::cerr << "# flags :   "  << std::endl;
            std::cerr << "#        | 2^0 (first bit)  : debug (y/no) "  << std::endl;
            std::cerr << "#        | 2^1 (second bit) : just count search space size"  << std::endl;
            std::cerr << "#        | 2^2 (third bit)  : strict normalization; factors with multiplicities equal to first degree partition entries are normalized "  << std::endl;
        
            std::cerr << "# example for three degree 13  (43222)-Shapes  in characteristic 11 : " <<  std::endl;
            std::cerr << "# 5  11  1   9 1   13   3   3 4 2 2 2 0    3 4 2 2 2 0     4 3 2 2 2 0 " <<  std::endl;
   
        }
      }
    bool dryRun = false;
    bool strictNormalization=false;
    
    int prime = 2;
    int extensionDegree = 1;
    int flags = 0;
    int debugLevel = 0;
    unsigned int modPolDegree = 0;
    std::vector<int>    modPolCoeffs;
    int currCoeff = 0;


    int branchValueCount = 0;
    int polDegree = 0; // todo: variablennamen verbessern

    std::cin >> std::ws >> flags;
   
    
    if (flags & 1  )
        debugLevel=1;

    dryRun =  (flags & 2 );

    strictNormalization =  (flags & 4 );


    DebugLogger::setLevel( debugLevel );
    
    DebugLogger::logStream() << "#I flags :  " << flags << std::endl;

    std::cin >> std::ws >> prime;

    std::cin >> std::ws >> extensionDegree;

    DebugLogger::logStream() << "#I debugLevel :  " << debugLevel << std::endl;

    DebugLogger::logStream() << "#I dryRun :  " << dryRun << std::endl;

    DebugLogger::logStream() << "#I prime :  " << prime << std::endl;

    DebugLogger::logStream() << "#I extension degree  :  " << extensionDegree << std::endl;

   

    if ( extensionDegree != 1 )
    {
           DebugLogger::logStream() << "#I galois fields  are currently not suppurted" << std::endl;
        assert(false);
    }
 
    
    //std::cin >> std::ws >> modPolDegree ;

    //DebugLogger::logStream() << "#I modPolDegree :  " << modPolDegree << std::endl;
    modPolDegree=extensionDegree;

  
    for (int currDegree=0; currDegree <= modPolDegree; currDegree++)
    {
              std::cin >> std::ws >> currCoeff;
              assert( std::abs( currCoeff )<prime );
              if (currDegree==modPolDegree)
                assert(currCoeff != 0 );
              modPolCoeffs.push_back(currCoeff);
    }
    
    DebugLogger::logStream() << "#I modPolynomial ";
    printUnivarPolynomial( DebugLogger::logStream(), modPolCoeffs);
  
 

    std::cin >> std::ws >> polDegree;
    DebugLogger::logStream() << "#I polynomial degree : " << polDegree << std::endl;

   

    std::cin >> std::ws >> branchValueCount;

    DebugLogger::logStream() << "#I branchValueCount = " << branchValueCount << std::endl;

    std::vector< RationalMapSearch::Shape >  shapeList;

     std::vector< int  >  normalizationExponents;

    for (int currShapePos=0; currShapePos< branchValueCount; currShapePos++)
    {
        std::vector<int>    partition;
        int exponent=0;
        int expSum=0;
        do
        {
              std::cin >> std::ws >> exponent;
              if (exponent!=0)
                partition.push_back(exponent);
              expSum += exponent;
            assert( expSum<=polDegree );
        }
        while (exponent !=0 );
        assert( expSum==polDegree );

        if ( currShapePos<3)
        {
            assert( partition.size()>0 );
            if (strictNormalization)
                normalizationExponents.push_back( partition[0] );
            else
                normalizationExponents.push_back( RationalMapSearch::NormalizationRule::dontcare );
        }
        
         RationalMapSearch::Shape shape =  RationalMapSearch::Shape(partition);
         DebugLogger::logStream() << "#I shape [" << currShapePos << "] :"  << shape << std::endl;

        shapeList.push_back( shape );
    }
    DebugLogger::logStream() << "#I shapeList constructed  " << std::endl;

    if (strictNormalization)
    {
        DebugLogger::logStream() << "#I strictNormalization  " << std::endl;
         DebugLogger::logStream() << "#I normalize factors with multiplicities ";
        for (size_t pos=0;pos< normalizationExponents.size();pos++)
        {
            if (pos>0) 
                 DebugLogger::logStream() << ",";
            DebugLogger::logStream() << normalizationExponents[pos] << " ";
        }
        DebugLogger::logStream()  << std::endl;
    }
     

 
    std::vector <int > reducedBranchValues;
 
    for (int currShapePos=3; currShapePos< branchValueCount; currShapePos++)
    {
            RationalMapSearch::HMSProblem::PolynomRepType     minimalPolynomial;
            int reducedBranchValue = 1;
            std::cin >> std::ws >> reducedBranchValue;
            reducedBranchValues.push_back(reducedBranchValue);
    }
    
  //////////////////////////
  


    bool  logStructure;

    const  RationalMapSearch::SearchOptions   searchOptions= 
                                            RationalMapSearch::SearchOptions( dryRun, 
                                                                              logStructure=false,     
                                                                              strictNormalization, 
                                                                              RationalMapSearch::OutputMode::GAPOutput );
    
    RationalMapSearch::HurwitzMapFinder    hmf;         
                                                        
    // gehoert 'prime' zu SearchOptions oder nicht? Eigentlich schon, aber dann auch der modPol, generator und extensionDegree.

    if (extensionDegree==1)
    {
        assert( modPolCoeffs.size()==2 );

        
        const  TPolRingType::RingType* field = new   TPolRingType::RingType(prime,0, prime - modPolCoeffs[0] % prime );
        const   TPolRingType * ring = new TPolRingType(*field);

       std::vector <    RationalMapSearch::HMSProblem::PolynomRepType  > minimalPolynomials;

        for (int currReducedBranchValuePos=0; currReducedBranchValuePos< reducedBranchValues.size() ; currReducedBranchValuePos++)
            {
                    RationalMapSearch::HMSProblem::PolynomRepType     minimalPolynomial;
                    assert(   reducedBranchValues[currReducedBranchValuePos] < field->getCardinality() &&  reducedBranchValues[currReducedBranchValuePos] >0 );
                    // set minimal polynomial to 'x-reducedBranchPoint';
                    DebugLogger::logStream() << "#I reduced branch value no " << currReducedBranchValuePos+4 << " : " << reducedBranchValues[currReducedBranchValuePos] << std::endl;
        
                    assert( extensionDegree == 1 ); //for extension field this part differs.
                    int fieldElement=field->generatorExponentToElem( reducedBranchValues[currReducedBranchValuePos] );
                    assert( std::abs( fieldElement ) < field->getCardinality() );
                    DebugLogger::logStream() << "# reduced branch value  fieldelem  = " << fieldElement << std::endl;
                    minimalPolynomial.push_back( - fieldElement );
                    minimalPolynomial.push_back(1);
                    minimalPolynomials.push_back( minimalPolynomial );
            }

        assert (normalizationExponents.size()==3 );
        RationalMapSearch::NormalizationRule infinity = RationalMapSearch::NormalizationRule(0, normalizationExponents[0], RationalMapSearch::NormalizationValue::infinity );
        RationalMapSearch::NormalizationRule zero = RationalMapSearch::NormalizationRule(1, normalizationExponents[1], RationalMapSearch::NormalizationValue::zero);
        RationalMapSearch::NormalizationRule one = RationalMapSearch::NormalizationRule(2, normalizationExponents[2], RationalMapSearch::NormalizationValue::one);

        //std::vector <RationalMapSearch::NormalizationRule> preRuleList =  {infinity,zero,one};
        std::vector <RationalMapSearch::NormalizationRule> preRuleList ;
        preRuleList.push_back(infinity);preRuleList.push_back(zero);preRuleList.push_back(one);
        RationalMapSearch::NormalizationRuleList nrl= RationalMapSearch::NormalizationRuleList(preRuleList, strictNormalization);

        RationalMapSearch::HMSProblem  hurwitzMapSearchProblem( shapeList,nrl, minimalPolynomials );

        hmf.finiteFieldSearch(hurwitzMapSearchProblem, searchOptions, *ring );  
    }
    else
        assert(false);


    return;
}


/// test 5 11 1   9 1  3 4 2 1 0  2 1 0 2 1 0 2 1 0 2 1 0

int main(int argc, char* argv[])
{
    simpleCommandLineInterface(argc, argv);

    return 0;
}
