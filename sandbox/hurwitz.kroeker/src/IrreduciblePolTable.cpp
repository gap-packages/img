

#include "IrreduciblePolTable.h"
/*

Folgender Plan:
stat isIrred nur pseudo is irred.
Dann für deg>3 ein iterator über (pseudo)irreduzible Polynome
beim Ersten durchgang (pseudo)irreducible test und eine  Liste aufbauen.
Es koennte aber auch sein, dass aufgrund der Speicherbeschaffenheit die Liste zu gross
ist, nicht in den Cache passt und ein test/ Pseudo-test schneller ist. Also zumindest beide Varianten ausprobieren.

Problem: je mehr faktoren vom gleichen irred-Grad in den ersten beiden Polynomen vorkommen können,
desto eher lohnt es sich, die echte irreduzible liste anzulegen.

- bereits verwendete Polynome aus den Listen nur für Grad 1 -2  streichen, ansonsten lohnt es sich nicht mehr
(es sei denn, es gibt zu viele 'false positives' )
Wenn man nicht alle verwendeten faktoren streicht, muss man bei einem Ergebniskandidaten 
auch gcd(A,B)=1 gcd(A,C)=1 und gcd(B,C )=1 und zusätzlich den Shape von (A) und (B).



Prototyp: erstelle in GAP eine Liste mit irreduziblen polynomen vom bestimmten Grad.


-probiere das 43222-Beispiel mit dem alten C++-Code zu finden, wobei nur Grad1-Faktoren zugelassen werden.


*/

namespace RationalMapSearch
{

    
   

    const IntFactorTable::IntType   IntFactorTable::MaxInt = 400 ; 

    IntFactorTable::IntFactorTable()    
    {
        isPrime_m = new bool[ IntFactorTable::MaxInt ];
    
        isPrime_m[0]=false;
        isPrime_m[1]=false;
        for (uint n = 2; n < IntFactorTable::MaxInt; n++)
            isPrime_m[n] = true;
            for (uint n = 2; n < IntFactorTable::MaxInt; n++)
            if (isPrime_m[n])
            {
                uint c=0;
                for (uint m = 2; (c = m * n) < IntFactorTable::MaxInt; m++)
                    isPrime_m[c] = false;
            }

        IntFactorVecType oneFactors(1,1);
        intFactorTable_m.insert(std::pair< IntFactorTable::IntType, IntFactorVecType > ( 1, oneFactors ) );
    }


    IntFactorTable::IntFactorVecType  IntFactorTable::computeIntFactors(IntType integer)
    {
        assert( integer < IntFactorTable::MaxInt );
        IntFactorVecType     res;
    
        for (IntType n = 2; n < IntFactorTable::MaxInt; n++)
        {   
            if (isPrime_m[n])
            {
                IntFactorTable::IntType c = 0;
                for (IntType m = 1;   (c = m * n) <= integer; m++)
                {
                    c = m*n;
                    if (c==integer) 
                    {
                        res.push_back(n);
                        break;
                    }
                }
            }
                 
        }
        return res;
    }
  
  

   void IntFactorTable::test()
            {
                   IntFactorTable   ift;
                   IntFactorTable::IntFactorVecType factorList = ift.getFactors(10);
                   std::cerr << "factorList.size()" << factorList.size();
                   std::cerr << "factorList[0]" << factorList[0];
                   assert( factorList.size()==2 );

                   factorList = ift.getFactors(7);
                   assert( factorList.size()==1 );
                   std::cerr << std::endl;
                   std::cerr << "IntFactorTable test passed" <<  std::endl;
            }


    IntFactorTable::IntFactorVecType IntFactorTable::getFactors(IntType integer)
            {
                if ( intFactorTable_m.find(integer ) == intFactorTable_m.end() )
                {
                    addTableEntry(integer);
                }
                return intFactorTable_m[integer];
            }

}