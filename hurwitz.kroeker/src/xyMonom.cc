
#include "xyMonom.h"
#include "parseTools.h"



using namespace std;


/** @brief reads a number from input stream (an exponent)
 is only for the case when in a input stream '^' was read - as next an exponent number is expected and
*/
unsigned int xyExpPair::extractExplicitExponent(stringstream&  sstream)
{
	string res="";
	
	int num=0;
	int sign=getSign(sstream);

	if (sign<0)
	{
		throw "negative exponents not allowed fom monoms! " ;
	}

	sstream >> ws;
	sstream >>num;
	if (sstream.fail())
	{
		throw "no exponent given" ;
	}
	return num;
}





/** @brief liest die Exponenten der Variablen (x,y).
 Nur positive Exponenten erlaubt.<br>

 Erlaubt ist also so etwas wie 'x' , 'y' , 'xy' 'x^2y 'xy^3' 'x*y'
 'x^2*y^5' aber nicht 'x^-1'!
 wenn x oder y, dann entweder ^oder + erlaubt
 einschränken: wenn *, dann nur noch x oder y erlaubt.
*/
xyExpPair::xyExpPair(std::stringstream&  sstream)
{
	#ifdef DEBUG
		std::cerr << "create xyExpPair" << std::endl;
		std::cerr << "xyExpPair stream: =" << sstream.str() << std::endl;
	#endif

	assert (!sstream.fail());

	x_exp=0;
	y_exp=0;

	char a;	
	bool isXExp=false;
	bool isYExp=false;

	sstream >> ws;
	
	while (!sstream.eof() )
	{
		assert (!sstream.fail());
		a=sstream.peek();
		#ifdef DEBUG
			std::cerr << "a " << a;
		#endif
		if (a=='*')
		{
			isXExp=false;
			isYExp=false;
			sstream >>a;
			sstream >> ws;
			continue;
		}
		if (a=='x')
		{
			x_exp=1; // der x-Exponent ist implizit  mindestens 1 , wenn kein '^' folgt
			isXExp=true;
			isYExp=false;
			sstream >>a;
			sstream >> ws;
			continue;
		}
		else 
		if (a=='y')
		{

			y_exp=1; // der y-Exponent ist implizit mindestens 1, wenn kein '^' folgt
			isYExp=true;
			isXExp=false;
			sstream >>a;
			sstream >> ws;
			continue;
		}
		else
		if (a=='^')
		{
			sstream >>a;
			if (isXExp)	
			{
				x_exp=extractExplicitExponent(sstream);
				//assert(x_exp!=0);  - kann auch 0 sein!
			}
			else if (isYExp)
			{
				y_exp=extractExplicitExponent(sstream);
				//assert(y_exp!=0);
			}
			else throw "[xy]-polynom: unknown variable";
			sstream >> ws;
			continue;
		}
		else 
		{
			sstream >> ws;
			if (sstream.eof() || a=='-' || a =='+' )		
				return;
			else
			{
				throw "error during extracting (x,y) exponents";
			}
		}
	}
}



