#ifndef COMPILE_FUNCTIONS
#define COMPILE_FUNCTIONS

#include <stdio.h>

#include <string.h>
#include <string>

/** \file CompileFunctions.h
*
* @brief contains help functions which are evaluated at compile-time

* contains help functions which are evaluated at compile-time<br>
* Limitation: the computation is recursive and there are some compiler-depending recursion limits, usually 256.
*
* @todo Dateiname irreführend -> ändern!
*/

//--------------------- Compile-Hilfsfunktionen --------------------

// Vorbelegung von Operationstabellen(z.B: Additionstabellen)  für std::endliche Körper:
// Auch wenn die Template-Berechnungen zur Compilezeit etwas bringen:
// die Definition eines statischen const array mit BOOST erwies sich als zu schwierig,
// (wegen der Syntax mussten Praeprozessorwiederholungen eingesetzt werden,
// diese sind aber meistens auf 256 wiederholungen beschraekt;// siehe TemplateExperiments.h)
// Daher sollten die const Definitionen fuer die Operationsstabellen mit einem Hilfsprogramm erstellt werden


/** @brief  pow2<x>::value computes \f$ \mbox {\large  $ 2^x $ } \f$ during compile time<br>
    the computation is recursive and there are some compiler-depending recursion limits, usually 256.
* 
*  @note the computation is recursive and there are some compiler-depending recursion limits, usually 256.
* @ingroup helpFunctions
*
* @author Jakob Kröker
* 
* @todo is dangerous, because enum value is limited!
*
*/
template <int NN>
struct pow2		
{
	enum	{value		= 2*pow2<NN - 1>::value    };
	enum	{valueMinusOne  = 2*pow2<NN - 1>::value -1};
};


/** @brief  pow2<0>::value  represents \f$ \mbox {\large  $ 2^0=1 $ } \f$ during compile time
*
* @ingroup helpFunctions
*
* @author Jakob Kröker
*
*/
template <>
struct pow2<0>	
{

	enum	{	value=1	};
};


/** @brief needbits<x>::value computes number of required bits to represent unsigned integer <b>x</b> during compile time
*
* @ingroup helpFunctions
*
* @author Jakob Kröker
*
*/
template <int NUM>
struct needbits 
{
   	enum	{	value 		= needbits<NUM/2>::value +1		};
	enum	{	valueplusone 	= needbits<NUM/2>::valueplusone +1	};
   	enum	{	doubledvalue 	= needbits<NUM/2>::doubledvalue +2	};
   
};


/** @brief needbits<1>::value computes number of required bits to represent integer <b>1</b> during compile time
*
* @ingroup helpFunctions
*
* @author Jakob Kröker
**/
template <>
struct needbits<1>
{ 
	enum	{	value =1	};
	enum	{	valueplusone = 2	};
	enum	{	doubledvalue = 2	};
 	
};

template <>
struct needbits<0>
{ 
	enum	{	value =1	};
	enum	{	valueplusone = 2	};
	enum	{	doubledvalue = 2	};
 	
};

/**  @brief nextpow2num<x>::value computes (at compile time) a power of 2 -value 
          such that  \f$ \mbox { \large $ 2^x < value=2^a $ } \f$ is valid,
         where x is the minimal number of bits to represent  the unsigned integer CHAR
* @ingroup helpFunctions
*
* @author Jakob Kröker
*
*/
template <int CHAR>
struct nextpow2num
{
	enum {value = pow2<  needbits<CHAR>::value  >::value };
};



// ----------Ende Compile-Hilfsfunktionen --------------------


#endif
