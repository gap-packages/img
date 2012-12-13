#ifndef TYPEDEFS_H
#define TYPEDEFS_H


/** @file typedefs.h 
*
*   @brief Contains \b  ulong64 and \b  long64 typedefs for Windows and Linux environment
           and P_or_QPolynom-enum.
*
* @todo P_or_QPolynom definition does not belong here!
*/


#if defined(_MSC_VER) || defined(__BORLANDC__)

	typedef unsigned __int64 ulong64;

	typedef   signed __int64  long64;

#else   //Linux-Environment:

	typedef unsigned long long ulong64;

	typedef signed long long  long64;

	typedef short scalarType;

	typedef short scalarType;

	typedef short epsScalarType;
	
#endif



#endif
