
#pragma once


#include <string>
#include <assert.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <string.h>


/** \file parseTools.h
 *   @brief  input file and input stream parsing utils
 * @todo die Datei  ist bestimmt noch voller Fehler... (fail-test, potenzielle Segfaults, string escaping, Kommentare in der Eingabedatei, etc.)
*/




//--------------------------Parser Kram ----------------------------------------------

/** @brief  jumps over next whitespaces series in the string 'str', starting at position 'stelle'. *
 The position is then set to the next first position with a non-whitespace char 
*/
inline void 		eatWS(char* str, unsigned int & stelle);


/** @brief 'eats' leading whitespaces in the string str, and returns the result string
 */
inline std::string 	eatWS(const std::string &str);


/** @brief extracts given char from a stream, throws an error, if failed
 */
template <class _istream>
char 		extractChar(char schar, _istream& sstream);

/*
template <class _istream>
char 		extractChar(char schar, std::string  str)
{
	//std::stringstream strstream(str);
	//extractChar(schar,strstream );
	//str=strstream.str();
	return schar;
}
*/

/** @brief extracts given char from a stream, throws an error, if failed
 */
template <class _istream>
char 		extractChar(_istream& sstream, char schar);


/** @brief extracts a '-' or '+' from sstream if possible, and returns -1 if '-' , otherwise 1.
 */
template <class _istream>
int 		getSign(_istream& sstream);


/** @brief if leading stream entries is a Macaulay-comment sequence, the function strips all chars until next line starts .

if leading stream entries are  \verbatim '--'  (Macaulay-comment) \endverbatim, strips all chars until next line starts (after \verbatim '\n' or '\r' \endverbatim).
 */
template <class _istream>
std::string 	stripComment(_istream& data, bool & isComment);


/** @brief jumps over Macaulay-comments  \verbatim '--' \endverbatim ,until a non-comment char occurs in the input stream
 */
template <class _istream>
void 		stripComments(_istream& data);


/** @brief try to read a comment from input stream. Comment starts with two  \verbatim '--' \endverbatim. Whether  comment reading was successful, return value 
 contains the extracted comment string, otherwise the return value is empty. 
 */
template <class _istream>
std::string 	readCurrentComment(_istream& data);


/** @brief seeks tho given braceLevel in the input stream if possible.
      Needs information about currentBraceLevel, which is also kept up to date */
template <class _istream>
void 		seekToBraceLevel(_istream & inputData, int &currentLevel, int braceLevel);


/** @brief extracts data until ',' of '}' at current brace level occur. Comments ( \verbatim '--' \endverbatim) are stripped
 */
inline std::string extractNextData(std::stringstream &Data);


/** @brief  extracts data starting with '{'  until of '}' at current brace level occurs.
 */
inline std::string extractNextBracedData(std::stringstream &Data);


/** @brief  Reads Macaulay parameter name in current line. Assumption: the input stream points 
            to first char of a parameter name (whitespaces are ignored). Parameter name ends with '='.*/
template <class _istream>
std::string 	readParamName(_istream & inputData);


/** @brief  Reads macaulay parameter value in current line. Assumption: the input stream points 
            to first char of a parameter name (whitespaces are ignored). Parameter name ends with '='. */
template <class _istream>
std::string 	readParamValue(_istream & inputData, std::string name);


/** @brief  count subgroups in braces,  example: count result for "{1,2,{..}}"-stream  is 1. 
	   Assumes, that comments are already stripped. Inside string are currently not handled.
*/
inline int 	countSubGroups(std::stringstream &polynomStream);


/** @brief  count groups in braces,  example: count result for "{1,2,{..}},{...},{...},{...}"-stream  is 4.
            Assumes, that comments are already stripped. Inside string are currently not handled.
*/
inline int 	countGroups(std::stringstream &polynomStream);


/** @brief  count elements in braces,  example: count result for "{1,2,{4,5}}"-stream  is 3. 
           Assumes, that comments are already stripped. Inside string are currently not handled.
*/

inline int 	countElements(std::stringstream &str);

//-------------------------------------------------------------------------------------------------------------

/** @brief reads scalar from string <b>strNum</b> with following format: \code <num>+<num>*e^1++<num>*e^2+... \endcode
 *  and converts result before returning in a element of Ring <b>ring1</b> 
 *  with maximal possible epsPrecision<=<b>_epsPrecision</b> for template class <b>RingType::ElementType</b>
 * 
 * from <b>kk</b> only one thing is important: the convert-function. 
 @todo evtl. weitere Verbesserung: ein Element kann sich selbst aus einem String einlesen.
Dabei ist aber zu beachten, dass hier eventuell  epsPrecision eingeschraenkt werden soll (Parameter _epsPrecision)
Eine weitere Moeglichkeit besteht darin, den Ring in dieser Funktion zu erzeugen. 
*/
template <class RingType>
typename RingType::ElementType	 parseNumber(RingType* ring1, std::string strNum, int _epsPrecision )
{
	std::stringstream ssNum;
	ssNum << strNum;

	// erzeute ein Element mit vorgegebener epsPrecsion.
	typename RingType::ElementType erg(_epsPrecision, "dummy");
	//=createNum<TNum>(_epsPrecision);
	
	assert(erg.getEpsPrecision()>=_epsPrecision);

	
	//    Algoritmus: solange nicht ganzen string ausgelesen:
	//	      solange weiterer Summand{
	//           -lese {Vorzeichen,SummandKoeff}
	//           -lese epspotenz(Summand)
	//           -addiere Summand zu erg, falls epspotenz<=_epsPrecision
	//	}
	//      -Wandle das Ganze in trage Ziffer in erg ein 

	while (!ssNum.fail() && !ssNum.eof())
	{
		int number=1;

		int sign=getSign(ssNum);
 		if (ssNum.peek()!='e')
			ssNum >> number;
		number=number*sign;

		////end  read coefficient
		//// read epsPrecision
		int precision=0;

		// wenn es eine epsZahl ist, muss jetzt '*' folgen
		if (!ssNum.eof() && ( ssNum.peek()=='*' || ssNum.peek()=='e')) // soll * weggelassen werden koennen???
		{
			// assert nicht gut... Programm sollte wenigstens zu Ende Laufen, wenn ein Polynom
			// nicht eingelesen werden kann. Oder zu Beginn sollten testweise alle Polynome eingelesen werden.
			if (ssNum.peek()=='*')
				extractChar('*',ssNum);
			extractChar('e',ssNum);
			if (ssNum.peek()=='p')
			{
				extractChar('p',ssNum);
				extractChar('s',ssNum);
			}
			if ( !ssNum.eof() && ssNum.peek()=='^')
			{
				extractChar('^',ssNum);
				ssNum >> precision;
			}
			else
			{
				precision=1;  //epsExponent ist voraussichtlich 1
			}

		}
		if (precision<=_epsPrecision)
		{
			//int number2=erg[precision];// get+set or addValue
			erg.setValue(precision, ring1->getField()->ConvertScalar( number) );
			//erg[precision]+=number;
		
		}
	}
	ssNum >> std::ws;
	if (ssNum.fail() || !ssNum.eof())
	{
		std::cerr << "error during parseNumber, wrong format ?" << std::endl;
		throw   "error during parseNumber, wrong format ?";
	}
	return erg;
}


/// monomgroup: monoms with dame degree, format (implicit) {coeff*x^deg*y^0,coeff*x^(deg-1)*y^0}. in Reality{coeff,coeff,...}
inline int getMonomDegree(std::string monomGroup)
{
	int monomDegree=0;
	for (size_t strPos=0;strPos<monomGroup.length();strPos++)
	{
		if (monomGroup[strPos]==',')
			monomDegree++;
	}
	return monomDegree;
}

inline int getMaxMonomDegree( std::string polynomString)
{
	#ifdef DEBUG
	std::cerr << "getMaxMonomDegree" << std::endl;
	#endif

	int MaxMonomDegree=0;

	std::stringstream 	  polynomStream1;
	polynomStream1 << polynomString;
	polynomStream1 >> std::ws;
	extractChar('{',polynomStream1);

	#ifdef DEBUG
		std::cerr << "polynomStream1" << polynomStream1.str() << std::endl;
	#endif

	while (!polynomStream1.eof() && polynomStream1.peek()=='{')
	{
		std::string monomGroup1=extractNextBracedData(polynomStream1);
		polynomStream1 >>std::ws;
		if ( polynomStream1.peek()!='}' )		
			extractChar(',',polynomStream1);
		polynomStream1 >>std::ws;
		// get degree of monoms in this grout : is equal to number of elements decreased by 1.
		int currMonomDegree=getMonomDegree(monomGroup1);
		if (currMonomDegree>MaxMonomDegree)
			(MaxMonomDegree=currMonomDegree);
	}
	#ifdef DEBUG
	std::cerr << "getMaxMonomDegree : finished" << std::endl;
	#endif
	return MaxMonomDegree;
}





inline void eatWS(char* str,unsigned int & stelle)
{
	while (isspace(str[stelle])&&(stelle<strlen(str))) 
	{
		stelle++;
	}
}


inline std::string ltrim(const std::string &str)
{
	std::string result="";
	size_t pos=0;
	for ( ; pos < str.length();pos++) 
	{
		if (!isspace(str[pos]))
			break;
	}

	for (; pos < str.length();pos++) 
	{
		result=result+str[pos];
	}
	return result;
}

inline std::string rtrim(const std::string &str)
{
	std::string result="";
	int pos=str.length()-1;
	for ( ; pos > 0 ; pos--) 
	{
		if (!isspace(str[pos]))
			break;
	}

	for (int pos2=0; pos2 <= pos ;pos2++) 
	{
		result=result+str[pos2];
	}
	return result;
}


inline std::string trim(const std::string &str)
{
	std::string res=ltrim(str);
	return rtrim(res);
}

inline std::string eatWS(const std::string &str)
{
	std::string result="";
	for (size_t pos=0; pos < str.length();pos++) 
	{
		if (!isspace(str[pos]))
			result=result+str[pos];
	}
	return result;
}

//--------------------------Allgemeiner Parser Kram-----------------------------



template <class _istream>
char extractChar(char schar, _istream& sstream)
{
	char a;
	sstream >> std::ws;
	//#ifdef DEBUG
		std::string sschara;
		sschara+=schar;
		std::string sa;
		
	//#endif
	
	if (sstream >> a)
	{
		
		#ifdef DEBUG
			sa+=a;
			std::cerr << "extract char: dest '" << sschara << "' ;real '" << sa << "'" << std::endl;
		#endif
 		if ( a==schar)
		{
			return a;
		}
		else 
		{
			std::cerr <<"extract char: dest " << schar << " real "<< a << std::endl;
		}
	}
	else 
	{
		if (sstream.eof())
			std::cerr << "sstream.eof" <<  std::endl;
		std::cerr << "failed read '" << sschara << "' !!!" <<  std::endl;
	}

	std::string error=std::string("extract char ") + schar + std::string("failed");
	throw error;
}

template <class _istream>
char extractChar(_istream& sstream,char schar)
{
	return extractChar(schar,sstream);
}


template <class _istream>
int getSign(_istream& sstream)
{
	char a;
	sstream >> std::ws;
	if (!sstream.eof() )
	{
		a=sstream.peek();

		if (a=='-')
		{
			sstream >>a;
			return -1;
		}
		if (a=='+')
		{
			sstream >>a;
		}
	}
	return 1;
}




template <class _istream>
std::string stripComment(_istream& data, bool & isComment)
{
	char a;
	std::stringstream scomment;
	isComment=false;

	assert(! data.fail() );
	data >> std::ws;
	assert(! data.fail() );

	if ( !data.eof()  )
	{
		data >> a;
		assert(! data.fail() );
		if ( data.eof() )
		{
			data.clear();
			//cerr << "putback1 " << a << endl;
			data.putback(a);
			isComment=false;
			assert(! data.fail() );
		}
		else
		{
			assert(! data.eof() );
			char b = data.peek();
			if ( data.eof() )
			{
				data.clear();
  				data.putback(a);
				isComment=false;
			}
			else
			{
			//	cerr << "peek b :" << b << endl;
				assert(! data.fail() );
		
				if (a=='-' && !data.eof() && b=='-')
				{
					data >> a;
					isComment=true;
					assert(! data.fail() );
				}
				else
				{
					assert(! data.eof() );
					assert(! data.fail() );
			//		cerr << "putback2 " << a << endl;
					data.putback(a);
					isComment=false;
					assert(! data.fail() );
				}
			}
		}
	}
	assert(! data.fail() );
	while (!data.eof() && isComment)
	{
		data >> std::noskipws >> a;
	//	std::cerr <<  a;
		if ( a=='\n' || a =='\r')
		{
			break;
		}
		scomment << a;
	//	cerr << "while scomment.str()" << scomment.str() << endl;
	}
	if ( data.fail() )
	{
		std::cerr  << " Failed strip comment.";
		throw "Failed strip comment.";
	}
	//cerr << "scomment.str()" << scomment.str() << endl;
	return scomment.str();
}


template <class _istream>
void stripComments(_istream& data)
{
	bool wasComment=true;
	while (wasComment)
	{
		stripComment(data, wasComment);
	}
	return;
}


template <class _istream>
std::string readCurrentComment(_istream& data)
{	bool wasComment;
	return stripComment(data,wasComment);
}


template <class _istream>
void seekToBraceLevel(_istream & inputData, int &currentLevel, int braceLevel)
{
	char a;
	while (!inputData.eof() && currentLevel!=braceLevel)
	{
		inputData >> a;

		if (!inputData.eof() && a=='-'  && inputData.peek()=='-')
		{
			inputData.putback(a);
			
			readCurrentComment(inputData);
		}

		if (a=='{')	currentLevel++;
		
		if (a=='}')	currentLevel--;

		if (a==';')	throw "seekToBraceLevel: Unexpected End of Value";
		
		if (a=='=')	throw "seekToBraceLevel: Unexpected Assignment";

		if (currentLevel==braceLevel)
			break;
	}
	if (currentLevel!=braceLevel)
		throw "seekToBraceLevel:failed seek to brace level";
}

/// @TODO Kommentieren !
inline std::string extractNextData(std::stringstream &Data)
{
	char a;
	std::stringstream result;

	int braceLevel=0;

	while (true )
	{
		#ifdef DEBUG
			std::cerr << "stripComments:  " <<  result.str() << std::endl;
		#endif 
		stripComments(Data);
		#ifdef DEBUG
			std::cerr << "stripComments finished" <<  result.str() << std::endl;
		#endif 

		if (Data.eof())
			throw "Unexpected eof in data";
		if (  Data.peek()=='}' )
		{
			if  (braceLevel==0)
				break;
		}
		if   ( Data.peek()==',' &&   braceLevel==0) 
			break;
		
		Data   >> a;
		result << a ;
		if (a=='{')
			braceLevel++;
		if (a=='}')
			braceLevel--;
		
		#ifdef DEBUG
			std::cerr << "extractNextData: data.str()" <<  result.str() << std::endl;
		#endif 
	}
	#ifdef DEBUG
			std::cerr << "extractNextData finished!" <<  result.str() << std::endl;
		#endif 
	return result.str();
}




inline std::string 	extractNextBracedData(std::stringstream &Data)
{
	char a;
	std::stringstream result;

	extractChar('{',Data);
	int braceLevel=1;

	result << '{';
	
	while (true )
	{
		stripComments(Data);
		if (Data.eof())
			throw "Unexpected eof in curly braced data";
		Data   >> a;
		result << a ;

		if (a=='}' && --braceLevel==0)
			break;

		if (a=='{')
			braceLevel++;
	}
	return result.str();
}




/** @brief extracts a parametername,  from a input stream
 
 \verbatim Grammatic: <name>=<value> ; <comment>
 <comment>={}\n | {--commentText}}n
 \endverbatim

@todo : maybe  using an free parser generator is an optimal and higly portable solution instead of writing unreadable, complex, errneous and unflexible hardcoded parser
*
* todo : exclude further illegal chars - ',' etc.
*/
template <class _istream>
std::string 	readParamName(_istream & inputData)
{
	char a;

	std::stringstream parameterName;
	bool beginCommentPossible=false;
	bool containsWSPossible=false;

	while (!inputData.eof() && !inputData.fail() )
	{
		a=inputData.peek();
		if (a=='\n')
		{
			std::cerr << "ReadParamName: unexpected end of ParamName";
			throw "ReadParamName: line break not allowed in parameter names";
		}
		if (a=='=')
		{
			break;
		}
		if (a==';')
		{
			std::cerr << "ReadParamName: unexpected End of ParamName";
			throw "ReadParamName: unexpected End of ParamName";
		}

		if (isspace(a))
		{
			inputData >> std::ws;
 			if ( parameterName.str().length()>0)
			{
				containsWSPossible=true;
			}
			continue;
		}
		if (a=='-')
		{
			if (!beginCommentPossible)
		 		beginCommentPossible=true;
			else
			{
				std::cerr << "ReadParamName: unexpected begin of a comment";
				throw "ReadParamName: unexpected begin of a comment";
			}
		}
		else
		{
			beginCommentPossible=false;
		}
	
		if (containsWSPossible)
		{
			std::cerr << "ReadParamName: Parameter name contains white spaces";
			throw "ReadParamName: Parameter name contains white spaces";
		}
		containsWSPossible = false;
		inputData >> std::noskipws >> a;
		parameterName << a;
			
	}
	if (inputData.eof() || inputData.fail() )
	{
		std::cerr << "ReadParamName: unexpected end of data";
		throw "ReadParamName: unexpected end of data";
	}
	return parameterName.str();
}

/// @todo eskaped anfuehrungszeichen werden von readStringValue nicht erkannt
template <class _istream>
std::string readStringValue(_istream & inputData, std::string paramName)
{
	std::stringstream 	parameterValue;
	char a;
	
	a=inputData.peek();
	assert (a=='\"');
	
	inputData >> a;
	
	parameterValue  <<  a;

	while (!inputData.eof() && !inputData.fail() )
	{
		a=inputData.peek();
		if (a=='\"')
		{
			inputData >> a;
			parameterValue  <<  a;
			return parameterValue.str();
		}
		inputData >> std::noskipws >> a;
		parameterValue  <<  a;
	}
	throw "readStringValue: unexpected end of data";
}

/** @brief extracts a parameter value  from a input stream<br>

 Grammatic: <br>
%<name%>=%<value%> ; %<comment%>
 %<comment%>={}\n | {--commentText}}n

@todo : maybe  using an free parser generator is an optimal and higly portable solution instead of writing 
unreadable, complex, errneous and unflexible hardcoded parser
*
*/
template <class _istream>
std::string readParamValue(_istream & inputData, std::string paramName)
{
	char a;

	std::stringstream 	parameterValue;
	bool 			beginCommentPossible = false;

	std::string wsCache="";
	
	inputData >> std::ws;
//	cerr << "readParamValue of " << paramName << endl;
	while (!inputData.eof() && !inputData.fail() )
	{

 		//cerr << "parameterValue.str() " <<  parameterValue.str() << endl;
		a=inputData.peek();
		if (a=='\"')
		{
			return readStringValue( inputData, paramName);
		}
		if (a=='\n')
		{
			beginCommentPossible=false;
			inputData >> std::noskipws >> a;
			if (isspace(a))
			{
				wsCache+= a;
			}
			else
			{
				parameterValue  << wsCache; wsCache="";
				parameterValue  <<  a;
			}
			continue;
		}
		if (a=='=')
		{
			std::string rest;
			inputData >> rest;
			//cerr << "current: " << parameterValue.str() ;
			
		//	cerr << "rest: " << rest ;
			throw "ReadParamValue: unexpected assignment in parameter value";
		}
		if (a==';')
			break;
		if (a=='-')
		{
				bool wasComment=false;
				stripComment(inputData,wasComment);
				if (wasComment)
					continue;
		}
		inputData >> std::noskipws >> a;
		if (isspace(a))
		{
			wsCache+= a;
		}
		else
		{
			parameterValue  << wsCache; wsCache="";
			parameterValue  <<  a;
		}
//		std::cerr << "parameterValue" << parameterValue.str() << endl;
	}
	if (inputData.eof() || inputData.fail() )
		throw "readParamValue: unexpected end of data";
	return parameterValue.str();
}



inline int countSubGroups(std::stringstream &polynomStream)
{
	std::string strCopy = polynomStream.str();

	int braceLevel=0;
	int subGroupCount=0;

	for (size_t strPos=0;strPos<strCopy.length();strPos++)
	{
		if (strCopy[strPos]=='{')
		{
			if (braceLevel==1)
				subGroupCount++;
			braceLevel++;
			
		}
		if (strCopy[strPos]=='}')
		{
			braceLevel--;
			
		}
		if (braceLevel<0)
			throw "countSubGroups: too much closing braces!";
	}
	if (braceLevel!=0)
			throw "countSubGroups: Error in curly brace nesting";
	return subGroupCount;
}


inline int countGroups(std::stringstream &polynomStream)
{
	std::string strCopy = polynomStream.str();

	int braceLevel=0;
	int groupCount=0;

	for (size_t strPos=0; strPos<strCopy.length(); strPos++)
	{
		if (strCopy[strPos]=='{')
		{
			if (braceLevel==0)
				groupCount++;
			braceLevel++;
			
		}
		if (strCopy[strPos]=='}')
		{
			braceLevel--;
			
		}
		if (braceLevel<0)
			throw "countSubGroups: too much closing braces!";
	}
	if (braceLevel!=0)
			throw "groupCount: Error in curly brace nesting";
	return groupCount;
}


inline int countElements(std::stringstream &str)
{
	std::string strCopy = str.str();

	int braceLevel=0;
	int elementCount=0;

	for (size_t strPos=0;strPos<strCopy.length();strPos++)
	{
		if (strCopy[strPos]=='{')
		{
			braceLevel++;
		}
		if (strCopy[strPos]==',')
		{
			if (braceLevel==1)
				elementCount++;
		}
		if (strCopy[strPos]=='}')
		{
			if (braceLevel==1)
				elementCount++;
			braceLevel--;
		}
		if (braceLevel<0)
			throw "countSubGroups: too much closing braces!";
	}
	if (braceLevel!=0)
			throw "countElements: Error in curly brace nesting";
	return elementCount;
}




