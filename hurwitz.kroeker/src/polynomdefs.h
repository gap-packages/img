

#pragma once

#include "CompileFunctions.h"

#include <assert.h>

enum P_or_QPolynom
    {
      PCoefficient,
      QCoefficient,
};


// siehe zum besseren Verstaendnis "Aspektorientierte Programmierung in C++: Teil 2, Multiple Aussichten"
// http://www.heise.de/ix/artikel/2001/09/142/09.shtml
// Traits heisst uebrigens Charakterzuege
// see also http://www.oonumerics.org/tmpw00/
// http://www.it.neclab.eu/~berti/generic/

/** 
@brief contains template parameter definition for fast_polynomXY (Traits), 
*     parametrized by characteristic(field) during compile time. The design was more try and error than
*     applying common metaprogramming concepts, but this schould be the right way
*
*
* @description
*
* der Aufbau ist wie folgt: 
* polynomdefs definiert (zur Kompilezeit) beispielsweise den inneren Datenaufbau von improved_polynom (polynom in zwei Variablen)
*<br>
* Erreichte Optimierung
* Durch die Shift-Operation in getPairIndex wird der Zugriff
* auf eine Polynomkoeffizient <br>
* (k*(x^a)*(y^b), a und b \in{0..DEGREE} )  <br>
* beschleunigt.
* <br>
* Nachteil: ist DEGREE keine Zweierpotenz, so muss für die Polynomkoeffizienten
* ein unnoetig groesseres Array definiert werden, in welchem nicht alle Eintraege belegt sind. 
* Dies führt wieder zu vermehrten Cache Misses
*
* improved_polynom haengt aber auch vom Datentyp des Polynomkoeffizienten ab, ich weiss
* momentan nur nicht, wie ich den Dateintyp als typedef-Variable in polynomdefs
* festlege
*
* ein Problem ist der dynamische Aspekt: wie 
* kann man 'improved_polynom' auch mit einem Objekt und nicht mit statischen 
* Daten initialiesieren?
*
* @param DEGREE ist der maximal verwendete Grad des Polynoms (in zwei Variablen.)
*
* @todo: Ueberpruefung auf DEGREE kleiner 128, wegen 	"static const short pshift" - wieso, DAS ist kein Problem, aber eine wellDefined-Funktion sollte programmiert werden.
*/
/// @todo Vorsicht beim Schiften! muss ich beim Shift die Variablen vorher in ein int konvertieren ja/nein?
/// @todo Korrektheit: Polynomdefs alle falsch, nur durch die Tatsache, dass nextpow2num(3)=4 ist (hoffentlich), passen die Größen doch und es gibt keine
///   Speicherschutzverletzung
//polynomdefs und MAtrixDefs sind genau gleich -> Eine der Definitionen koennte wegfallen!
template <int DEGREE>
class polynomdefs
{
	///due to compile problems with some compilers enums are used instead of regular member variables
private:
	enum { 
		pshift_m = needbits<DEGREE>::value
	};
public:
	enum {	/// maximal zulaessige Grad eines Polynommonoms mit der Einstellung DEGREE
		maxdegree_m	= nextpow2num<DEGREE>::value - 1
	};
	int getSize()  const { return size_m; }
	enum {/// Anzahl bytes, die von einem arrayindex der Polynomkoeffizientenliste belegt werden
		/// @todo wieso nicht (DEGREE <<needbits<DEGREE>::value |degree) +1 ? nextpow2num ist aber nicht falsch (nur zu viel) und belegt eine gerade Anzahl von Plätzen.
		size_m		= (nextpow2num<DEGREE>::value)*(nextpow2num<DEGREE>::value)
	};
	inline  static bool wellDefined()
	{
		return (DEGREE<128);
	};
	/// compute a*(2^pshift)+b
	inline  static short getPairIndex(const short a, const short b) 
	{
		#ifdef COUNT
			bitwiseShift	+= 1;
			bitwiseOR	+= 1;
		#endif	
		
		//return (((int)a<<pshift)|(int)b);
		return ( ((int)a<<needbits<DEGREE>::value) | (int)b );
	}
};

/// @brief polynomdefsNoShift: bei dieser Datenstruktur liegen die Polynomkoeffizienten noch nicht hintereinander.
/// @note es ging bei dieser Definition um den Test, ob die Schiftoptimierungen überhaupt noch was an Performance bringen, 
/// nachdem in dem fast_frommer-Code auf die meisten Polynomkoeffizienten nicht mehr über Addressberechnung, sondern
/// sequentiell zugegriffen wird. Fazit: nein, kaum noch ein Performance-Unterschied.
template <int DEGREE>
class polynomdefsNoShift
{
	///due to compile problems with some compilers enums are used instead of regular member variables
private:
	 
public:
	enum {	/// maximal zulaessige Grad eines Polynommonoms mit der Einstellung DEGREE
		maxdegree_m	= DEGREE
	};

	enum {/// Anzahl Bytes, die von einem arrayindex der Polynomkoeffizientenliste belegt werden
		size_m	=  (DEGREE+1)*(DEGREE+1)
	};
	int getSize()  const { return size_m; }
	/// compute a*(2^pshift)+b
	inline  static short getPairIndex(const short a, const short b) 
	{
		return ( ( (int)(a+b)*(DEGREE+1)) + (int)b );
	}
};

class dynamicQuadraticMatrixDefsNoShift
{
private:
	int 	size_m;	///< groesse des zu reservierenden Speicherblcks

	dynamicQuadraticMatrixDefsNoShift();
 
public:
	const int 	dim_m;///< anzahl der spalten/zeilen
	

	inline int getSize()  const { return size_m; }

	dynamicQuadraticMatrixDefsNoShift(int dimension):	size_m((dimension+1)*(dimension+1) ) ,
										dim_m(dimension+1)	
										
	{
	}

	inline  short getPairIndex(const short a, const short b) const
	{
		return (  (int)(a*dim_m) + (int)b );
	}
};


class dynamicPolynomdefsNoShift
{

private:
	  int 	size_m;
	dynamicPolynomdefsNoShift();
public:
	 int 	maxdegree_m;
	 int 	maxdegreePlusOne_m;
	

	inline int getSize()  const { return size_m; }

	dynamicPolynomdefsNoShift(int degree):	size_m((degree+1)*(degree+1) ),
								maxdegree_m(degree),
								maxdegreePlusOne_m(degree+1)
	{
	}

	inline  short getPairIndex(const short a, const short b) const
	{
		return ( ( (int)(a+b)*maxdegreePlusOne_m) + (int)b );
	}
};


/// @note für den Test: die offsets sollten 0,1,3,6,10,15,... u.s.w. sein.
class dynamicPolynomdefsNoShiftNoMemoryHoles
{
private:
	/// kann man daraus ein const int array machen?
	int* 	offsets_m;
 	int 	size_m;
	dynamicPolynomdefsNoShiftNoMemoryHoles();
	
	

public:
	 int 	maxdegree_m;
	 int 	maxdegreePlusOne_m;
	
	// wegen short getIndex()
	bool wellDefined()
	{
		assert(size_m<32000);
	}

	inline int getSize()  const { return size_m; }

	int computeSize(int degree)
	{
		int size=0;
		for (int currDegree=0;currDegree<=degree; currDegree++)
		{
			offsets_m[currDegree] = size;
			for (int yExp=0;yExp<=currDegree; yExp++)
				size++;
		}
		return size;
	}

	dynamicPolynomdefsNoShiftNoMemoryHoles(int degree):	
							offsets_m(new int[degree+1]),
							size_m(computeSize(degree)),		
							maxdegree_m(degree),
							maxdegreePlusOne_m(degree+1)
	{
	}

	inline  short getPairIndex(const short a, const short b) const
	{

		int offset=0;
		int monomDegree = a+b;
		for (int currDegree=0; currDegree< monomDegree; currDegree++)
		{
			for (int yExp=0;yExp<=currDegree; yExp++)
				offset++;
		}
		assert( offsets_m[ monomDegree ] == offset );
		offset = offset + b;
		return offset;
	}
};

/** @brief polynom traits*/
/// @todo zusaetzlich mit einem Datentyp parametrisieren? - nö, Typ ist der Exponent
/// @note diese Definitionen sind für die Verwendung als Matrix eher untauglich:
///       man Bedenke, dass wenn eine  NxN-Matrix gespeichert werden soll, als DEGREE-Parameter (N+N=2N) übergeben werden muss.
template <int DEGREE>
class polynomdefsNew
{
	///due to compile problems with some compilers enums are used instead of regular member variables
private:
	enum { 
		pshift_m = needbits<DEGREE>::value
	};
public:
	enum {	/// maximal zulaessige Grad eines Polynommonoms mit der Einstellung DEGREE
		maxdegree_m	= DEGREE
	};

	enum {/// Anzahl bytes, die von einem arrayindex der Polynomkoeffizientenliste belegt werden
		size_m		= (nextpow2num<DEGREE>::value)*(nextpow2num<DEGREE>::value)
	};

	inline static unsigned short getSize()
	{
		return (nextpow2num<DEGREE>::value)*(nextpow2num<DEGREE>::value);
	}
	/// compute a*(2^pshift)+b
	inline  static size_t getPairIndex(const unsigned short a, const unsigned short b) 
	{
		#ifdef COUNT
			bitwiseShift	+= 1;
			bitwiseOR	+= 1;
		#endif	
		// Vorsicht bein shiften! muss ich pshift sicherheitshalber in ein int konvertieren ja/nein?
		//return (((int)a<<pshift)|(int)b);
		return ( ((int)(a+b)<<needbits<DEGREE>::value) | (int)b );
	}
};



/** @brief polynompair traits 
 @todo Codereduzierung: polynompairdefskoennte wegfallen, wenn pair eingesetzt wird. */
template <int DEGREE>
class polynompairdefs
{
	///due to compile problems with some compilers enums are used instead of regular member variables
private:
	enum { 
		pshift_m = needbits<DEGREE>::valueplusone
	};
public:
	enum {	/// maximal zulaessige Grad eines Polynommonoms mit der Einstellung DEGREE
		maxdegree_m	= DEGREE
	};

	enum {/// Anzahl bytes, die von einem arrayindex der Polynomkoeffizientenliste belegt werden
		size_m		= 2*(nextpow2num<DEGREE>::value) * (nextpow2num<DEGREE>::value) 
	};

	int getSize()  const { return size_m; }
	/// compute a*(2^pshift)+b
	inline  static short 	getPairIndex(const short a, const short b) 
	{
		#ifdef COUNT
			bitwiseShift 	+= 1;
			bitwiseOR	+= 1;
		#endif	
		// Vorsicht bein shiften! muss ich pshift sicherheitshalber in ein int konvertieren ja/nein?
		//return (((int)a<<pshift)|(int)b);
		return ( ((int)(a+b) << (needbits<DEGREE>::valueplusone)) | (int)(b<<1) );
	}/// compute a*(2^pshift)+b


	inline  static short 	getGroupIndex(const short a) 
	{
		#ifdef COUNT
			bitwiseShift	+= 1;
			bitwiseOR	+= 1;
		#endif	
		// Vorsicht bein shiften! muss ich pshift sicherheitshalber in ein int konvertieren ja/nein?
		//return (((int)a<<pshift)|(int)b);
		return ( (int) a << (needbits<DEGREE>::valueplusone) );
	}
};



//maxdegree mei Matrixdef
// haengt nicht vom unteren Datenyp der MAtrix ab, aber legt wohl die 
// maximale groesse fest
/**
@brief 	contains template parameter definition for matrix classes (simulated by fast_polynomXY), 
*	parametrized by char(field) during compile time
*/
template <short CHAR>
class matrixdefs
{
private:
	enum {
		pshift_m = needbits<CHAR>::value
	};
public:
	///Statt Variablen wurden enums verwendet, weil z.B: der Intel Compiler damit nicht zurechtkommt
	enum {
		maxdegree_m = nextpow2num<CHAR>::value - 1
	};
	enum {/// Anzahl bytes, die von einem arrayindex der Polynomkoeffizientenliste belegt werden
		size_m	= (nextpow2num<CHAR>::value)*(nextpow2num<CHAR>::value)
	};
	int getSize()  const { return size_m; }
		/// compute a*(2^pshift)+b
		inline static  short getPairIndex(short a, short b) 
		{
			#ifdef COUNT
				bitwiseShift+=1;
				bitwiseOR+=1;
			#endif
			// Vorsicht bein shiften! muss ich pshift sicherheitshalber in ein int konvertieren ja/nein?
			return ( ((int)a<<needbits<CHAR>::value) | (int)b );

		}
};


/**
* @brief contains template parameter definition for fast_Ring class, 
*	 parametrized by static field characteristik during compile time
*
*
*
* @todo Parametrize with EPSPRECISION ?
*/
template <unsigned int CHAR, unsigned short EPSPREC>
class kdefs_zahl_x
{
	public:
	/*	static  const  short charakteristik;
		static  const  short epsPrecision; // umbenennen in base_num_epsPrecision???
	*/
			
	enum {		
		charakteristik_m = CHAR	
	};
	enum {		
		charakteristik_minus_one_m = CHAR-1
	};
	enum {
		epsPrecision_m	 = EPSPREC
	};
};

/* // this implementation does not compile with some compilers, therefore the 'enum' approach is used - leave it as syntax example!
template <int CHAR>
const  short kdefs_zahl_x<CHAR>::charakteristik(CHAR);

template <int CHAR>
const  short kdefs_zahl_x<CHAR>::epsPrecision(EPSPRECISION);
*/

