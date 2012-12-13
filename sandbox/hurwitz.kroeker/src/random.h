#ifndef RANDOM1_H
#define RANDOM1_H




/** @file random.h
*
* @ingroup RandomGenerator

* @brief random number generator of L'Ecuyer with Bays-Durham shuffle and added safeguards (numeric recipes)
*
* @todo seed ist fuer die random-Funktion  nicht dokumentiert.
*/

#include <stdint.h>
#include <assert.h>

/* note #undef's at end of file */
/*
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


// Long period (> 2.0e18) random number generator of L'Ecuyer with
// Bays-Durham shuffle and added safeguards.  Returns a uniform
// random deviate between 0.0 and 1.0 (exclusive of the endpoints).

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=(float)(AM*iy)) > RNMX) return (float) RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
*/


// Long period (> 2.0e18) random number generator of L'Ecuyer with
// Bays-Durham shuffle and added safeguards.  Returns a uniform
// random deviate between 0.0 and 1.0 (exclusive of the endpoints).

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


/**
* @brief Long period (> 2.0e18) random number generator of L'Ecuyer with<br>
* Bays-Durham shuffle and added safeguards.  Returns a uniform<br>
* random deviate between 0.0 and 1.0 (exclusive of the endpoints).
*
* @note probably cannot be inlined because of the static variables.
*
* @ingroup RandomGenerator
*/
double ran2(long *idum);

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/// @todo ein Random-Objekt sollte den seed selbst speichern.

/** @brief returns a random value between zero and <b>max</b>,  Long period (> 2.0e18) (see ran2())
*
* @ingroup RandomGenerator
*/
inline unsigned short random (long *seed, unsigned short max) 
{
	#ifdef DEBUG 
		std::cerr << "----------- " << std::endl;
		std::cerr << "Random call " << std::endl;
		std::cerr << "seed "<<  *seed << std::endl;
		std::cerr << "max "<<  max << std::endl;
	#endif 
	
	//float zahl = 1.0; //Attention: float for ran1, double for ran2 !!!
	double zahl = 1.0;
	int counter=0;
	while (zahl == 1.0)
	{
		zahl = ran2(seed);
		counter++;
	}

	assert(counter<=1);
	// Zahl ist zwischen 0 (eingeschl.) und 1 (ausgeschlossen)
	zahl = zahl * (max+1);
	// Zahl ist zwischen 0 (eingeschl.) und max+1 (ausgeschlossen)

	//also muss Zahl nur noch abgerundet werden
	#ifdef SAFE
		assert((unsigned short)(zahl - 0.0) == (unsigned short)(zahl  ) );
		assert((unsigned short)(zahl  )<=max );
	#endif

	#ifdef DEBUG 
		
		std::cerr << "random "<<  (unsigned short)(zahl - 0.0) << std::endl;
		std::cerr << "Random call " << std::endl;
		std::cerr << "----------- " << std::endl;
	#endif 
	return  (unsigned short) (zahl - 0.0);
}


inline uint32_t randomUInt32 (long *seed, uint32_t max) 
{
	#ifdef DEBUG 
		std::cerr << "----------- " << std::endl;
		std::cerr << "Random call " << std::endl;
		std::cerr << "seed "<<  *seed << std::endl;
		std::cerr << "max "<<  max << std::endl;
	#endif 
	
	//float zahl = 1.0; //Attention: float for ran1, double for ran2 !!!
	double zahl = 1.0;
	int counter=0;
	while (zahl == 1.0)
	{
		zahl = ran2(seed);
		counter++;
	}

	assert(counter<=1);
	// Zahl ist zwischen 0 (eingeschl.) und 1 (ausgeschlossen)
	zahl = zahl * (max+1);
	// Zahl ist zwischen 0 (eingeschl.) und max+1 (ausgeschlossen)

	//also muss Zahl nur noch abgerundet werden
	#ifdef SAFE
		assert((unsigned short)(zahl - 0.0) == (unsigned short)(zahl  ) );
		assert((unsigned short)(zahl  )<=max );
	#endif

	#ifdef DEBUG 
		
		std::cerr << "random "<<  (unsigned short)(zahl - 0.0) << std::endl;
		std::cerr << "Random call " << std::endl;
		std::cerr << "----------- " << std::endl;
	#endif 
	return  (uint32_t) (zahl - 0.0);
}


static long CF_MM_s1 = 1;
static long CF_MM_s2 = 1;

#define MODMULT(a, b, c, m, s) q = s/a; s = b*(s-a*q)-c*q; if (s < 0) s += m;

static double combinedLCG(void)
{
    long q, z;

    MODMULT(53668, 40014, 12211, 2147483563L, CF_MM_s1)
    MODMULT(52774, 40692, 3791, 2147483399L, CF_MM_s2)
    z = CF_MM_s1 - CF_MM_s2;
    if (z < 1)
        z += 2147483562;
    return z * 4.656613e-10;
}


static void initLCG(long init_s1, long init_s2)
{
    CF_MM_s1 = init_s1;
    CF_MM_s2 = init_s2;
}


#endif
