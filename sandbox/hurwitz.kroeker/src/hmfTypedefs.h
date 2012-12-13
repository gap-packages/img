#pragma once

#include <assert.h>
#include <iostream>
#include <cstdio>
#include <string>

#include "fast_Ring.h"
#include "fastNumber.h"
#include "typedefs.h"
#include "polynomdefs.h"


#include "random.h"
#include "basicNumber.h"
#include "fastNumber.h"

#include "polynomialRing.h"

#include <algorithm>
#include <map>

typedef std::map<std::string, int64_t> HashMapType;

#ifndef HURWITZ_MAXCHAR
    #define HURWITZ_MAXCHAR  123
#endif

#ifndef HURWITZ_SCALARTYPE
    #define HURWITZ_SCALARTYPE int8_t
#endif

typedef number_eps0< HURWITZ_MAXCHAR, HURWITZ_SCALARTYPE >  CoeffType   ;

#include "polynom.h"

typedef fast_Ring< CoeffType, kdefs_zahl_x < HURWITZ_MAXCHAR ,0 > >   defined_Field_Type;

typedef UnivariatePolynomialRing< polynomx<defined_Field_Type::ElementType > ,defined_Field_Type > TPolRingType;

// TPolFactorPowerType::first: polynomial, TPolFactorPowerType::second: exponent
typedef   std::pair<TPolRingType::Element , uint>         TPolFactorPowerType;
