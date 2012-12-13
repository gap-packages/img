#############################################################################
##
#W hurwitz.gd                                               Laurent Bartholdi
##                                                               Jakob Kröker
##
#H   @(#)$Id$
##
#Y Copyright (C) 2012, Laurent Bartholdi
##
#############################################################################
##
##  Solutions to the Hurwitz problem
##
#############################################################################


# depends on :   package 'Float',  "hurwitz/gap/utils",  "hurwitz/gap/padicLift"
#



Hurwitz@FR := rec();

DeclareGlobalFunction( "Hurwitz@FR.");
Hurwitz@FR.Tests    := rec();
Hurwitz@FR.Internal := rec();


#DeclareGlobalFunction( "Hurwitz@FR.Internal");
DeclareGlobalFunction( "Hurwitz@FR.Internal.");


# hack:  ( switch HURWITZMAPSEARCHBIN for debug/development).
# Problem will disappear in case build system is ready ( TODO )
# HURWITZ_MAP_SEARCH_BIN@FR := Filename( [Directory("/home/kroeker/rationalFunctionSearch/c-program/bin")],  "hurwitzMapSearchForGAP.mathpc26" );

HURWITZ_MAP_SEARCH_BIN@FR := Filename(  DirectoriesPackagePrograms("fr") ,  "hurwitzMapSearch" );


########################################## FIND HURWITZ MAP OVER A FINITE FIELD #####################################################


# create a representation for multiplicity structure of a polynomial
# Parameter: integer partition. 
DeclareGlobalFunction("Shape@FR");
Hurwitz@FR.Shape := Shape@FR;
DeclareGlobalFunction( "Hurwitz@FR.Shape" );

DeclareGlobalFunction("IsShape@FR");
Hurwitz@FR.IsShape := IsShape@FR;
#DeclareGlobalFunction( "Hurwitz@FR.IsShape" );


# computes the shape of an univariate polynomial. 
# A shape is here a desc-ordered list of root multiplicities.
# Parameters: (polynomial, [expected degree] )
# the optinal parameter 'expected degree' is required to determine the shape correctly if the polynomial has infinity root factor.
DeclareOperation("ComputeShape@FR", [ IsPolynomial, IsInt ] );
Hurwitz@FR.ComputeShape := ComputeShape@FR;
DeclareGlobalFunction( "Hurwitz@FR.ComputeShape" );


# get the multiplicity of a univariate polynomial root
# Parameters: ( polynomial  over a finite field,  root, [ poldegree] ) 
# poldegree is passed to get the correct multiplicity of the infinity root.
DeclareOperation( "RootMultiplicity@FR", [ IsUnivariatePolynomial, IsObject, IsInt ] );
Hurwitz@FR.RootMultiplicity := RootMultiplicity@FR;
DeclareGlobalFunction( "Hurwitz@FR.RootMultiplicity" );

# find a solution for a Hurwitz map problem over a finite field. 
#Parameters: ( prime field, permutations, criticalValues )
#      or :  ( prime field, shapes,       criticalValues, strictNormalization(bool) )
# preconditions: product of the permutations is =1; all shapes have same degree; number of shapes/permutations and criticalValues matches.
DeclareOperation( "FindHurwitzMapModPrime@FR", [ IsPrimeField, IsList, IsList ] );
Hurwitz@FR.FindHurwitzMapModPrime := FindHurwitzMapModPrime@FR;
DeclareGlobalFunction( "Hurwitz@FR.FindHurwitzMapModPrime" );


# compute the search space size for a Hurwitz map search problem over a given finite field.
# Parameters: ( prime field, permutations, criticalValues )
DeclareOperation( "HurwitzMapSearchSpaceSize@FR", [ IsPrimeField, IsList, IsList, IsBool ] );
Hurwitz@FR.HurwitzMapSearchSpaceSize := HurwitzMapSearchSpaceSize@FR;
DeclareGlobalFunction( "Hurwitz@FR.HurwitzMapSearchSpaceSize" );




##################### LIFT FINITE FIELD HURWITZ MAP TO RATIONALS/COMPLEX NUMBERS ###########################


# create a Hurwitz map search problem 
# Parameters: ( partitions, criticalValues, strictNormalization (bool) )
# or 
# Parameters: ( partitions, criticalValues, normalizationRules)
# expect first critival values to be infinity, zero, one  and the following to be 'rational number approximations'
# 'rational number approximations' :  pairs of real and imaginary parts  of an rational approximation .
 DeclareOperation( "HurwitzMapSearchProblem@FR", [IsList, IsList, IsBool] );
Hurwitz@FR.HurwitzMapSearchProblem := HurwitzMapSearchProblem@FR;

# Parameter:  ( hurwitzMapSearchProblem, polynomial list W_i mod prime , finiteField, liftOptions )
# for liftOptions see LiftOptions@FR
DeclareGlobalFunction( "ApproxComplexHurwitzMaps@FR");
Hurwitz@FR.ApproxComplexHurwitzMaps := ApproxComplexHurwitzMaps@FR;


########## functions required for customizing lift: 
########## e.g using  normalization different from default ( preimage( inf,0, 1) = (inf,0,1) )

# given a rational approximation of a complex root a+ib (a pair of real and imaginary part approximations ),
# create a minimal polynomial for roots [a+ib, a-ib]  over integers  using the second parameter 'variable' as indeterminate.
DeclareGlobalFunction( "RationalMinPolyFromRootApprox@FR");
Hurwitz@FR.RationalMinPolyFromRootApprox := RationalMinPolyFromRootApprox@FR;


# create a normalization rule 
# Parameters_: (polynomialId, multiplicity, rootValue)
DeclareGlobalFunction( "NormalizationRule@FR" );
Hurwitz@FR.NormalizationRule := NormalizationRule@FR;  


DeclareGlobalFunction( "IsNormalizationRule@FR" );
Hurwitz@FR.NormalizationRule := NormalizationRule@FR;  

# Parameters ( polTuple, finiteField, HurwitzMapSearchProblem )
# rename to PolSet HurwitzMap or ReducedHurwitzMap ?
DeclareGlobalFunction( "HurwitzMapLifter@FR" );
Hurwitz@FR.HurwitzMapLifter := HurwitzMapLifter@FR;



MakeImmutable( Hurwitz@FR.Tests );
MakeImmutable( Hurwitz@FR.Internal );

MakeImmutable( Hurwitz@FR );
MakeReadOnlyGlobal("Hurwitz@FR");
