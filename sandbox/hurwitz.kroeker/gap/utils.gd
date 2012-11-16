#############################################################################
##
#W hurwitzUtils                                                  Jakob Kröker
##                                                               
##
#H   @(#)$Id$
##
#Y Copyright (C) 2012, Laurent Bartholdi
##
#############################################################################
##
##
##  Implements following set of functions:
##
##  -list operation (flatten list)
##  -get/set polynomial coefficients
##  -evaluate polynomials
##  -derive polynomial lists (jacobian) 
##  -coercing (nested lists) of polynomials and scalars to different rings 
##   or fields.
##  -utils for factoring univariate polynomials
##
##############################################################################

# todo : define DeclareGlobalRecordFunction, DeclareGlobalRecordOperation!


BindGlobal("@FR@Utils" , rec() ) ;
DeclareGlobalFunction( "@FR@Utils.");
@FR@Utils.Tests    := rec();
@FR@Utils.Internal := rec();


# todo : replace Null@FR with Null@FR. 
# Introduce Null@FR to  indicate an uninitialized pointer/variable
if not IsBound( Null@FR ) then
    #BindGlobal("Null@FR", MakeImmutable( []) );
    BindGlobal("Null@FR", MakeImmutable( rec(value := "uninitialized data" ) ));
fi;

Assert( 0, false = IsMutable(Null@FR) );




DeclareGlobalFunction("InstallGlobalRecordFunctionOrMethod@FR");
DeclareGlobalFunction("InstallGlobalRecordFunction@FR");
DeclareGlobalFunction("InstallGlobalRecordFunctionEx@FR");
DeclareGlobalFunction("ReInstallGlobalRecordFunction@FR");

DeclareGlobalFunction("InstallGlobalRecordOperation@FR");
DeclareGlobalFunction("ReInstallGlobalRecordOperation@FR");


DeclareGlobalFunction("InstallGlobalRecordMethodEx@FR");
DeclareGlobalFunction("InstallGlobalRecordMethod@FR");
DeclareGlobalFunction("InstallGlobalRecordOtherMethod@FR");



#################################### LIST OPERATIONS ##################################

# flattens a list inplace. Example: [ [1,[2] ],[1] ] changes to [1,[2], 1] changes to [1,2,1] .
DeclareOperation("FlattenList@FR", [ IsList ] );
@FR@Utils.FlattenList:=FlattenList@FR;

DeclareOperation("FirstElement@FR", [ IsList ] );
@FR@Utils.FirstElement:=FirstElement@FR;

DeclareOperation("LastElement@FR", [ IsList ] );
@FR@Utils.LastElement:=LastElement@FR;



#################################### POLYNOMIAL DIFFERENTIATION ############################


# Jacobian: compute ( d [fktlist_i] / d [indeterminants_j] ). Parameters:  ( fktlist, indeterminants )
# todo: maybe add indeterminate information to the returned result.
DeclareGlobalFunction("Jacobian@FR");
@FR@Utils.Jacobian:=Jacobian@FR;

#################################### GET/SET COEFFICIENTS ##################################


# checks if an object is a monomial.
# Note: existing operation 'IsMonomial' probably means 'IsMonomialGroup'
DeclareOperation( "IsMonomial@FR",  [IsObject] ); 
# DeclareProperty( "IsMonomial@FR",  IsObject );
@FR@Utils.IsMonomial:=IsMonomial@FR;

# get coefficient for a specified monomial. Parameters: ( polynomial, monomial )
DeclareOperation( "MonomialCoefficient@FR", [ IsPolynomial, IsPolynomial] );
@FR@Utils.MonomialCoefficient:=MonomialCoefficient@FR;

# get coefficients [of specific monomials]. Parameters: ( polynomial, [monomials] )
# To get all nonzero coefficients just omit the second parameter.
# return value: coefficients
DeclareOperation("Coefficients@FR", [ IsPolynomial, IsList ] );
@FR@Utils.Coefficients:=Coefficients@FR;

# return value: [ coefficients, corresponding monomials ]
DeclareOperation("CoefficientsEx@FR", [ IsPolynomial, IsList ] );
@FR@Utils.CoefficientsEx:=CoefficientsEx@FR;

#################################### EVALUATE POLYNOMIALS #################################

# substitute indeterminates in the tensor by corresponding values. 
#Parameter: (tensor, indeterminates, values)
DeclareGlobalFunction("EvalPolynomialTensor@FR");
@FR@Utils.EvalPolynomialTensor:=EvalPolynomialTensor@FR;

DeclareGlobalFunction("SUBSTITUTE_POLYNOMIAL_COEFFICIENTS@FR");
@FR@Utils.SUBSTITUTE_POLYNOMIAL_COEFFICIENTS:=SUBSTITUTE_POLYNOMIAL_COEFFICIENTS@FR;

#################################### POLYNOMIAL AND SCALAR COERCION ########################

# note: coercion works only for some special cases and therefore  should probably not be exported. 
# (swith naming to all upper case? )

# coerce a scalar to a specific ring
#parameters: (scalar, destRing)
DeclareGlobalFunction("CoerceScalar@FR");
@FR@Utils.CoerceScalar:=CoerceScalar@FR;

# coerce a polynomial to a specific ring
# parameters: (polynomial, destRing)
DeclareGlobalFunction( "CoercePolynomial@FR" );
@FR@Utils.CoercePolynomial:=CoercePolynomial@FR;

# coerce a nested list of polynomials or scalars to a specific ring
# parameters: ( tensor,  destRing ). 
# limitations: accepts only scalars or polynomial as basic elements.
DeclareGlobalFunction( "CoerceTensor@FR" );
@FR@Utils.CoerceTensor:=CoerceTensor@FR;

DeclareSynonym("PromoteScalarTensor@FR", CoerceTensor@FR);
DeclareSynonym("CoerceScalarTensor@FR", CoerceTensor@FR);
DeclareSynonym("CoercePolynomialTensor@FR", CoerceTensor@FR);
@FR@Utils.PromoteScalarTensor:=PromoteScalarTensor@FR;
@FR@Utils.CoerceScalarTensor:=CoerceScalarTensor@FR;
@FR@Utils.CoercePolynomialTensor:=CoercePolynomialTensor@FR;


#################################### POLYNOMIAL PROPERTIES ##################################

DeclareOperation( "CountPolynomialVariables@FR", [IsPolynomial] );
DeclareSynonym( "IndeterminateNumber@FR", CountPolynomialVariables@FR);
DeclareSynonym( "NumberOfIndeterminates@FR", CountPolynomialVariables@FR);    

#DeclareSynonym( "IndeterminateNumber@FR",  IndeterminateNumberOfUnivariateRationalFunction);
@FR@Utils.IndeterminateNumber:=IndeterminateNumber@FR;
@FR@Utils.CountPolynomialVariables:=CountPolynomialVariables@FR;


#DeclareGlobalFunction("Degree@FR");
#@FR@Utils.Degree:=Degree@FR;

################################# POLYNOMIAL FACTORS AND PRODUCTS ##########################


DeclareOperation("IsPower@FR", [IsList] );
@FR@Utils.IsPower:=IsPower@FR;

# Parameters: (base, exponent)
DeclareGlobalFunction( "CreatePower@FR" );
@FR@Utils.CreatePower:=CreatePower@FR;

# todo: introduce 'Power' and 'Product' objects 

# note:  in following power means a pair [ base, exponent ].
#        and a product means a product (list) of powers.
 
# return distinct monic factors (no constants) of a polynomial.
DeclareOperation("DistinctMonicFactors@FR", [IsPolynomial ]);
@FR@Utils.DistinctMonicFactors:=DistinctMonicFactors@FR;

# computes the value of a product
DeclareGlobalFunction("PRODUCT_VALUE@FR");
@FR@Utils.PRODUCT_VALUE := PRODUCT_VALUE@FR;


# factors a homogenized univariate polynomial over rationals or finite fields into  power factors with distinct bases. 
DeclareOperation("UNIQUE_PRODUCT_HOMOGENIZED@FR", [ IsPolynomial ]);

# factors an univariate polynomial over rationals or finite fields into  power factors with distinct bases. 
# see also 'Factors'
DeclareOperation("UNIQUE_PRODUCT@FR", [ IsPolynomial ]);
@FR@Utils.UNIQUE_PRODUCT := UNIQUE_PRODUCT@FR;




DeclareOperation("UNIQUE_PRODUCT_1@FR",[ IsPolynomial ]);  # just a different implementation.


# removes constant factors from a list of powers.
DeclareGlobalFunction("REMOVE_CONSTANT_FACTORS@FR");  
@FR@Utils.REMOVE_CONSTANT_FACTORS := REMOVE_CONSTANT_FACTORS@FR;

# sort a list of powers   by exponent desc
DeclareGlobalFunction("SORT_POWERS_BY_EXPONENT@FR");
@FR@Utils.PRODUCT_VALUE := PRODUCT_VALUE@FR;




#E hurwitzUtils.gd . . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
