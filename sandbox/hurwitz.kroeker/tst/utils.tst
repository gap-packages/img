#############################################################################
##
#W  utils.tst                  FR Package               Jakob Kroeker
##
#H  @(#)$Id$
##
#Y  Copyright (C) 2012,  Laurent Bartholdi
##
#############################################################################
##
##  This file tests the polynomials and list utils for hurwitz package
##
#############################################################################

# following lines generated with "str:= @FR@Utils.CreateTestString(true); Print(str);"
#
#
#  @FR@Utils.Tests.TEST_FLATTEN_LIST : 
#
gap>     Assert( 0, [  ] = FlattenList@FR( [  ] ) );;
gap>     Assert( 0, [ 1, 2, 1 ] = FlattenList@FR( [ 1, [ 2, 1 ] ] ) );;
gap>     Assert( 0, [ 1, 2, [ 1 ] ] = FlattenList@FR( [ 1, [ 2, [ 1 ] ] ] ) );;
gap>     Assert( 0, [ [ 1 ], 1 ] = FlattenList@FR( [ [  ], [ [ 1 ] ], 1 ] ) );;
#
#
#  @FR@Utils.Tests.TEST_IS_MONOMIAL : 
#
gap> 
gap>     rng := PolynomialRing( ZmodnZ( 11 ), [ "x", "y" ] );;
gap>     indet := IndeterminatesOfPolynomialRing( rng );;
gap>     x := indet[1];;
gap>     y := indet[2];;
gap>     Assert( 0, IsMonomial@FR( x ) );;
gap>     Assert( 0, IsMonomial@FR( x * y ) );;
gap>     Assert( 0, not IsMonomial@FR( 2 * x * y ) );;
gap>     Assert( 0, not IsMonomial@FR( x + y ) );;
gap>     Assert( 0, not IsMonomial@FR( 3 ) );;
gap>     Assert( 0, not IsMonomial@FR( rng ) );;
#
#
#  @FR@Utils.Tests.TEST_MONOMIAL_COEFFICIENT : 
#
gap> 
gap>     rng := PolynomialRing( ZmodnZ( 11 ), [ "x", "y" ] );;
gap>     indeterminates := IndeterminatesOfPolynomialRing( rng );;
gap>     x := indeterminates[1];;
gap>     y := indeterminates[2];;
gap>     polynomial := (x ^ 4 - 4) ^ 3 * (4 * y ^ 2 + 2);;
gap>     Assert( 0, Z( 11 ) ^ 4 = MonomialCoefficient@FR( polynomial, x ^ 4 * y ^ 2 ) );;
gap>     Assert( 0, Zero( Z( 11 ) ) = MonomialCoefficient@FR( polynomial, x ^ 42 * y ^ 2 ) );;
gap>     Assert( 0, Z( 11 ) ^ 2 = MonomialCoefficient@FR( polynomial, x ^ 0 * y ^ 0 ) );;
#
#
#  @FR@Utils.Tests.TEST_COEFFICIENTS : 
#
gap> 
gap>     rng := PolynomialRing( ZmodnZ( 11 ), [ "x", "y" ] );;
gap>     indeterminates := IndeterminatesOfPolynomialRing( rng );;
gap>     x := indeterminates[1];;
gap>     y := indeterminates[2];;
gap>     polynomial := (x ^ 4 - 4) ^ 3 * (4 * y ^ 2 + 2);;
gap>     Assert( 0, [ Z( 11 ) ^ 4, Zero( Z( 11 ) ) ] = Coefficients@FR( polynomial, [ x ^ 4 * y ^ 2, x ^ 42 * y ^ 2 ] ) );;
#
#
#  @FR@Utils.Tests.TEST_JACOBIAN : 
#
gap> 
gap>     rng := PolynomialRing( Rationals, 2 );;
gap>     ind := IndeterminatesOfPolynomialRing( rng );;
gap>     x := ind[1];;
gap>     y := ind[2];;
gap>     scalar := 5 / 3;;
gap>     pol := scalar * x;;
gap>     jacobian := Jacobian@FR( [ pol, y ^ 2 ], ind );;
gap>     Assert( 0, jacobian = [ [ Derivative( pol, x ), Derivative( pol, y ) ], [ Derivative( y ^ 2, x ), Derivative( y ^ 2, y ) ] ] );;
#
#
#  @FR@Utils.Tests.TEST_COERCE_SCALAR : 
#
gap> 
gap>     scalar := 1 / 3;;
gap>     dstRing := Integers;;
gap>     dstRing := GF( 11 );;
gap>     CoerceScalar@FR( scalar, dstRing );;
gap>     dstRing := ZmodnZ( 11 );;
gap>     Assert( 0, One( dstRing ) * scalar = CoerceScalar@FR( scalar, dstRing ) );;
gap>     scalar := 23;;
gap>     dstRing := Integers;;
gap>     Assert( 0, One( dstRing ) * scalar = CoerceScalar@FR( scalar, dstRing ) );;
gap>     dstRing := GF( 11 );;
gap>     Assert( 0, One( dstRing ) * scalar = CoerceScalar@FR( scalar, dstRing ) );;
gap>     dstRing := ZmodnZ( 121 );;
gap>     Assert( 0, One( dstRing ) * scalar = CoerceScalar@FR( scalar, dstRing ) );;
#
#
#  @FR@Utils.Tests.TEST_COERCE_POLYNOMIAL : 
#
gap> 
gap>     rng := PolynomialRing( Rationals, 1 );;
gap>     ind := IndeterminatesOfPolynomialRing( rng );;
gap>     x := ind[1];;
gap>     scalar := 5 / 3;;
gap>     pol := scalar * x;;
gap>     baseField := ZmodnZ( 11 );;
gap>     dstRng := PolynomialRing( baseField, 1 );;
gap>     coercedPol := CoercePolynomial@FR( pol, dstRng );;
gap>     dstInd := IndeterminatesOfPolynomialRing( dstRng );;
gap>     expectedResult := dstInd[1] * Z( 11 ) ^ 6;;
gap>     Assert( 0, coercedPol = expectedResult );;
gap>     baseField := ZmodnZ( 121 );;
gap>     dstRng := PolynomialRing( baseField, 1 );;
gap>     coercedPol := CoercePolynomial@FR( pol, dstRng );;
gap>     dstInd := IndeterminatesOfPolynomialRing( dstRng );;
gap>     expectedResult := dstInd[1] * ZmodnZObj( 42, 121 );;
gap>     Assert( 0, coercedPol = expectedResult );;
gap>     CoerceScalar@FR( scalar, dstRng );;
#
#
#  @FR@Utils.Tests.TEST_EVAL_POLYNOMIAL_TENSOR : 
#
gap> 
gap>     rng := PolynomialRing( Rationals, 2 );;
gap>     ind := IndeterminatesOfPolynomialRing( rng );;
gap>     x := ind[1];;
gap>     y := ind[2];;
gap>     mat := [ [ 1 / 3 * x ^ 0, x ^ 0, x + y ] ];;
gap>     dstRng := ZmodnZ( 121 );;
gap>     evaluatedTensor := EvalPolynomialTensor@FR( mat, [ x, y ], [ ZmodnZObj( 1, 121 ), ZmodnZObj( 2, 121 ) ] );;
gap>     EvalPolynomialTensor@FR( mat, [ x, y ], [ ZmodnZObj( 1, 121 ), ZmodnZObj( 2, 121 ) ] );;
#
#
#  @FR@Utils.Tests.TEST_SUBSTITUTE_POLYNOMIAL_COEFFICIENTS : 
#
gap> 
gap>     rng := PolynomialRing( Rationals, 3 );;
gap>     ind := IndeterminatesOfPolynomialRing( rng );;
gap>     a := ind[1];;
gap>     b := ind[2];;
gap>     PREV_ITER_POLY_WARN := ITER_POLY_WARN;;
gap>     ITER_POLY_WARN := false;;
gap>     iterRng := PolynomialRing( rng, 2 );;
gap>     ITER_POLY_WARN := PREV_ITER_POLY_WARN;;
gap>     iterInd := IndeterminatesOfPolynomialRing( iterRng );;
gap>     x := iterInd[1];;
gap>     y := iterInd[2];;
gap>     pol := a * b * x + b * y;;
gap>     dstRng := PolynomialRing( Rationals, 2 );;
gap>     dstFam := FamilyObj( One( dstRng ) );;
gap>     result := SUBSTITUTE_POLYNOMIAL_COEFFICIENTS@FR( pol, ind, [ 2, 1, 0 ], dstFam );;
gap>     Assert( 0, CoercePolynomial@FR( result, iterRng ) = CoercePolynomial@FR( 2 * x + y, iterRng ) );;
#
#
#  @FR@Utils.Tests.TEST_COUNT_POLYNOMIAL_VARIABLES : 
#
gap> 
gap>     rng := PolynomialRing( ZmodnZ( 11 ), [ "x", "y" ] );;
gap>     indeterminates := IndeterminatesOfPolynomialRing( rng );;
gap>     x := indeterminates[1];;
gap>     y := indeterminates[2];;
gap>     Assert( 0, CountPolynomialVariables@FR( y ) = 1 );;
gap>     Assert( 0, CountPolynomialVariables@FR( x * y ) = 2 );;
gap>     Assert( 0, CountPolynomialVariables@FR( x + y ) = 2 );;
#
#
#  @FR@Utils.Tests.TEST_DISTINCT_MONIC_FACTORS : 
#
gap> 
gap>     rng := PolynomialRing( ZmodnZ( 11 ), [ "x" ] );;
gap>     indeterminates := IndeterminatesOfPolynomialRing( rng );;
gap>     x := indeterminates[1];;
gap>     pol := 4 * (x - 3) ^ 10 * (3 * x - 2) ^ 3;;
gap>     result := DistinctMonicFactors@FR( pol );;
gap>     Assert( 0, result = [ x - 3, x - 8 ] );;
gap>     pol := 4 * x ^ 0;;
gap>     result := DistinctMonicFactors@FR( pol );;
gap>     Assert( 0, Size( result ) = 0 );;
#
#
#  @FR@Utils.Tests.TEST_PRODUCT_VALUE : 
#
gap> 
gap>     rng := PolynomialRing( ZmodnZ( 11 ), [ "x" ] );;
gap>     indeterminates := IndeterminatesOfPolynomialRing( rng );;
gap>     x := indeterminates[1];;
gap>     product := [ [ 2, 3 ] ];;
gap>     Assert( 0, 2 ^ 3 = PRODUCT_VALUE@FR( product ) );;
gap>     product := [ [ x - 3, 3 ] ];;
gap>     Assert( 0, (x - 3) ^ 3 = PRODUCT_VALUE@FR( product ) );;
gap>     product := [ [ x - 3, 3 ], [ x, 2 ] ];;
gap>     Assert( 0, (x - 3) ^ 3 * x ^ 2 = PRODUCT_VALUE@FR( product ) );;
gap>     product := [  ];;
gap>     Assert( 0, 1 = PRODUCT_VALUE@FR( product ) );;
#
#
#  @FR@Utils.Tests.TEST_UNIQUE_PRODUCT : 
#
gap> 
gap>     rng := PolynomialRing( ZmodnZ( 11 ), [ "x" ] );;
gap>     indeterminates := IndeterminatesOfPolynomialRing( rng );;
gap>     x := indeterminates[1];;
gap>     pol := (x - 3) ^ 3;;
gap>     result := UNIQUE_PRODUCT@FR( pol );;
gap>     Assert( 0, result = [ [ x - 3, 3 ] ] );;
gap>     pol := 3 * (x - 3) ^ 3;;
gap>     result := UNIQUE_PRODUCT@FR( pol );;
gap>     Assert( 0, result = [ [ x - 3, 3 ], [ One( rng ) * 3, 1 ] ] );;
gap>     pol := (x - 3) ^ 3 * x ^ 2;;
gap>     result := UNIQUE_PRODUCT@FR( pol );;
gap>     expectedProduct := [ [ x, 2 ], [ x - 3, 3 ] ];;
gap>     Assert( 0, expectedProduct = result );;
gap>     pol := x ^ 0;;
gap>     result := UNIQUE_PRODUCT@FR( pol );;
gap>     expectedProduct := [  ];;
gap>     Assert( 0, expectedProduct = result );;
gap>     pol := 5 * x ^ 0;;
gap>     result := UNIQUE_PRODUCT@FR( pol );;
gap>     expectedProduct := [ [ One( rng ) * 5, 1 ] ];;
gap>     Assert( 0, expectedProduct = result );;
#
#
#  @FR@Utils.Tests.TEST_SORT_POWERS_BY_EXPONENT : 
#
gap> 
gap>     factors := [ [ 3, 2 ], [ 3, 1 ], [ 4, 2 ], [ 3, 3 ] ];;
gap>     sortedFactors := SORT_POWERS_BY_EXPONENT@FR( factors );;
gap>     expectedResult := [ [ [ 3, 1 ] ], [ [ 3, 2 ], [ 4, 2 ] ], [ [ 3, 3 ] ] ];;
gap>     Assert( 0, expectedResult = sortedFactors );;
gap>     factors := [  ];;
gap>     sortedFactors := SORT_POWERS_BY_EXPONENT@FR( factors );;
gap>     expectedResult := [  ];;
gap>     Assert( 0, expectedResult = sortedFactors );;






#E utils.tst . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
