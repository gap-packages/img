#############################################################################
##
#W  padicLift.tst                  FR Package               Jakob Kroeker
##
#H  @(#)$Id$
##
#Y  Copyright (C) 2012,  Laurent Bartholdi
##
#############################################################################
##
##  This file tests the padicLift
##
#############################################################################

#  following lines generated with "str:= @PadicLift.CreateTestString(true); Print(str);"
#
#
# @PadicLift.Tests.TEST_LIFT_OPTIONS : 
#
gap> 
gap>     liftOptions := LiftOptions@FR(  );;
gap>     liftOptions.setMaxLiftDepth( 22 );;
gap>     Assert( 0, liftOptions.maxLiftDepth(  ) = 22 );;
gap>     liftOptions.setMaxLatticeDim( 3 );;
gap>     Assert( 0, liftOptions.maxLatticeDim(  ) = 3 );;
gap>     liftOptions.setVerboseLevel( 2 );;
gap>     Assert( 0, liftOptions.verboseLevel(  ) = 2 );;
gap>     liftOptions.setVerbosePairing( false );;
gap>     Assert( 0, liftOptions.verbosePairing(  ) = false );;
gap>     liftOptions.setInitialLatticeDim( 4 );;
gap>     Assert( 0, liftOptions.initialLatticeDim(  ) = 4 );;
gap>     liftOptions.setInitialLiftDepth( 0 );;
gap>     Assert( 0, liftOptions.initialLiftDepth(  ) = 0 );;
gap>     liftOptions.setMaxPairingTolerance( 0.1 );;
gap>     Assert( 0, liftOptions.maxPairingTolerance(  ) = 0.1 );;
gap>     CHECK_LIFT_OPTIONS@FR( liftOptions );;
#
#
# @PadicLift.Tests.TEST_LLL : 
#
gap> 
gap>     mat := [ [ 1, 2 ], [ 2, 1 ] ];;
gap>     lllResult := FPLLLReducedBasis( mat );;
gap>     Assert( 0, lllResult = [ [ 1, -1 ], [ 1, 2 ] ] );;
#
#
# @PadicLift.Tests.TEST_JENKINS_TRAUB_USAGE : 
#
gap> 
gap>     rng := PolynomialRing( Rationals, [ "x" ] );;
gap>     indeterminates := IndeterminatesOfPolynomialRing( rng );;
gap>     x := indeterminates[1];;
gap>     FZ1 := 33 * x ^ 3 + 19 * x ^ 2 - 81 * x - 4;;
gap>     roots := RootsByJenkinsTraub@FR( FZ1, 16 );;
gap>     roots := RootsByJenkinsTraub@FR( FZ1, 320 );;
gap>     roots := RootsByJenkinsTraub@FR( FZ1, 330 );;
gap>     rootCalculator := CreateJenkinsTraubWrapper@FR( 16 );;
gap>     roots := rootCalculator.computeRoots( FZ1 );;
gap>     roots := rootCalculator.computeRoots( FZ1 );;
gap>     roots := rootCalculator.computeRoots( FZ1 );;
#
#
# @PadicLift.Tests.TEST_LIFT_STEP_1@FR : 
#
gap> 
gap>     rng := PolynomialRing( Rationals, [ "x" ] );;
gap>     ind := IndeterminatesOfPolynomialRing( rng );;
gap>     x := ind[1];;
gap>     FZ := 33 * x ^ 3 + 19 * x ^ 2 - 81 * x - 4;;
gap>     ideal := Ideal( rng, [ FZ ] );;
gap>     jac := Jacobian@FR( [ FZ ], ind );;
gap>     solution := [ Z( 11 ) ^ 0 ];;
gap>     gens := GeneratorsOfTwoSidedIdeal( ideal );;
gap>     Assert( 0, IsZero( Value( FZ, ind, solution ) ) );;
gap>     Assert( 0, IsZero( EvalPolynomialTensor@FR( gens, ind, solution ) ) );;
gap>     solution := QuadraticLiftStep@FR( gens, jac, ind, solution );;
gap>     Assert( 0, IsZero( EvalPolynomialTensor@FR( gens, ind, solution ) ) );;
#
#
# @PadicLift.Tests.TEST_BLACKBOX_LIFT_STEP_1 : 
#
gap> 
gap>     rng := PolynomialRing( Rationals, [ "x" ] );;
gap>     ind := IndeterminatesOfPolynomialRing( rng );;
gap>     x := ind[1];;
gap>     FZ := 33 * x ^ 3 + 19 * x ^ 2 - 81 * x - 4;;
gap>     ideal := Ideal( rng, [ FZ ] );;
gap>     jac := Jacobian@FR( [ FZ ], ind );;
gap>     solution := [ Z( 11 ) ^ 0 ];;
gap>     gens := GeneratorsOfTwoSidedIdeal( ideal );;
gap>     Assert( 0, IsZero( Value( FZ, ind, solution ) ) );;
gap>     Assert( 0, IsZero( EvalPolynomialTensor@FR( gens, ind, solution ) ) );;
gap>     evalIdealGens := function ( point )
gap>           return EvalPolynomialTensor@FR( gens, ind, point );;
gap>       end;;
gap>     jacobianAt := function ( point )
gap>           return EvalPolynomialTensor@FR( jac, ind, point );;
gap>       end;;
gap>     solution := BlackBoxQuadraticLiftStep@FR( evalIdealGens, jacobianAt, solution );;
gap>     Assert( 0, IsZero( EvalPolynomialTensor@FR( gens, ind, solution ) ) );;
#
#
# @PadicLift.Tests.TEST_LIFT_STEP_2 : 
#
gap> 
gap>     problem := CREATE_RATIONAL_TEST_PROBLEM@FR(  );;
gap>     gens := GeneratorsOfTwoSidedIdeal( problem.ideal );;
gap>     jac := Jacobian@FR( gens, problem.indeterminates );;
gap>     Assert( 0, IsZero( EvalPolynomialTensor@FR( gens, problem.indeterminates, problem.solution ) ) );;
gap>     solution := QuadraticLiftStep@FR( gens, jac, problem.indeterminates, problem.solution );;
gap>     Assert( 0, IsZero( EvalPolynomialTensor@FR( gens, problem.indeterminates, solution ) ) );;
#
#
# @PadicLift.Tests.TEST_PADIC_LIFT : 
#
gap> 
gap>     problem := CREATE_RATIONAL_TEST_PROBLEM@FR(  );;
gap>     solution := PadicLift@FR( problem.ideal, problem.solution, 3 );;
gap>     gens := GeneratorsOfTwoSidedIdeal( problem.ideal );;
gap>     Assert( 0, IsZero( EvalPolynomialTensor@FR( gens, problem.indeterminates, solution ) ) );;
#
#
# @PadicLift.Tests.TEST_BLACKBOX_PADIC_LIFT : 
#
gap> 
gap>     problem := CREATE_RATIONAL_TEST_PROBLEM@FR(  );;
gap>     gens := GeneratorsOfTwoSidedIdeal( problem.ideal );;
gap>     evalIdealGens := function ( point )
gap>           return EvalPolynomialTensor@FR( gens, problem.indeterminates, point );;
gap>       end;;
gap>     jac := Jacobian@FR( gens, problem.indeterminates );;
gap>     jacobianAt := function ( point )
gap>           return EvalPolynomialTensor@FR( jac, problem.indeterminates, point );;
gap>       end;;
gap>     solution := BlackBoxPadicLift@FR( evalIdealGens, jacobianAt, problem.solution, 3 );;
gap>     Assert( 0, IsZero( evalIdealGens( solution ) ) );;
#
#
# @PadicLift.Tests.TEST_LLL_REDUCTION : 
#
gap> 
gap>     problem := CREATE_RATIONAL_TEST_PROBLEM@FR(  );;
gap>     liftResult := PadicLift@FR( problem.ideal, problem.solution, 3 );;
gap>     nextLiftResult := PadicLift@FR( problem.ideal, problem.solution, 4 );;
gap>     gens := GeneratorsOfTwoSidedIdeal( problem.ideal );;
gap>     Assert( 0, IsZero( EvalPolynomialTensor@FR( gens, problem.indeterminates, liftResult ) ) );;
gap>     Assert( 0, IsZero( EvalPolynomialTensor@FR( gens, problem.indeterminates, nextLiftResult ) ) );;
gap>     reductionOpts := LiftOptions@FR(  );;
gap>     LLL_REDUCTION_ATTEMPT@FR( problem.unknowns[1], problem.indeterminates, liftResult, nextLiftResult, reductionOpts );;
#
#
# @PadicLift.Tests.TEST_COMPUTE_MINIMAL_POLYNOMIAL : 
#
gap> 
gap>     problem := CREATE_RATIONAL_TEST_PROBLEM@FR(  );;
gap>     options := LiftOptions@FR(  );;
gap>     unknown := problem.indeterminates[1];;
gap>     minimalPolynomialVariable := Indeterminate( Rationals );;
gap>     liftAndLLLRes := ComputeMinimalPolynomialEx@FR( problem.ideal, problem.solution, unknown, minimalPolynomialVariable, options );;
gap>     unknown := problem.indeterminates[2];;
gap>     liftAndLLLRes := ComputeMinimalPolynomialEx@FR( problem.ideal, problem.solution, unknown, minimalPolynomialVariable, options );;
#
#
# @PadicLift.Tests.TEST_COMPUTE_MINIMAL_POLYNOMIALS : 
#
gap> 
gap>     liftAndLLLOptions := LiftOptions@FR(  );;
gap>     problem := CREATE_RATIONAL_TEST_PROBLEM@FR(  );;
gap>     x := problem.indeterminates[1];;
gap>     y := problem.indeterminates[2];;
gap>     liftAndLLLRes := ComputeMinimalPolynomials@FR( problem.ideal, problem.solution, problem.unknowns, liftAndLLLOptions );;
gap>     expectedUnknowns := [ [ x, -11 * x ^ 2 - 21 * x - 1 ], [ y, y - 1 ] ];;
gap>     expectedMergedLiftInfo := rec(
gap>         dataType := "LiftInfo",
gap>         maxLatticeDimension := 3,
gap>         maxLiftDepth := 3,
gap>         requiredLatticeDimension := 3 );;
gap>     Assert( 0, liftAndLLLRes.unknowns = expectedUnknowns );;
gap>     Assert( 0, liftAndLLLRes.mergedLiftInfo = expectedMergedLiftInfo );;
#
#
# @PadicLift.Tests.TEST_COMPATIBILITY_ROWS_VALID : 
#
gap> 
gap>     matrix := [ [ 0, 1 ], [ 2, 0 ], [ 0, 3 ] ];;
gap>     Assert( 0, COMPATIBILITY_ROWS_VALID@FR( matrix, false ) );;
gap>     Assert( 0, COMPATIBILITY_ROWS_VALID@FR( matrix, true ) );;
gap>     matrix := [ [ 2, 1 ], [ 2, 0 ], [ 0, 3 ] ];;
gap>     Assert( 0, COMPATIBILITY_ROWS_VALID@FR( matrix, false ) );;
gap>     matrix := [ [ 2, 1 ], [ 2, 0 ], [ 0, 3 ] ];;
gap>     Assert( 0, not COMPATIBILITY_ROWS_VALID@FR( matrix, true ) );;
gap>     matrix := [ [ 2, 1 ], [ 0, 0 ], [ 0, 3 ] ];;
gap>     Assert( 0, not COMPATIBILITY_ROWS_VALID@FR( matrix, true ) );;
gap>     Assert( 0, not COMPATIBILITY_ROWS_VALID@FR( matrix, false ) );;
#
#
# @PadicLift.Tests.TEST_IS_VALID_ROOT_COMPATIBILITY : 
#
gap> 
gap>     logger := function ( a, b )
gap>           return;;
gap>       end;;
gap>     matrix := [ [ 1, 2 ], [ 1, 4 ], [ 5, 6 ] ];;
gap>     Assert( 0, false = IS_VALID_ROOT_COMPATIBILITY@FR( matrix, 6, logger ) );;
gap>     matrix := [ [ 1, 2 ], [ 3, 4 ], [ 5, 6 ] ];;
gap>     Assert( 0, true = IS_VALID_ROOT_COMPATIBILITY@FR( matrix, 6, logger ) );;
gap>     Assert( 0, true = IS_VALID_ROOT_COMPATIBILITY@FR( matrix, 6, logger ) );;
gap>     matrix := [ [ 1, 0 ], [ 3, 0 ], [ 2, 0 ] ];;
gap>     Assert( 0, false = IS_VALID_ROOT_COMPATIBILITY@FR( matrix, 3, logger ) );;
gap>     matrix := [ [ 1, 0 ], [ 0, 3 ], [ 2, 0 ] ];;
gap>     Assert( 0, true = IS_VALID_ROOT_COMPATIBILITY@FR( matrix, 3, logger ) );;
#
#
# @PadicLift.Tests.TEST_COMPUTE_ROOT_COMPATIBILITY : 
#
gap> 
gap>     SetFloats( MPC, 1000 );;
gap>     firstPolRoots := [ 0.03, 34.0, 10.0 ];;
gap>     secondPolRoots := [ 5.03, 4.0, 1.0 ];;
gap>     combinedPolRoots := [ 4.03, 11.0, 39.02 ];;
gap>     operation := function ( a, b )
gap>           return a + b;;
gap>       end;;
gap>     opts := LiftOptions@FR(  );;
gap>     opts.setMaxPairingTolerance( 0.001 );;
gap>     compatibility := COMPUTE_HURWITZ_ROOT_COMPATIBILITY@FR( firstPolRoots, secondPolRoots, combinedPolRoots, operation, opts.maxPairingTolerance(  ), opts.logger );;
gap>     Assert( 0, compatibility = fail );;
gap>     opts := LiftOptions@FR(  );;
gap>     opts.setMaxPairingTolerance( 0.02 );;
gap>     opts.setVerbosePairing( false );;
gap>     compatibility := COMPUTE_HURWITZ_ROOT_COMPATIBILITY@FR( firstPolRoots, secondPolRoots, combinedPolRoots, operation, opts.maxPairingTolerance(  ), opts.logger );;
gap>     Assert( 0, compatibility = [ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ] );;
gap>     firstPolRoots := [ 4.0, 10.0 ];;
gap>     secondPolRoots := [ 5.0 ];;
gap>     combinedPolRoots := [ 9.0, 15.0 ];;
gap>     compatibility := ComputeRootCompatibilityEx@FR( firstPolRoots, secondPolRoots, combinedPolRoots, operation, opts.maxPairingTolerance(  ), opts.logger );;
gap>     Assert( 0, compatibility = [ [ 1 ], [ 2 ] ] );;
#
#
# @PadicLift.Tests.TEST_COMPUTE_APPROX_IDEAL_POINTS : 
#
gap> 
gap>     TestHelper := function ( problem )
gap>           local  opts, gens, result, errorTolerance, evaluation, evaluationAbs, max, root;;
gap>           opts := LiftOptions@FR(  );;
gap>           result := ComputeApproxIdealPoints@FR( problem.ideal, problem.solution, opts );;
gap>           gens := GeneratorsOfTwoSidedIdeal( problem.ideal );;
gap>           errorTolerance := 1.e-14;;
gap>           for root  in result.approxIdealElems  do
gap>               evaluation := EvalPolynomialTensor@FR( gens, problem.indeterminates, root );;
gap>               evaluationAbs := List( evaluation, function ( n )
gap>                       return AbsoluteValue( n );;
gap>                   end );;
gap>               max := Maximum( evaluationAbs );;
gap>               Assert( 0, max < errorTolerance );;
gap>           od;;
gap>           return;;
gap>       end;;
gap>     TestHelper( CREATE_RATIONAL_TEST_PROBLEM@FR(  ) );;
gap>     TestHelper( CREATE_SYMM_TEST_PROBLEM@FR(  ) );;
#
#
# @PadicLift.Tests.TEST_COMPUTE_HURWITZ_APPROX_IDEAL_POINT : 
#
gap> 
gap>     problem := CREATE_RATIONAL_TEST_PROBLEM@FR(  );;
gap>     opts := LiftOptions@FR(  );;
gap>     result := COMPUTE_APPROX_HURWITZ_IDEAL_POINTS@FR( problem.ideal, problem.solution, opts );;
gap>     gens := GeneratorsOfTwoSidedIdeal( problem.ideal );;
gap>     errorTolerance := 1.e-14;;
gap>     for root  in result.approxIdealElems  do
gap>         evaluation := EvalPolynomialTensor@FR( gens, problem.indeterminates, root );;
gap>         evaluationAbs := List( evaluation, function ( n )
gap>                 return AbsoluteValue( n );;
gap>             end );;
gap>         max := Maximum( evaluationAbs );;
gap>         Assert( 0, max < errorTolerance );;
gap>     od;;
#
#
# @PadicLift.Tests.TEST_COERCE_POLYNOMIAL_TO_COMPLEX_RING : 
#
gap> 
gap>     rng := PolynomialRing( Rationals, 1 );;
gap>     indeterminates := IndeterminatesOfPolynomialRing( rng );;
gap>     x := indeterminates[1];;
gap>     pol := x ^ 2 + 3;;
gap>     dstrng := PolynomialRing( MPC_PSEUDOFIELD, 1 );;
gap>     coercedPol := CoercePolynomialTensor@FR( pol, dstrng );;
gap>     dstInd := IndeterminatesOfPolynomialRing( dstrng );;
gap>     expectedResult := dstInd[1] ^ 2 + 3.0_c;;
gap>     Assert( 0, coercedPol = expectedResult );;



#E padicLift.tst . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
