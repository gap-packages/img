# is there a class of ideals ?
# todo: - add blackbox functionality  - partly done
# Q: does setFloats have global impact? 



InstallGlobalFunction( RootsByJenkinsTraub@FR ,
function ( polynomial, decimalPrecision)
    local bitPrecision, complexUnivariatePolynomialRing, coeffData, conversionFactor,
    controlDecimalPrecision, coeffDataCopy, fam, complexPol, pos ;

    if not IsUnivariatePolynomial(polynomial) then 
    	Error("JenkinsTraub: first parameter is not an univariate polynomial!" );
    fi;
    
    if not IsPosInt(decimalPrecision) then 
    	Error("JenkinsTraub: second parameter decimalPrecision is not a positive integer!" );
    fi;

    # determine bitPrecision from decimalPrecision
    ##################################   
      SetFloats( MPC, 100 );    

        #conversionFactor :=  3.32192809488736; 
        conversionFactor :=  Log(10.0_c)/Log(2.0_c);
    
        bitPrecision := decimalPrecision*conversionFactor; 
        # bitPrecision := convertFloatToIntSimple( bitPrecision );
        bitPrecision :=  Int(RealPart( bitPrecision) );
          
    
        conversionFactor :=  Log(2.0_c)/Log(10.0_c);
        controlDecimalPrecision := bitPrecision*conversionFactor;
        while controlDecimalPrecision<decimalPrecision*1.0 do
            bitPrecision := bitPrecision + 1;
            controlDecimalPrecision := bitPrecision*conversionFactor;
        od;
    ##################################   
        
    SetFloats( MPC, bitPrecision );    

    complexUnivariatePolynomialRing := PolynomialRing( MPC_PSEUDOFIELD,1 ); # 1 indeterminate
 
    fam := FamilyObj( One(complexUnivariatePolynomialRing) ); 
        

    # coerce polynomial to  complexUnivariatePolynomialRing.      
    if not FamilyObj(polynomial) = fam then
      complexPol := CoercePolynomialTensor@FR(polynomial,complexUnivariatePolynomialRing);
    fi;
    ##################################

    return RootsFloat(complexPol);;
end
);


InstallGlobalFunction( CreateJenkinsTraubWrapper@FR ,
function( decimalPrecision )

  local rootCalculator, bitPrecision, complexUnivariatePolynomialRing,  fam ;

    rootCalculator := rec();

    rootCalculator.getDecimalPrecision := function()
        return decimalPrecision;
    end;

     rootCalculator.decimalToBitPrecision := function( decimalPrecision )
         local localBitPrecision,  conversionFactor,
        controlDecimalPrecision ;

        SetFloats( MPC, 100 );    
    
        conversionFactor :=  Log(10.0_c)/Log(2.0_c);
    
        localBitPrecision := decimalPrecision*conversionFactor; 
        localBitPrecision := Int ( RealPart( localBitPrecision ) );
    
        conversionFactor :=  Log(2.0_c)/Log(10.0_c);
        controlDecimalPrecision := localBitPrecision*conversionFactor;
        while controlDecimalPrecision<decimalPrecision*1.0 do
            localBitPrecision := localBitPrecision + 1;
            controlDecimalPrecision := localBitPrecision*conversionFactor;
        od;
        return localBitPrecision;
    end;

    bitPrecision := rootCalculator.decimalToBitPrecision(decimalPrecision);

    rootCalculator.BitPrecision := function()
        return bitPrecision;
    end;


    SetFloats( MPC, bitPrecision );    

    complexUnivariatePolynomialRing := PolynomialRing( MPC_PSEUDOFIELD,1 ); # 1 indeterminate
    fam := FamilyObj( One(complexUnivariatePolynomialRing) ); 
    
    rootCalculator.getDstPolynomialFam := function()
        return fam;
    end;
    
    
    rootCalculator.getPolynomialRing :=function( )  
        return complexUnivariatePolynomialRing ;
    end;


    rootCalculator.computeRoots := function (polynomial)
        local complexPoly;
        complexPoly := polynomial;
        # only convert if not already converted. Todo: this should be checked in 'CoercePolynomialTensor'. 
        # Reason for this special case: coercion fails, if polynomial is already in complex ring.
        if not FamilyObj(One(complexPoly)) = fam then
            complexPoly := CoercePolynomialTensor@FR( polynomial, complexUnivariatePolynomialRing );
        fi;
        return RootsFloat(complexPoly);
    end;
    return Immutable(rootCalculator);
end
);



InstallGlobalFunction( QuadraticLiftStep@FR ,
function ( gens, jacobian, indeterminates, solutionApprox )
    local nextChar, higherSolutionApprox, JacobianAtSolution,  rightHandSide, correction, idealRing;
	Assert(0, IsZero( EvalPolynomialTensor@FR(gens, indeterminates, solutionApprox)) );

	if not IsZero( EvalPolynomialTensor@FR(gens, indeterminates, solutionApprox) ) then
		Error("solution does not belong to ideal");
	return [];
	fi;
 
    # improve padic approximation
        nextChar := ( Characteristic(solutionApprox) )^2;
        higherSolutionApprox := PromoteScalarTensor@FR( solutionApprox, ZmodnZ( nextChar ) ) ;
        JacobianAtSolution := EvalPolynomialTensor@FR( jacobian, indeterminates, higherSolutionApprox);
    
        if not Size(indeterminates) = Rank(PromoteScalarTensor@FR(JacobianAtSolution,Integers)) then
            Error("Jacobian is not invertible");
            return [];
        fi;
            
        rightHandSide := EvalPolynomialTensor@FR(gens, indeterminates, higherSolutionApprox);
    
        correction := -(JacobianAtSolution^-1) *rightHandSide;
        higherSolutionApprox := higherSolutionApprox + correction;

    # result check
        Assert(0, IsZero( EvalPolynomialTensor@FR(gens,indeterminates,higherSolutionApprox)) );

    return higherSolutionApprox;
end
);


# blackbox : ideal can check if point belongs to ideal; 
# -for a point a jacobian can be computed.
InstallGlobalFunction( BlackBoxQuadraticLiftStep@FR ,
function ( evalIdealGens, computeJacobianAt, solutionApprox )
    local nextChar, higherSolutionApprox, JacobianAtSolution,  rightHandSide, correction;
   
    # parameter consistency check
    
        if not IsZero( evalIdealGens(  solutionApprox) ) then
            Error("solution does not belong to ideal");
            return [];
        fi;
 
    # improve padic approximation
        nextChar := ( Characteristic(solutionApprox) )^2;
        higherSolutionApprox := PromoteScalarTensor@FR( solutionApprox, ZmodnZ( nextChar ) ) ;
        JacobianAtSolution := computeJacobianAt( higherSolutionApprox );
       
        if not Maximum(Size(JacobianAtSolution),Size(JacobianAtSolution[1])) = Rank(PromoteScalarTensor@FR(JacobianAtSolution,Integers)) then
            Error("Jacobian is not invertible");
            return [];
        fi;
            
        rightHandSide := evalIdealGens(  higherSolutionApprox );
    
        correction := -(JacobianAtSolution^-1) *rightHandSide;
        higherSolutionApprox := higherSolutionApprox + correction;

    # result check
        Assert(0, IsZero( evalIdealGens(  higherSolutionApprox) )  );

    return higherSolutionApprox;
end
);


InstallGlobalFunction( PadicLift@FR ,
function( ideal,  solutionPoint , numLiftDepth )
    
    local charSolution, laring , indeterminates , gens , JacobianOfIdeal, currLiftDepth, localSolution;
    charSolution := Characteristic(solutionPoint);
    Assert( 0, not IsZero( charSolution ));
    Assert( 0, IsOne( Size( Set( Factors(charSolution) ))) ); # 

    #Assert(0, IsIdeal(ideal) ); # todo: how to check ? 
    Assert( 0, IsZero( Characteristic (ideal)) );
    laring := LeftActingRingOfIdeal( ideal );
  
    Assert(0, IsPolynomialRing(laring) );
    indeterminates := IndeterminatesOfPolynomialRing(laring);
  
    gens := GeneratorsOfTwoSidedIdeal(ideal);   
    Assert(0, IsZero(EvalPolynomialTensor@FR(gens, indeterminates, solutionPoint)));
    
    JacobianOfIdeal := Jacobian@FR( gens, indeterminates);

    currLiftDepth := 0;
    localSolution := solutionPoint;
    while currLiftDepth < numLiftDepth do 
      localSolution := QuadraticLiftStep@FR( gens, JacobianOfIdeal,indeterminates, localSolution);
      currLiftDepth := currLiftDepth + 1;
    od;
    return localSolution;
end
);



InstallGlobalFunction( BlackBoxPadicLift@FR ,
 function( evaluateIdealGens, jacobianAt,  solutionPoint , numLiftDepth)
    
    local charSolution, currLiftDepth, localSolution;
    charSolution := Characteristic(solutionPoint);
    Assert( 0,not IsZero( charSolution ));
    Assert( 0, IsOne( Size( Set( Factors(charSolution) ))) ); # 
   
    Assert(0, IsZero(evaluateIdealGens( solutionPoint)) );

    currLiftDepth := 0;
    localSolution := solutionPoint;
    while currLiftDepth < numLiftDepth do 
      localSolution := BlackBoxQuadraticLiftStep@FR( evaluateIdealGens, jacobianAt, localSolution);
      currLiftDepth := currLiftDepth + 1;
    od;
    return localSolution;
end
);



InstallGlobalFunction( CHECK_LIFT_OPTIONS@FR ,
function(liftOptions)
    Assert(0, liftOptions.verbose()=true or liftOptions.verbose()=false);
    Assert(0, liftOptions.verbosePairing()=true or liftOptions.verbosePairing()=false);
    Assert(0, liftOptions.verboseLevel() in Integers);
    Assert(0,  liftOptions.decimalPrecision() in PositiveIntegers);
    Assert(0,  liftOptions.minColumnNormDistanceFactor() in PositiveIntegers);
    Assert(0,  liftOptions.initialLiftDepth() in NonnegativeIntegers);
    Assert(0,  liftOptions.initialLatticeDim() in PositiveIntegers);
    Assert(0,  liftOptions.maxLiftDepth() in PositiveIntegers or liftOptions.maxLiftDepth()=infinity );
    Assert(0,  liftOptions.maxLatticeDim() in PositiveIntegers or liftOptions.maxLatticeDim()=infinity );
    Assert(0,  liftOptions.latticeDimIncrementFkt(5) > 5 );
    Assert(0,  IsFloat( liftOptions.maxPairingTolerance() ) and AbsoluteValue( liftOptions.maxPairingTolerance() ) =  liftOptions.maxPairingTolerance()  );
end
);

InstallGlobalFunction(  CREATE_EMPTY_LOGGER_FKT@FR, 
function()
	return  function (level, message) end;
end
);

	

# optional: maybe split lift options and pairing options...
InstallGlobalFunction ( CREATE_LIFT_OPTIONS@FR ,
 function(optionData)
	local privateData, liftOptions;

	privateData :=  optionData ;
	
	liftOptions := rec( );

	liftOptions.decimalPrecision := function() return  privateData.rootCalculator.getDecimalPrecision(); end;
	liftOptions.setDecimalPrecision := function(precision)  
	    Assert(0, precision in PositiveIntegers);
        privateData.rootCalculator :=	privateData.rootCalculatorConstructor( precision );
        # todo: improvement: RootCalculator supports 'setDecimalPrecision'.
		#Error(" please call setRootCalculator instead: e.g.  setRootCalculator( CreateJenkinsTraubWrapper@FR( decimalPrecision );");
		
	end;
	
	liftOptions.rootCalculator := function() return  privateData.rootCalculator; end;
	liftOptions.setRootCalculator := function(rootCalculator)  
		local rnames;
		rnames := RecNames(rootCalculator);
		if IsRecord(rootCalculator) and 
 		   "computeRoots" in rnames and
 		   "decimalPrecision" in rnames 
 		then 
			privateData.rootCalculator := rootCalculator;
		else
			Error("rootCalculator does not match interface (\"ComputeRoots\"(polynomial), \"decimalPrecision\"() ");
		fi;
	end;
	
	
	liftOptions.maxLiftDepth := function() return  privateData.maxLiftDepth; end;
	liftOptions.setMaxLiftDepth := function(liftDepth)  
		if not liftDepth in NonnegativeIntegers and not liftDepth=infinity then 
			Error(" maxLiftDepth has to be a nonnegative int or infinity");	
		fi;	
		privateData.maxLiftDepth := liftDepth;
	end;
	
	
	liftOptions.maxLatticeDim := function() return  privateData.maxLatticeDim; end;
	liftOptions.setMaxLatticeDim := function(maxLatticeDim)  
		if not  maxLatticeDim in PositiveIntegers and not maxLatticeDim=infinity then 
			Error(" maxLiftDepth has to be a positive int or infinity");	
		fi;	
		privateData.maxLatticeDim := maxLatticeDim;
	end;
	
	
	liftOptions.verbose := function() return  privateData.verbose; end;
	liftOptions.setVerbose := function(verbose)  
		if not verbose=true and not verbose=false then 
			Error(" setVerbose to true or to false ");	
		fi;	
		privateData.verbose := verbose;
	end;
	
	
	liftOptions.verbosePairing := function() return  privateData.verbosePairing; end;
	liftOptions.setVerbosePairing := function(verbosePairing)  
		if not verbosePairing=true and not verbosePairing=false then 
			Error(" setVerbosePairing to true or to false ");	
		fi;	
		privateData.verbosePairing := verbosePairing;
	end;
	
	
	liftOptions.verboseLevel := function() return  privateData.verboseLevel; end;
	liftOptions.setVerboseLevel := function(level)  
		if not level in NonnegativeIntegers then 
			Error(" verbose level not a nonnegative integer ");	
		fi;	
		if level>0 then
			privateData.verbose := true;
		fi;
		privateData.verboseLevel := level;
	end;
	
	
	liftOptions.minColumnNormDistanceFactor := function() return  privateData.minColumnNormDistanceFactor; end;
	liftOptions.setMinColumnNormDistanceFactor := function(factor)  
		if not factor>=1 then 
			Error(" expected min column norm distance factor >=1 ");	
		fi;	
		privateData.minColumnNormDistanceFactor := factor;
	end;
	
	
	liftOptions.initialLiftDepth := function() return  privateData.initialLiftDepth; end;
	liftOptions.setInitialLiftDepth := function(depth)  
		if not depth in NonnegativeIntegers then 
			Error(" initial lift depth not an integer ");	
		fi;	
		privateData.initialLiftDepth := depth;
	end;
	
	
	liftOptions.initialLatticeDim := function() return  privateData.initialLatticeDim; end;
	liftOptions.setInitialLatticeDim := function(latticeDim)  
		if not latticeDim in PositiveIntegers then 
			Error(" initial lattice dimension not a positive integer ");	
		fi;	
		privateData.initialLatticeDim := latticeDim;
	end;
	
	
	liftOptions.maxPairingTolerance := function() return  privateData.maxPairingTolerance; end;
	liftOptions.setMaxPairingTolerance := function( pairingTolerance )  
		if not  IsFloat(pairingTolerance) or AbsoluteValue(pairingTolerance)<>pairingTolerance or IsZero(pairingTolerance) then 
			Error(" root pairing tolerance has to be a positive floating number ");	
		fi;	
		privateData.maxPairingTolerance :=pairingTolerance;
	end;
	
	
	liftOptions.latticeDimIncrementFkt := function(val) return  privateData.latticeDimIncrementFkt(val); end;
	 liftOptions.setLatticeDimIncrementFkt := function(incrementFunction)
	 	if not IsFunction(incrementFunction) then 
	 		Error("set latticeDimIncrementFkt: expected a function accepting an integer");
	 	fi;
	 	privateData.latticeDimIncrementFkt := incrementFunction;
	 end;     
		                            
	
	liftOptions.clone 	:= function()
		local newLiftOptions;
		newLiftOptions := CREATE_LIFT_OPTIONS@FR( ShallowCopy( privateData )  ) ;
		return newLiftOptions;
	end;

	liftOptions.print := function()
		local name,recNames;
		Info(InfoFR,1 ,"LiftOptions object: \n");
		for name in RecNames(privateData) do
			if IsRecord( privateData.(name) ) then 
				recNames :=  ShallowCopy(String( RecNames(privateData.(name)) ));
				Assert(0, recNames [1]='[');
				Assert(0, recNames [Size(recNames)]=']');
				recNames[1] := '(';
				recNames[Size(recNames)] := ')';
				Info(InfoFR,1, Concatenation("\t", name, " \t := rec", recNames, ";\n" ));
			else
				Info(InfoFR,1, Concatenation("\t", name, " \t := ", String(privateData.(name)), ";\n" ));
			fi;
		od;
			Info(InfoFR,1 ,"end; \n");
	end;
	
	liftOptions.Setters := function()
		local name,replacedName;
		Info(InfoFR,1,"Set functions:\n");
		 for name in RecNames(liftOptions) do
			 
			 replacedName := ReplacedString(name,"set","");
			 
			 if  replacedName<>name  then 
			 	Info(InfoFR,1, Concatenation("\t", name, "(); \n" ));
			 fi;
		od;
	 end;
	 
 	liftOptions.Getters := function()
		local name,replacedName;
		Info(InfoFR,1,"Get functions:\n");
		 for name in RecNames(liftOptions) do
			 replacedName := ReplacedString(name,"set","");
			 if  Size(replacedName)=Size(name)  then 
			 	Info(InfoFR,1, Concatenation("\t", name, " (); \n" ));
			 fi;
		od;
	 end;
	 
	privateData.logger := function (opts, level, message)
	if opts.verbose() then
	    if  opts.verboseLevel() >= level then
	    if level<1 then level:=1; fi;
		Info(InfoFR,level,message);    
		Info(InfoFR,level,"\n");
	    fi;
	fi;
	end;
	
	liftOptions.logger := function (level, message)
	 	privateData.logger( liftOptions, level, message);
	end;
	
	liftOptions.dataType:="LiftOptions";
	
	CHECK_LIFT_OPTIONS@FR( liftOptions );
	
	liftOptions := Immutable( liftOptions );
	
	

	return liftOptions;
end
);


# todo: improve interface: optionsdata knows, which root calculator to use. 
InstallGlobalFunction( LiftOptions@FR , function()
local optionData, objectifiedData, privateData, liftOptions,getPrivateData;


	optionData:= rec( );
	 optionData.latticeDimIncrementFkt := function(latticeDim)
                                        return latticeDim+1;
		                            end;
	optionData.maxLiftDepth := infinity;
	optionData.maxLatticeDim := infinity;
	optionData.verbose := false;

	optionData.minColumnNormDistanceFactor := 100;
	optionData.initialLatticeDim := 1;
	optionData.initialLiftDepth := 0;

    optionData.rootCalculatorConstructor := CreateJenkinsTraubWrapper@FR ;
	optionData.verbosePairing := false;
	optionData.rootCalculator := CreateJenkinsTraubWrapper@FR( 16 );
	optionData.maxPairingTolerance := 0.1;
	optionData.verboseLevel := 0;
	
    	return CREATE_LIFT_OPTIONS@FR( optionData );
end
);





InstallMethod( IsLiftOptions, "", [IsRecord],
function(record)
	if not "type" in RecNames(record) then
		return false;
	fi;

	return record.dataType="LiftOptions";
end
);

# sandbox: 
# @PadicLift.Tests.TEST_LIFT_OPTIONS := TEST_LIFT_OPTIONS@FR;
InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_LIFT_OPTIONS", 
#DeclareGlobalFunction( "@PadicLift\.Tests\.TEST_LIFT_OPTIONS"); 
#InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_LIFT_OPTIONS", @PadicLift\.Tests\.TEST_LIFT_OPTIONS);
#InstallGlobalFunction( @PadicLift\.Tests\.TEST_LIFT_OPTIONS, 
 function()
	local liftOptions;
	liftOptions :=  LiftOptions@FR();

    #bla := Immutable(rec());
    #bla.bla := 5;
    #Error("bla");
    
	liftOptions.setMaxLiftDepth(22);
	Assert(0, liftOptions.maxLiftDepth()=22);
			
	liftOptions.setMaxLatticeDim(3);
	Assert(0, liftOptions.maxLatticeDim() = 3 );
	
	liftOptions.setVerboseLevel(2);
	Assert(0, liftOptions.verboseLevel()=2);
	
	liftOptions.setVerbosePairing(false);
	Assert(0, liftOptions.verbosePairing() = false );
	
	liftOptions.setInitialLatticeDim(4);
	Assert(0, liftOptions.initialLatticeDim() = 4 );
	
	liftOptions.setInitialLiftDepth(0);
	Assert(0, liftOptions.initialLiftDepth() = 0 );
	
	liftOptions.setMaxPairingTolerance(0.1);
	Assert(0, liftOptions.maxPairingTolerance() = 0.1 );
	
	CHECK_LIFT_OPTIONS@FR( liftOptions );
		
end);	



# input is a Integer Matrix where the rows form the Basis of the lattice.
# DeclareGlobalFunction("ROW_NORMS@FR");
InstallGlobalFunction( ROW_NORMS@FR ,
function(mat)
    local MM;
    if not IsMatrix(mat) then 
    	Error("ROW_NORMS@FR: parameter is not a matrix!");
    fi;
    
    MM  := PromoteScalarTensor@FR( mat, Rationals);
    MM := MM*TransposedMat(MM);
    return  List( [1..Size(MM)], n->MM[n][n] );
end
);

# todo: why did you copy code? => because TransposedMat could became costly for big LLL Matrices.
# DeclareGlobalFunction("COLUMN_NORMS@FR");
InstallGlobalFunction( COLUMN_NORMS@FR ,
function(mat)
    local MM;
    if not IsMatrix(mat) then 
    	Error("ColumnNorms: parameter is not a matrix!");
    fi;

    MM  := PromoteScalarTensor@FR( mat, Rationals);
    MM := TransposedMat(MM)*MM;
    return  List( [1..Size(MM)], n->MM[n][n] );
end
);



# DeclareGlobalFunction("NORMALIZED_ROW_NORMS@FR");
InstallGlobalFunction(  NORMALIZED_ROW_NORMS@FR ,
function (mat)
    local rowNormlist, minM, maxM, normalizedRowNormList, result,pos;
    if not IsMatrix(mat) then 
    	Error(" NORMALIZED_ROW_NORMS: parameter is not a matrix!");
    fi;
    
    rowNormlist := ROW_NORMS@FR( mat );
    minM := Minimum( rowNormlist );
    maxM := Maximum( rowNormlist );
    normalizedRowNormList :=  List( [ 1..Size(rowNormlist) ], pos-> rowNormlist[pos]*RealPart(1.0)/minM );
    #normalizedRowNormList :=  List( [ 1..Size(rowNormlist) ], pos-> rowNormlist[pos]/minM );
    result := rec();
    result.unchanged := rowNormlist;
    result.normalized := normalizedRowNormList;
    result.max := maxM;
    result.min := minM;
    result.dataType := "NormalizedRowNorms";
    return Immutable(result);  
end
);




InstallGlobalFunction( CREATE_LIFT_INFO@FR ,
function( maxLiftDepth, maxLatticeDimension, requiredLatticeDimension, minLiftDepth )
    local liftInfo;
    liftInfo := rec();
    liftInfo.dataType := "LiftInfo";
    liftInfo.minLiftDepth := minLiftDepth;
    liftInfo.maxLiftDepth := maxLiftDepth;
    liftInfo.maxLatticeDimension := maxLatticeDimension;
    liftInfo.requiredLatticeDimension := requiredLatticeDimension;
    return Immutable(liftInfo);
end
);


InstallGlobalFunction( MERGE_LIFT_INFO@FR ,
function( liftInfo1, liftInfo2 )
    local minLiftDepth,maxLiftDepth, maxLatticeDimension, requiredLatticeDimension;

    maxLiftDepth := Maximum ( liftInfo1.maxLiftDepth,liftInfo2.maxLiftDepth );
    minLiftDepth := Minimum ( liftInfo1.minLiftDepth,liftInfo2.minLiftDepth );
    Info(InfoFR,2, Concatenation ( "minLiftDepth: = ", String(minLiftDepth) , "\n" ) );
    maxLatticeDimension := Maximum ( liftInfo1.maxLatticeDimension, liftInfo2.maxLatticeDimension );
    requiredLatticeDimension := Null@FR;
    if ( not liftInfo1.requiredLatticeDimension=Null@FR  and not liftInfo2.requiredLatticeDimension=Null@FR ) then 
        requiredLatticeDimension :=  Maximum ( liftInfo1.requiredLatticeDimension, liftInfo2.requiredLatticeDimension );
    fi;
    return CREATE_LIFT_INFO@FR( maxLiftDepth, maxLatticeDimension, requiredLatticeDimension,minLiftDepth );
end
);



InstallGlobalFunction( LLL_INPUT_FROM_LIFT@FR ,
function( unknown, indeterminates, liftResult, currentLatticeDim )
    local liftResultOverIntegers, M, sM, idx, result;

    liftResultOverIntegers := PromoteScalarTensor@FR( liftResult, Rationals ); 
    M :=  EvalPolynomialTensor@FR( [ List( [0..currentLatticeDim-1], i->unknown^i) ], indeterminates, liftResultOverIntegers );
    Append( M[1], [ Characteristic(liftResult) ] );
    
    sM := List( [1..currentLatticeDim+1] , n-> List( [1..currentLatticeDim],l->0));
    
    # write kernel(M) : (each 'sM' column is a kernel element)
    for idx in [1..Size(M[1])-1] do
        sM[1][idx ] := -M[1][idx+1];
        sM[idx+1][idx] := M[1][1];
    od;
    # remove last 'sM' row and transpose the result.
    result := TransposedMat( List( [1..Size(sM)-1], n->sM[n] ) );
   
    return  PromoteScalarTensor@FR( result, Rationals );
end
);





# try to find for a given lift the minimal polynomial in variable 'unknown' by guessing its degree (heurustic method)
# latticeBasisNormList is evaluated by the heuristic method
InstallGlobalFunction( LLL_REDUCTION_ATTEMPT@FR ,
function (unknown, indeterminates, liftResult, nextLiftResult, reductionOpts ) 
     
   local reducedLiftResult, currentLatticeDim, lastColumnNormMin, LLLInput, bvec, basisNormRecord, nextLiftResultOverInts;
  

    reducedLiftResult := rec();
    reducedLiftResult.foundMinPolyCandidate := false;
    reducedLiftResult.latticeBasisNormList := [];
    reducedLiftResult.latticeBasis := Null@FR;
    reducedLiftResult.minPolynomial := Null@FR;
    reducedLiftResult.liftInfo := Null@FR;
    reducedLiftResult.currentLatticeDim   :=-1;
    reducedLiftResult.dataType := "ReducedPadicLiftResult";
 
    currentLatticeDim  := reductionOpts.initialLatticeDim();
    
    lastColumnNormMin := -1;
    nextLiftResultOverInts :=  PromoteScalarTensor@FR(nextLiftResult,Integers);

    while currentLatticeDim <= reductionOpts.maxLatticeDim()  do 
        reductionOpts.logger(1, Concatenation("# currentLatticeDim: ", String(currentLatticeDim) ) );
        LLLInput := LLL_INPUT_FROM_LIFT@FR(unknown, indeterminates, liftResult, currentLatticeDim );
        reducedLiftResult.latticeBasis := FPLLLReducedBasis( LLLInput );
     
        # test, if a solution have been found in this step (reducedLiftResult.foundMinPolyCandidate):
            bvec :=   EvalPolynomialTensor@FR(  List([0..currentLatticeDim-1],n->unknown^n ),  indeterminates, nextLiftResultOverInts ) ;
    
            basisNormRecord := NORMALIZED_ROW_NORMS@FR( reducedLiftResult.latticeBasis );
            if  IsZero( PromoteScalarTensor@FR( bvec*( reducedLiftResult.latticeBasis [1] ), nextLiftResult[1] ) ) and 
                EuclideanQuotient( basisNormRecord.max, basisNormRecord.min )> reductionOpts.minColumnNormDistanceFactor()  then 
                 reducedLiftResult.foundMinPolyCandidate := true;
            fi;
        #TODO: sometimes first condition ( "IsZero (PromoteScalarTensor@FR( bvec*firstBasisRow, nextLiftResult[1] ) )" ) passes, 
        #      but we do not have a solution. Due to HC if we will use a higher lift, this could be detected at the end.
        #
        #TODO: instead of hardcoding, parametrize lllReduction with a stop condition for increasing lattice dimension: 
        #       is it sufficient for a generic stop condition to pass as input previous and current latticeBasis ?
        
        Append (reducedLiftResult.latticeBasisNormList , [ basisNormRecord ]) ;
        reductionOpts.logger( 2, Concatenation("column norms: ", String( basisNormRecord.normalized ) ));
        if   lastColumnNormMin=basisNormRecord.min  then 
              reductionOpts.logger(1, " lattice dimension increase: stop condition triggered. ");
               break;
        fi;

       if      reducedLiftResult.foundMinPolyCandidate  then  
                reductionOpts.logger(1, Concatenation("found minpoly candidate; lattice dimension: ", String(currentLatticeDim) ) );
               break;
        fi;

        lastColumnNormMin := basisNormRecord.min;
        Assert(0,  (reductionOpts.latticeDimIncrementFkt ( currentLatticeDim )) > currentLatticeDim or currentLatticeDim = infinity );
        currentLatticeDim := reductionOpts.latticeDimIncrementFkt( currentLatticeDim );
    od;
    # todo: Typ für Rückgabe einfuehren.
    reducedLiftResult.currentLatticeDim   := currentLatticeDim;
    reductionOpts.logger( 3, Concatenation("LLLReductionAttempt result:", String(reducedLiftResult)));     
    return reducedLiftResult;
end
);



InstallGlobalFunction( LATTICE_BASIS_TO_POLYNOMIAL@FR,
function (latticeBasis, variable)
    local localVar, nrows, pol;
    localVar := Null@FR;
    if variable = Null@FR then 
         localVar := Indeterminate(Rationals);
    else
       # TODO: ensure that variable is an indeterminate  ( either of rationals or of Integers ) 
       localVar :=variable;    
    fi;
    nrows := Size( latticeBasis );
    pol :=  List( [0..nrows-1] , exp->localVar^exp) *( PromoteScalarTensor@FR( latticeBasis[1], localVar) );
    return pol;
end
);


# PadicLift.ComputeMinimalPolynomal
# optional: it is thinkable that opts.maxLiftDepth() is depending on the unknown (if there is some apriori knowledge)
InstallGlobalFunction(  ComputeMinimalPolynomialEx@FR ,
function( ideal, solution, unknown, minimalPolynomialVariable, liftOptions  )

    local  gens, jacobianOfIdeal,  currLiftDepth, indeterminates, liftResult, nextLiftResult, 
    localLiftOptions,  minimalPolynomialCandidateFactors, idealRing, reducedLiftResult ;
 
    CHECK_LIFT_OPTIONS@FR (liftOptions);

    idealRing := LeftActingRingOfIdeal(ideal);
    Assert(0, IsPolynomialRing(idealRing) );
    Assert(0, idealRing = RightActingRingOfIdeal(ideal) ) ;
    indeterminates := IndeterminatesOfPolynomialRing(idealRing);
    Assert(0, unknown in idealRing);

    liftOptions.logger(2, "ComputeMinimalPolynomial@FR" );
    
    # assert( idealRing === ring unknown); TODO: how to check?

    gens := GeneratorsOfTwoSidedIdeal( ideal );    
    Assert(0, IsZero( EvalPolynomialTensor@FR( gens, indeterminates, solution ) ) );

    jacobianOfIdeal := Jacobian@FR ( gens, indeterminates) ;

   
    reducedLiftResult := Null@FR;

    currLiftDepth := 0;
    liftResult :=  solution ;
    
    nextLiftResult :=    QuadraticLiftStep@FR( gens,  jacobianOfIdeal, indeterminates,  liftResult);

    # increase lift depth and perform LLL until a solution is found or maxLiftDepth is reached.
    while currLiftDepth <= liftOptions.maxLiftDepth()  do 
        liftOptions.logger(1, Concatenation("#\n # currLiftDepth: ", String(currLiftDepth) ));
        
        # perform LLL only if (currLiftDepth >= startingLiftDepth ). 
        # The condition is useful in case minimalLiftDepth (=startingLiftDepth ) is known (e.g. from similar previous computations )
        if ( currLiftDepth >= liftOptions.initialLiftDepth() ) then   
        
            reducedLiftResult := LLL_REDUCTION_ATTEMPT@FR(   unknown, indeterminates, liftResult,  nextLiftResult, liftOptions );
            
            if  reducedLiftResult.foundMinPolyCandidate  then  
                liftOptions.logger(1, Concatenation("#FinalLiftDepth: " ,String (currLiftDepth) ) );
                break;
            fi;
        fi;
        currLiftDepth := currLiftDepth+1;    
        liftResult := nextLiftResult;
        nextLiftResult :=     QuadraticLiftStep@FR(  gens,  jacobianOfIdeal, indeterminates, liftResult);
    od;
  
    if reducedLiftResult=Null@FR or not reducedLiftResult.foundMinPolyCandidate   then  
        Info(InfoFR,1, "failed to compute minimal polynomial");
        return fail;
    fi;
  
    reducedLiftResult.minPolynomial :=  LATTICE_BASIS_TO_POLYNOMIAL@FR( reducedLiftResult.latticeBasis, minimalPolynomialVariable );

    liftOptions.logger(1, Concatenation("---------------polynomial candidate degree: ", String(Degree(reducedLiftResult.minPolynomial)))  );
    minimalPolynomialCandidateFactors :=   Factors( reducedLiftResult.minPolynomial) ;

    if ( Size( minimalPolynomialCandidateFactors) >1) then
        liftOptions.logger(1, "----------------lattice dimension too big: reducing lattice dimension ");
       localLiftOptions :=  liftOptions.clone();
       localLiftOptions.setInitialLatticeDim( localLiftOptions.initialLatticeDim() - Size(minimalPolynomialCandidateFactors )+1 );
       return ComputeMinimalPolynomialEx@FR( ideal,  solution, unknown, minimalPolynomialVariable, localLiftOptions );
    fi;

    reducedLiftResult.unknown := unknown;
    reducedLiftResult.liftInfo := CREATE_LIFT_INFO@FR( currLiftDepth, reducedLiftResult.currentLatticeDim, (Degree (reducedLiftResult.minPolynomial) + 1),currLiftDepth );
   
    return Immutable(reducedLiftResult);
end
);


InstallGlobalFunction(  ComputeMinimalPolynomial@FR ,
function( ideal, solution, unknown,  liftAndLLLOptions)
	return ComputeMinimalPolynomial@FR( ideal, solution, unknown, unknown,   liftAndLLLOptions);
end
);


InstallGlobalFunction( ComputeMinimalPolynomials@FR ,
function( solutionIdeal,  solutionPoint, unknowns,  computeOptions )

    local unknown, liftResult, minimalPolynomialsData, mergedLiftInfo, 
          minPolVar, optsCopy, idealRing, indeterminates, unknownIdx;

    Assert(0,  LeftActingRingOfIdeal (solutionIdeal)=RightActingRingOfIdeal (solutionIdeal) );
    idealRing  :=  LeftActingRingOfIdeal (solutionIdeal);
    indeterminates := IndeterminatesOfPolynomialRing(idealRing);

    CHECK_LIFT_OPTIONS@FR (computeOptions);
  
    Assert(0, Characteristic(solutionPoint)>0 );
  
    mergedLiftInfo := CREATE_LIFT_INFO@FR(0,0,0,0);
    
    mergedLiftInfo := Null@FR;

    minimalPolynomialsData := rec(); 
    minimalPolynomialsData.dataType:= "PadicLift.MinimalPolynomials";
    minimalPolynomialsData.unknowns := [] ; # TODO maybe wanna to use a Hashtable in unknowns.
    minimalPolynomialsData.liftInfo := [] ;

    for unknownIdx in [1..Size(unknowns)] do
        unknown:=unknowns[unknownIdx];
      
        Info(InfoFR,2, Concatenation("------------------lifting variable ", String(unknownIdx),"(",String(Size(unknowns)),") -----------------------") );
        if Size(ExtRepPolynomialRatFun(unknown))=2 then
            minPolVar := unknown; # use unknown as variable for minimal polynomial.
        else
            # unknown variable is composed and cannot be used as variable for minimal polynomial.
            minPolVar :=  Indeterminate( Rationals ) ;; #maybe Integers are sufficient.
        fi;           
    
        # heuristic: adjust lift options. TODO: parametrise 'ComputeMinimalPolynomials' with heuristic.
            optsCopy :=  computeOptions.clone();
            
           if not mergedLiftInfo = Null@FR then
            if optsCopy.initialLiftDepth() < mergedLiftInfo.maxLiftDepth then 
                optsCopy.setInitialLiftDepth( mergedLiftInfo.maxLiftDepth );
            fi;
    
            if optsCopy.initialLatticeDim()  < mergedLiftInfo.requiredLatticeDimension then 
                optsCopy.setInitialLatticeDim ( mergedLiftInfo.requiredLatticeDimension);
            fi;
          fi;
        
        liftResult := ComputeMinimalPolynomialEx@FR( solutionIdeal,  solutionPoint, unknown, minPolVar, optsCopy );
        if liftResult=fail then
	        return fail;	
        fi;
          if not mergedLiftInfo = Null@FR then
             mergedLiftInfo := MERGE_LIFT_INFO@FR(  mergedLiftInfo, liftResult.liftInfo );
         else
                 mergedLiftInfo :=  liftResult.liftInfo ;
         fi;
        Append( minimalPolynomialsData.unknowns , [ [ unknown, liftResult.minPolynomial ] ] );
        Append( minimalPolynomialsData.liftInfo , [ liftResult.liftInfo ] );
    
    od;
    minimalPolynomialsData.mergedLiftInfo := mergedLiftInfo;
    return Immutable(minimalPolynomialsData);
end
);


# todo : ADJUST_PAIRING_TOLERANCE@FR also parametrizable
# tolerance is after adjusting smaller or equal to the 1/3 of the minimal distance between two roots in rootList
InstallGlobalFunction( ADJUST_PAIRING_TOLERANCE@FR ,
function (tolerance, rootList)

    local numRoots, col, row, localTolerance;

    localTolerance := tolerance;

    numRoots := Size(rootList);
    
  for row in [1..numRoots] do
    for col in [(row+1)..numRoots] do
         if AbsoluteValue(   (rootList[row] - rootList[col]) )/3.0 < localTolerance then 
            localTolerance := AbsoluteValue(rootList[row] - rootList[col])/3.0;
        fi;
    od;
    od;
    return localTolerance;
end
);


# each row should contain at least one entry (exact=false)  or exact one entry (exact=true)
InstallGlobalFunction( COMPATIBILITY_ROWS_VALID@FR ,
function(compatibiltyMatrix, exact)
    local rowSums, entry, l;
    
      rowSums := List([1..Size(compatibiltyMatrix)], i-> Number(  compatibiltyMatrix[i] , function(l) return l>0; end ) );
     for entry in rowSums do
        if entry>1 and exact then
            return false;
        fi;
        if entry<1 then
          return false;
        fi;
    od;
    return true;
end
);


# each row and each column should contain at least one entry (exact=false) or exact one entry (exact=true)
InstallGlobalFunction(  IS_VALID_ROOT_COMPATIBILITY@FR,
function( matrix, combinedRootsCount, logger )

 local mathchedRoots;
 Assert(0, Characteristic(matrix)=0);
 
 mathchedRoots := Set( FlattenList@FR(matrix));
 SubtractSet( mathchedRoots, [0] );
 
 # for each combined root there should be a existing compatibility:
  if Size(mathchedRoots)<>combinedRootsCount then
        logger( 0, "--------------root compatibility warning: Size(mathchedRoots)<>combinedRootsCount, problem with error tolerance?" );
        return  false ;
    fi;

     if not COMPATIBILITY_ROWS_VALID@FR( matrix, false) or   
        not COMPATIBILITY_ROWS_VALID@FR( TransposedMat(matrix),false)  then
             logger(0,"-------------root compatibility warning: compatibility not given;  problem with error tolerance ?");
        return false ;
    fi;    
    return true;
end
);


InstallGlobalFunction(  COMPUTE_HURWITZ_ROOT_COMPATIBILITY@FR ,
function( firstPolRoots, secondPolRoots, combinedPolRoots, operation, maxTolerance, logger)
    local localTolerance, numRoots, compatibiltyMatrix, 
          extendedCompatibilityMatrix, row, col, i, rowSums, entry;
    localTolerance := maxTolerance;
    
    Assert(0, Size(firstPolRoots)>=Size(secondPolRoots) );

    if not Size(firstPolRoots)>=Size(secondPolRoots) or 
       not Size(combinedPolRoots) = Maximum( Size(firstPolRoots), Size(secondPolRoots) ) then
          return fail;
    fi;
  
    numRoots := Size(firstPolRoots);
    compatibiltyMatrix := List( [1..numRoots] ,n-> List([1..Size(secondPolRoots)], l->0)
                                );
    extendedCompatibilityMatrix := List( [1..numRoots] ,n-> List([1..Size(secondPolRoots)], l->0)
                                );
    
    # tolerance is after adjusting smaller or equal to the minimal distance between two roots for each  root list
    localTolerance := ADJUST_PAIRING_TOLERANCE@FR( localTolerance, firstPolRoots );
    localTolerance := ADJUST_PAIRING_TOLERANCE@FR( localTolerance, secondPolRoots );
    localTolerance := ADJUST_PAIRING_TOLERANCE@FR( localTolerance, combinedPolRoots );


    if IsZero( localTolerance) then
         logger( 0,   "COMPUTE_HURWITZ_ROOT_COMPATIBILITY@FR: pairing tolerance is zero ");
         return fail;
    fi;
   
    
    for row in [1..numRoots] do
    for col in [1..Size(secondPolRoots)] do
    for i in [1..numRoots] do
        if  AbsoluteValue( operation (firstPolRoots[row], secondPolRoots[col] )- combinedPolRoots[i] ) <localTolerance then
            compatibiltyMatrix[row][col] := 1;
            extendedCompatibilityMatrix[row][col] := i;
        fi;
    od;
    od;
    od;
    
    if  not Rank( compatibiltyMatrix) = Size(secondPolRoots) then 
        return fail;
    fi;

    rowSums := List( [1..Size(compatibiltyMatrix)], i-> Sum( compatibiltyMatrix[i]) );
    for entry in rowSums do
        if not entry=1 then
          return fail;
        fi;
    od;

    return  compatibiltyMatrix;
end
);


InstallGlobalFunction(  ComputeRootCompatibilityEx@FR ,
function( firstPolRoots, secondPolRoots, combinedPolRoots, operation, maxTolerance, logger)
    local localTolerance,  compatibiltyMatrix, combinedRootsMatched,  
          row, col, i,   simpleCompatibiltyMatrix;
          
    localTolerance := maxTolerance;
    
    if  not Size(combinedPolRoots) >= Maximum( Size(firstPolRoots), Size(secondPolRoots) ) then
            logger( 1, "ComputeRootCompatibility@FR: Error: Size(combinedPolRoots)<Maximum( Size(firstPolRoots), Size(secondPolRoots) )" ); 
          return fail;
    fi;

    compatibiltyMatrix := List( [1..Size(firstPolRoots)] ,n-> List([1..Size(secondPolRoots)], l->0)
                                );

    simpleCompatibiltyMatrix := List( [1..Size(firstPolRoots)] ,n-> List([1..Size(secondPolRoots)], l->0)
                                );
    
    localTolerance := ADJUST_PAIRING_TOLERANCE@FR( localTolerance, firstPolRoots );
    localTolerance := ADJUST_PAIRING_TOLERANCE@FR( localTolerance, secondPolRoots );
    localTolerance := ADJUST_PAIRING_TOLERANCE@FR( localTolerance, combinedPolRoots );


    if IsZero( localTolerance) then
        logger( 0,   "ComputeRootCompatibility@FR: error tolerance is zero ");
        return fail;
    fi;
   
    combinedRootsMatched := List( [1..Size(combinedPolRoots)] ,n->0);
    for row in [ 1..Size(firstPolRoots)   ] do
    for col in [ 1..Size(secondPolRoots)  ] do
    for i   in [ 1..Size(combinedPolRoots)] do
        if  AbsoluteValue( operation (firstPolRoots[row], secondPolRoots[col] )- combinedPolRoots[i] ) <localTolerance then
            combinedRootsMatched[i] := 1;
            compatibiltyMatrix[row][col] := i;
            simpleCompatibiltyMatrix[row][col] := 1;
        fi;
    od;
    od;
    od;   

    logger( 2, "compatibiltyMatrix");
    logger( 2, String(compatibiltyMatrix) );

    if not IS_VALID_ROOT_COMPATIBILITY@FR( compatibiltyMatrix, Size(combinedPolRoots), logger  ) then 
        logger( 0, "--------------ComputeRootCompatibility@FR:   probably a problem with error tolerance....." );
        return fail;
    fi;

    return  compatibiltyMatrix;
end
);


InstallGlobalFunction(  ComputeRootCompatibility@FR ,
function( firstPolRoots, secondPolRoots, combinedPolRoots, operation, maxTolerance)
	return ComputeRootCompatibilityEx@FR(  firstPolRoots, secondPolRoots, combinedPolRoots, 
	                                   operation, maxTolerance, CREATE_EMPTY_LOGGER_FKT@FR() 
	                                );
end
);

 
InstallGlobalFunction( IDEAL_POINTS_APPROXIMATION@FR,
function( minPolyData, approxSolutions, mergedLiftInfo  )
    local approxSolutionData,indeterminates, errorList,root, error, unknownMinPolyData;
    
    approxSolutionData:=rec();
    approxSolutionData.approxIdealElems  := Immutable(approxSolutions);
    approxSolutionData.minPolyData := Immutable(minPolyData);
    approxSolutionData.mergedLiftInfo := Immutable(mergedLiftInfo);
   
    
    errorList := [];
    indeterminates := List(minPolyData.unknowns, i->i[1] );
    for unknownMinPolyData in minPolyData.unknowns do
        for root in  approxSolutions do
            error := EvalPolynomialTensor@FR( unknownMinPolyData[2],  indeterminates, root);
            Append(errorList,[ AbsoluteValue( error) ] );
        od;
    od;
    approxSolutionData.residue := Maximum( errorList );
    
    approxSolutionData.dataType := Immutable("IdealPointsApprox"); 
    return Immutable(approxSolutionData);
end
);

# todo: it might be that the minimal polynomials are already computed and one does not want to compute them again. redesign.
InstallGlobalFunction( ComputeApproxIdealPoints@FR ,
function( inputIdeal,  solutionPoint , opts)
   
    local   minimalPolynomialsData, mergedLiftInfo, rootListList, 
            operation, operationInputList,  operationUsedList, 
            unknown, newUnknown,   unknownIdx, referenceRoots,  unknownRoots,  
            preApproxSolutions, tmppreApproxSolutions, 
            compatibilityResult,  compMatrix,   row, col, entry, entryCopy, 
            currentCoordinatePaired, idealRing,  indeterminates, approxSolutionData;
  
    idealRing  :=  LeftActingRingOfIdeal (inputIdeal);
    Assert(0, idealRing=RightActingRingOfIdeal (inputIdeal) );
    indeterminates := IndeterminatesOfPolynomialRing(idealRing);

    
    minimalPolynomialsData := ComputeMinimalPolynomials@FR( inputIdeal, solutionPoint, indeterminates, opts);
    if fail=minimalPolynomialsData then
    	return fail;
    fi;
  
    mergedLiftInfo :=  minimalPolynomialsData.mergedLiftInfo;

    opts.setInitialLiftDepth(  mergedLiftInfo.maxLiftDepth+1 ); # is a heuristic. could be suboptimal for generic problems. 
    opts.setInitialLatticeDim ( mergedLiftInfo.requiredLatticeDimension) ;

    
    opts.logger(1, "------------------------pairing part---------------------------") ;
    if not mergedLiftInfo. requiredLatticeDimension=0  then 
    
        # compute roots for each minimalPolynomial ( unknowns[i][2] ) 
        rootListList := List([ 1..Size(indeterminates)] , unknownIdx->opts.rootCalculator().computeRoots( minimalPolynomialsData.unknowns[unknownIdx][2]) );  
      
        operationInputList := List( [1..Characteristic(solutionPoint)-1], fieldNonzeroElem-> function(a,b) return a + fieldNonzeroElem*b; end ) ; 
        operationUsedList  := []; #debugging

        unknown := indeterminates[1];

        referenceRoots := rootListList[1];

       preApproxSolutions := List( [1..Size(referenceRoots)], n-> [[ referenceRoots[n] ]] );

       for unknownIdx in [2..Size(indeterminates)] do
        
            if opts.verbosePairing() then
                opts.logger(1, Concatenation("unknownIdx", String(unknownIdx)) );
            fi;
            currentCoordinatePaired := false;
            for operation in operationInputList do
                    newUnknown :=  operation ( unknown , indeterminates[unknownIdx] );
                    opts.logger(2, Concatenation("newUnknown: ", String(newUnknown  ) ) ) ;

                    opts.setInitialLatticeDim ( 1+ Size(preApproxSolutions) );
                    
                    # adjust 'maxLatticeDim': the worst situtation would be if each root in  preApproxSolutions is compatible with each root in 'rootListList[unknownIdx]'
                    opts.setMaxLatticeDim ( 1+ Size(preApproxSolutions)*Size( rootListList[unknownIdx] ) );
                  
                    opts.logger(1, Concatenation("opts.maxLatticeDim: ", String(opts.maxLatticeDim ) ) ) ;
                    compatibilityResult := ComputeMinimalPolynomials@FR( inputIdeal, solutionPoint, [newUnknown], opts);
                    
                    if fail=compatibilityResult then 
                      continue;  
                    fi;
                    
                    opts.logger(1, Concatenation(" ----------------pairing variable ",String(unknownIdx) ) );
                    unknownRoots := opts.rootCalculator().computeRoots( compatibilityResult.unknowns[1][2]);
                
             
                    compMatrix := ComputeRootCompatibilityEx@FR( referenceRoots, 
                                                                   rootListList[unknownIdx] ,  
                                                                   unknownRoots  ,
                                                                   operation, 
                                                                   opts.maxPairingTolerance(), 
                                                                   opts.logger );
                    if not fail=compMatrix then
                            opts.logger(2,  "---------------------------compatibility matrix---------------------------------");
                            opts.logger(2, String(compMatrix) );
                        tmppreApproxSolutions := List([ 1..Size(unknownRoots)], n->[] );
                        
                        for row in [1..Size(compMatrix) ] do
                        for col in  [1..Size(compMatrix[1]) ] do
                            if compMatrix[row][col]>0 then
                                for entry in preApproxSolutions[row] do
                                    entryCopy := ShallowCopy(entry);
                                    Append( entryCopy,[ rootListList[unknownIdx][col] ] );
                                    Append( tmppreApproxSolutions[ compMatrix[row][col] ], [ entryCopy ]  );
                                od;
                            fi;
                        od;
                        od;
                        preApproxSolutions := tmppreApproxSolutions;
                        Append( operationUsedList, [operation] );       
                        referenceRoots := unknownRoots;
                        unknown := newUnknown;
                        currentCoordinatePaired := true;
                        mergedLiftInfo := MERGE_LIFT_INFO@FR(  compatibilityResult.mergedLiftInfo , mergedLiftInfo);
                        opts.logger(1, " ----------------pairing success---------------------------\n");
                        opts.logger(1, Concatenation("unknownIdx: ", String(unknownIdx) )  );
                        break;
                    fi;
            od;
            if not currentCoordinatePaired then
                opts.logger(0, Concatenation("pairing failed for indeterminate ", String(unknownIdx) ));
                return fail;
            fi;
        od;
    fi;
    opts.logger (1, " ---------------- All variables paired !---------------------------\n");
    # todo: save input parameters in the result or not?
     # debugInfo := Immutable (rec ( operationsUsedForPairing := operationUsedList ) );
     approxSolutionData := IDEAL_POINTS_APPROXIMATION@FR( minimalPolynomialsData, FlattenList@FR( preApproxSolutions ) , mergedLiftInfo );
    
    return approxSolutionData;
  
end
);



#  COMPUTE_APPROX_HURWITZ_IDEAL_POINTS@FR: 
# -see also 'ComputeApproxIdealPoints@FR' .
# may run faster than the generic version, but not succeed for all cases!
# precondition: assumes that number of solutions of the first indeterminate (Degree of its minimal polymomial) 
# is the same as the number of   all paired coordinates.
# thus wont work for each situation, but may work for HurwitzMapSearch problems !
#
InstallGlobalFunction( COMPUTE_APPROX_HURWITZ_IDEAL_POINTS@FR,
function( inputIdeal, solutionPoint , opts)

    local     minimalPolynomialsData, mergedLiftInfo,  pairedRootRootList, rootListList,
     operation,  operationInputList, operationUsedList,   unknown, unknownIdx, 
     compatibilityResult,    compMatrix, modDstRootList, roots, approxIdealElems,
     paired, idealRing, indeterminates, approxSolutionData;


    Assert(0, LeftActingRingOfIdeal (inputIdeal)=RightActingRingOfIdeal (inputIdeal) );   
    idealRing  :=  LeftActingRingOfIdeal (inputIdeal);
    indeterminates := IndeterminatesOfPolynomialRing(idealRing);

    
    minimalPolynomialsData := ComputeMinimalPolynomials@FR( inputIdeal, solutionPoint, indeterminates, opts);
    if minimalPolynomialsData=fail then
    	Info(InfoFR,1, "failed to compute minimal polynomials");
    	return fail;
    fi;
    mergedLiftInfo :=  minimalPolynomialsData.mergedLiftInfo;
    
    #opts.setInitialLiftDepth( mergedLiftInfo.maxLiftDepth );
    opts.setInitialLiftDepth( mergedLiftInfo.minLiftDepth );
    opts.setInitialLatticeDim( mergedLiftInfo.requiredLatticeDimension ) ;
    opts.setMaxLatticeDim ( mergedLiftInfo.requiredLatticeDimension );
    # opts.setMaxLatticeDim ( mergedLiftInfo.requiredLatticeDimension^2 ); leads to memory error!

    pairedRootRootList := List( [1..Size(indeterminates)], n->0) ;

    if not mergedLiftInfo. requiredLatticeDimension=0  then 
    
        # compute roots for each minimalPolynomial ( unknowns[i][2] ) 
        rootListList := List([ 1..Size(indeterminates)] , unknownIdx->opts.rootCalculator().computeRoots( minimalPolynomialsData.unknowns[unknownIdx][2]) );  
      
        pairedRootRootList[1] := rootListList[1];

        # todo: is a+c*b sufficient or is it also required c*a+d*b?
        operationInputList := List( [1..Characteristic(solutionPoint)-1], pos-> function(a,b) return a+pos*b; end ) ; 

        # for debugging:
        operationUsedList := [];

       for unknownIdx in [2..Size(indeterminates)] do
            opts.logger (1, Concatenation(" ---------------- (special) Pairing variable ", String(unknownIdx))) ;
            paired := false;
            for operation in operationInputList do
                    unknown :=  operation ( indeterminates[1] , indeterminates[unknownIdx] );
                    compatibilityResult := ComputeMinimalPolynomials@FR( inputIdeal,  solutionPoint, [unknown], opts);
                    if fail=compatibilityResult then 
                    	continue;
                    fi;
             
               
                    mergedLiftInfo := MERGE_LIFT_INFO@FR( compatibilityResult.mergedLiftInfo, mergedLiftInfo );
                    roots := opts.rootCalculator().computeRoots( compatibilityResult.unknowns[1][2]);
                    if not Size(roots) = Size(rootListList[1]) then 
                        continue;
                    fi;
                 
                    compMatrix := COMPUTE_HURWITZ_ROOT_COMPATIBILITY@FR( rootListList[1], 
                    							     rootListList[ unknownIdx ],  
                    							     roots  ,
                    							     operation, 
                    							     opts.maxPairingTolerance(),
                    							     opts.logger  );
                      if not fail=compMatrix then
                        if opts.verbosePairing() then
                             opts.logger (1, "compMatrix");
                             opts.logger (1, compMatrix);
                        fi;
                      
                        modDstRootList := compMatrix*TransposedMat( [ rootListList[unknownIdx] ] );
                        Append( operationUsedList, [operation] );       
                        Assert(0, Size(TransposedMat( modDstRootList ))=1);
                        pairedRootRootList[unknownIdx] := TransposedMat(modDstRootList)[1];    
                        opts.logger (1, " ----------------Pairing success---------------------------\n");
                        opts.logger(1, Concatenation("paired unknownIdx: ", String(unknownIdx) ) ) ;
                        paired := true;
                        break;
                    fi;
            od;
            if not paired then
                Error (Concatenation("pairing failed for unknownIdx ", String(unknownIdx) ));
            fi;
        od;
    fi;
     opts.logger (1, " ---------------- All variables paired !---------------------------\n");
    # check:
    for roots in pairedRootRootList do
        Assert(0, Size(roots) = mergedLiftInfo.requiredLatticeDimension-1);
    od;
    
    # compose approximate ideal elements from all coordinates. 
    approxIdealElems := List( [1..Size( pairedRootRootList[1] )], rootIdx-> List( [1..Size(indeterminates)], unknownIdx-> pairedRootRootList[unknownIdx][rootIdx] )
                );
                
    #debugInfo := (rec ( operationsUsedForPairing := operationUsedList ) );
    approxSolutionData := IDEAL_POINTS_APPROXIMATION@FR( minimalPolynomialsData, Immutable(approxIdealElems), mergedLiftInfo );
    
    return approxSolutionData;
end
);

##########################################################################################################################################################################


InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_LLL", 
function()
    local mat,lllResult;
    mat:=[[1,2],[2,1]];
    lllResult:= FPLLLReducedBasis(mat);
    Assert(0, lllResult=[ [ 1, -1 ], [ 1, 2 ] ] );
end
);


InstallGlobalFunction( CREATE_FINITE_TEST_PROBLEM@FR ,
function()
    local  rng, indeterminates,x,y, FZ1,FZ2, ideal, solutionOverFiniteField, expectedResult, problem ;

    rng := PolynomialRing( ZmodnZ(11)  ,["x","y"] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
    y := indeterminates[2];
    FZ1 := 33*x^3+19*x^2-81*x-4;
    FZ2 := y-1;
    ideal := Ideal(rng,[FZ1,FZ2]);
    solutionOverFiniteField := [ Z(11)^0, Z(11)^0 ];

    problem := rec();
    problem.ideal := ideal;
    problem.indeterminates := indeterminates;   
    problem.solution := solutionOverFiniteField;      
    problem.unknowns := indeterminates;
    return problem;
end
);


InstallGlobalFunction( CREATE_RATIONAL_TEST_PROBLEM@FR ,
 function()
    local  rng, indeterminates,x,y, FZ1,FZ2, ideal, solutionOverFiniteField, expectedResult, problem ;

    rng := PolynomialRing( Rationals  ,["x","y"] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
    y := indeterminates[2];
    FZ1 := 33*x^3+19*x^2-81*x-4;
    FZ2 := y-1;
    ideal := Ideal(rng,[FZ1,FZ2]);
    solutionOverFiniteField := [ Z(11)^0, Z(11)^0 ];

    problem := rec();
    problem.ideal := ideal;
    problem.indeterminates := indeterminates;   
    problem.solution := solutionOverFiniteField;      
    problem.unknowns := indeterminates;
    return problem;
end
);


InstallGlobalFunction( CREATE_SYMM_TEST_PROBLEM@FR ,
function()
    local  rng, indeterminates,x,y, FZ1,FZ2, ideal, solutionOverFiniteField, expectedResult, problem ;

    rng := PolynomialRing( Rationals  ,["y","x"] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[2];
    y := indeterminates[1];
    FZ1 := 33*x^3+19*x^2-81*x-4;
    FZ2 := y-1;
    ideal := Ideal(rng,[FZ1,FZ2]);
    solutionOverFiniteField := [ Z(11)^0, Z(11)^0 ];

    problem := rec();
    problem.ideal := ideal;
    problem.indeterminates := indeterminates;   
    problem.solution := solutionOverFiniteField;      
    problem.unknowns := indeterminates;
    return problem;
end
);



InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_JENKINS_TRAUB_USAGE", 
function()
    local rng, indeterminates, x, y, FZ1, bitPrecision, roots, rootCalculator;
  
    rng := PolynomialRing( Rationals  ,["x"] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];

      
    FZ1 := 33*x^3+19*x^2-81*x-4;

    roots:= RootsByJenkinsTraub@FR(FZ1,16);
    
    roots:= RootsByJenkinsTraub@FR(FZ1,320);
    roots:= RootsByJenkinsTraub@FR(FZ1,330);


    rootCalculator := CreateJenkinsTraubWrapper@FR(16);

    roots := rootCalculator.computeRoots(FZ1);
    
    roots := rootCalculator.computeRoots(FZ1);
    roots := rootCalculator.computeRoots(FZ1);

end
);



InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_LIFT_STEP_1@FR", 
function()
    local  rng,jac,ind,x,FZ,ideal,finiteField,solution,gens;
    rng := PolynomialRing( Rationals  ,["x"] );
    ind := IndeterminatesOfPolynomialRing(rng);
    x := ind[1];
    FZ := 33*x^3+19*x^2-81*x-4;
    ideal := Ideal(rng,[FZ]);
    jac := Jacobian@FR( [FZ] ,ind );
    solution := [Z(11)^0];
    gens := GeneratorsOfTwoSidedIdeal( ideal );
    Assert(0, IsZero( Value(FZ,ind,solution)) );
    Assert(0, IsZero( EvalPolynomialTensor@FR(gens,ind,solution)) );
    solution := QuadraticLiftStep@FR( gens, jac, ind,  solution);
    Assert(0, IsZero( EvalPolynomialTensor@FR(gens,ind,solution)) );
end
);


InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_BLACKBOX_LIFT_STEP_1", 
function()
    local  rng,jac,ind,x,FZ,ideal,finiteField,solution,gens, evalIdealGens, jacobianAt;
    rng := PolynomialRing( Rationals  ,["x"] );
    ind := IndeterminatesOfPolynomialRing(rng);
    x := ind[1];
    FZ := 33*x^3+19*x^2-81*x-4;
    ideal := Ideal(rng,[FZ]);

    jac := Jacobian@FR( [FZ] ,ind );
    solution := [Z(11)^0];
    gens := GeneratorsOfTwoSidedIdeal( ideal );

    Assert(0, IsZero(Value( FZ,ind,solution)) );
    Assert(0, IsZero(EvalPolynomialTensor@FR( gens,ind,solution)) );

    evalIdealGens := function(point)
        return EvalPolynomialTensor@FR( gens, ind, point);
    end;

    jacobianAt := function( point )
        return EvalPolynomialTensor@FR( jac, ind, point);
    end;

    solution := BlackBoxQuadraticLiftStep@FR( evalIdealGens, jacobianAt, solution);
    Assert(0, IsZero(EvalPolynomialTensor@FR(gens,ind,solution)) );
end
);


InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_LIFT_STEP_2", 
function()
    local  problem, jac,solution, gens ;
   
    problem :=  CREATE_RATIONAL_TEST_PROBLEM@FR();

    gens := GeneratorsOfTwoSidedIdeal( problem.ideal );
    jac := Jacobian@FR( gens , problem.indeterminates );
    Assert(0, IsZero( EvalPolynomialTensor@FR(gens, problem.indeterminates, problem.solution)) );
    solution := QuadraticLiftStep@FR( gens, jac,  problem.indeterminates, problem.solution);
    Assert(0, IsZero( EvalPolynomialTensor@FR(gens, problem.indeterminates, solution)) );
end
);


InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_PADIC_LIFT", 
function()
    local  problem, solution, gens;

    problem :=  CREATE_RATIONAL_TEST_PROBLEM@FR();
    solution := PadicLift@FR( problem.ideal,  problem.solution, 3);
    gens := GeneratorsOfTwoSidedIdeal( problem.ideal );
    Assert(0, IsZero(EvalPolynomialTensor@FR(gens, problem.indeterminates, solution)) );
    # return solution;
end
);


 
InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_BLACKBOX_PADIC_LIFT", 
function()
    local  problem, solution, gens, jac, jacobianAt, evalIdealGens;

    problem :=  CREATE_RATIONAL_TEST_PROBLEM@FR();
    gens := GeneratorsOfTwoSidedIdeal( problem.ideal );

    evalIdealGens := function(point)
        return EvalPolynomialTensor@FR( gens, problem.indeterminates,  point );
    end;

    jac := Jacobian@FR( gens , problem.indeterminates );

    jacobianAt := function( point )
        return EvalPolynomialTensor@FR( jac, problem.indeterminates, point );
    end;

    solution := BlackBoxPadicLift@FR( evalIdealGens, jacobianAt, problem.solution, 3 );
   
    Assert(0, IsZero(evalIdealGens(solution)) );
    # return solution;
end
);



InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_LLL_REDUCTION", 
function()
  local  problem, solution, liftResult, nextLiftResult, gens, reductionOpts;

    problem :=  CREATE_RATIONAL_TEST_PROBLEM@FR();
     
    liftResult := PadicLift@FR( problem.ideal,  problem.solution, 3);
    nextLiftResult := PadicLift@FR( problem.ideal,  problem.solution, 4);

    gens := GeneratorsOfTwoSidedIdeal( problem.ideal );
    Assert(0, IsZero( EvalPolynomialTensor@FR(gens, problem.indeterminates, liftResult)) );
    Assert(0, IsZero( EvalPolynomialTensor@FR(gens, problem.indeterminates, nextLiftResult)) );
    
    reductionOpts := LiftOptions@FR();
    LLL_REDUCTION_ATTEMPT@FR ( problem.unknowns[1], problem.indeterminates, liftResult, nextLiftResult, reductionOpts );
end
);


InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_COMPUTE_MINIMAL_POLYNOMIAL", 
function()
  local  problem, options, unknown, minimalPolynomialVariable, liftAndLLLRes;
   
    problem :=  CREATE_RATIONAL_TEST_PROBLEM@FR();
      
    options := LiftOptions@FR();
    unknown := problem.indeterminates[1];
    minimalPolynomialVariable := Indeterminate(Rationals);

    liftAndLLLRes := ComputeMinimalPolynomialEx@FR ( problem.ideal,  problem.solution, unknown,  minimalPolynomialVariable, options );

    unknown := problem.indeterminates[2];
    liftAndLLLRes := ComputeMinimalPolynomialEx@FR (problem.ideal,  problem.solution, unknown, minimalPolynomialVariable,  options );
end
);


InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_COMPUTE_MINIMAL_POLYNOMIALS", 
function()
  local  indeterminates,x,y,FZ1,FZ2, ideal, solutionOverFiniteField, liftAndLLLOptions,  
   expectedUnknowns, expectedMergedLiftInfo, problem, unknowns, liftAndLLLRes ;

    liftAndLLLOptions := LiftOptions@FR();

    problem :=  CREATE_RATIONAL_TEST_PROBLEM@FR();
    x :=  problem.indeterminates[1];
    y :=  problem.indeterminates[2];

    liftAndLLLRes := ComputeMinimalPolynomials@FR ( problem.ideal, problem.solution,    problem.unknowns ,  liftAndLLLOptions );
    expectedUnknowns :=   [ [ x, -11*x^2-21*x-1 ], [ y, y-1 ] ];
    expectedMergedLiftInfo := rec( dataType:="LiftInfo", maxLatticeDimension := 3, maxLiftDepth := 3, requiredLatticeDimension := 3 );
    
    Assert( 0, liftAndLLLRes.unknowns=expectedUnknowns );
    Assert( 0, liftAndLLLRes.mergedLiftInfo=expectedMergedLiftInfo );
    
end
);


InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_COMPATIBILITY_ROWS_VALID", 
function()
    local matrix;
	matrix := [[0,1],[2,0],[0,3]];
	Assert(0, COMPATIBILITY_ROWS_VALID@FR(matrix,false));
	Assert(0, COMPATIBILITY_ROWS_VALID@FR(matrix,true));
	matrix := [[2,1],[2,0],[0,3]];
	Assert(0, COMPATIBILITY_ROWS_VALID@FR(matrix,false));
	matrix := [[2,1],[2,0],[0,3]];
	Assert(0, not COMPATIBILITY_ROWS_VALID@FR(matrix,true));
	matrix := [[2,1],[0,0],[0,3]];
	Assert(0, not COMPATIBILITY_ROWS_VALID@FR(matrix,true));
	Assert(0, not COMPATIBILITY_ROWS_VALID@FR(matrix,false));
end
);


InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_IS_VALID_ROOT_COMPATIBILITY", 
function()

   local logger,matrix;
    logger := function(a,b) end;
   
    matrix:= [[1,2],[1,4],[5,6]];
    Assert(0, false=IS_VALID_ROOT_COMPATIBILITY@FR(matrix,6,logger) );
    matrix:= [[1,2],[3,4],[5,6]];
    Assert(0, true=IS_VALID_ROOT_COMPATIBILITY@FR(matrix,6,logger) );
    Assert(0, true=IS_VALID_ROOT_COMPATIBILITY@FR(matrix,6,logger) );
    
     matrix:= [[1,0],[3,0],[2,0]];
     Assert(0, false=IS_VALID_ROOT_COMPATIBILITY@FR(matrix,3,logger) );
     
      matrix:= [[1,0],[0,3],[2,0]];
     Assert(0, true=IS_VALID_ROOT_COMPATIBILITY@FR(matrix,3,logger) );
  
end
);



InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_COMPUTE_ROOT_COMPATIBILITY", 
function()
    local firstPolRoots, secondPolRoots, combinedPolRoots, operation, opts, compatibility;
    # Probleme, wenn SetFloats nicht aufgerufen worden ist...   
    SetFloats(MPC,1000);  

    firstPolRoots:=[ 0.03, 34.0, 10.0 ];
    secondPolRoots:=[ 5.03, 4.0, 1.0 ];   
    combinedPolRoots := [ 4.03, 11.0, 39.02 ];
    
    operation := function(a,b) return a+b; end;
    opts := LiftOptions@FR();
    opts.setMaxPairingTolerance ( 0.001 );
    compatibility := COMPUTE_HURWITZ_ROOT_COMPATIBILITY@FR ( firstPolRoots, secondPolRoots, combinedPolRoots, operation, opts.maxPairingTolerance(), opts.logger );
    Assert(0, compatibility = fail);

    opts:=LiftOptions@FR();
    opts.setMaxPairingTolerance ( 0.02 );      
    opts.setVerbosePairing (false );
    compatibility := COMPUTE_HURWITZ_ROOT_COMPATIBILITY@FR ( firstPolRoots, secondPolRoots, combinedPolRoots, operation, opts.maxPairingTolerance(), opts.logger );
    Assert(0, compatibility= [[ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ] );

    firstPolRoots  := [  4.0, 10.0 ];
    secondPolRoots := [ 5.0 ];   
    combinedPolRoots := [ 9.0, 15.0 ];

    compatibility := ComputeRootCompatibilityEx@FR ( firstPolRoots, secondPolRoots, combinedPolRoots, operation, opts.maxPairingTolerance(), opts.logger );
    Assert(0, compatibility = [[1],[2]]);
end
);


InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_COMPUTE_APPROX_IDEAL_POINTS", 
function()
    local TestHelper;
    TestHelper := 
    function(problem)
        local    opts, gens, result, errorTolerance, evaluation, evaluationAbs, max , root;
        
        opts := LiftOptions@FR();

        result := ComputeApproxIdealPoints@FR( problem.ideal,   problem.solution, opts);

        gens := GeneratorsOfTwoSidedIdeal( problem.ideal );

        errorTolerance := 1.e-14;

        for root in result.approxIdealElems do
            evaluation := EvalPolynomialTensor@FR( gens,problem.indeterminates, root ) ;
            #evaluationAbs := List([1..Size(evaluation)], n-> AbsoluteValue( evaluation[n]) );
            evaluationAbs := List( evaluation, n-> AbsoluteValue( n) );
            max := Maximum(evaluationAbs);
            Assert(0, max<errorTolerance );
        od;
    end;
    # problem:=  CREATE_RATIONAL_TEST_PROBLEM@FR();
    TestHelper(  CREATE_RATIONAL_TEST_PROBLEM@FR() );
    TestHelper(  CREATE_SYMM_TEST_PROBLEM@FR() );
end
);



InstallGlobalRecordFunction@FR ( ["@PadicLift","Tests"], "TEST_COMPUTE_HURWITZ_APPROX_IDEAL_POINT", 
function()
    local   problem, opts, gens, result, errorTolerance, evaluation, evaluationAbs, max, root ;

    problem :=  CREATE_RATIONAL_TEST_PROBLEM@FR();
    
    opts := LiftOptions@FR();

    result := COMPUTE_APPROX_HURWITZ_IDEAL_POINTS@FR(problem.ideal, problem.solution, opts);

    gens := GeneratorsOfTwoSidedIdeal( problem.ideal );

    errorTolerance := 1.e-14;

    for root in result.approxIdealElems do
        evaluation := EvalPolynomialTensor@FR( gens,problem.indeterminates, root ) ;
        evaluationAbs := List( evaluation, n-> AbsoluteValue( n) );
        max := Maximum(evaluationAbs);
        Assert(0, max<errorTolerance );
    od;
end
);



# todo: depends on package 'float' - write this test elsewhere !
InstallGlobalRecordFunction@FR (["@PadicLift","Tests"], "TEST_COERCE_POLYNOMIAL_TO_COMPLEX_RING",
 function()
    local rng, indeterminates, x, pol, dstrng, coercedPol,dstInd, expectedResult;
    rng := PolynomialRing(Rationals,1);
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
    pol := x^+2+3;

    dstrng := PolynomialRing( MPC_PSEUDOFIELD,1 ); # 1 indeterminate
    
    coercedPol := CoercePolynomialTensor@FR(pol, dstrng);
    
    dstInd := IndeterminatesOfPolynomialRing(dstrng);
    expectedResult := dstInd[1]^2+3.0_c;
    Assert(0, coercedPol=expectedResult);
    
end
);



 
InstallGlobalRecordFunction@FR (["@PadicLift"], "CreateTestString",
function(prefix)
    return @FR@Utils.Internal.CreateTestString("@PadicLift.Tests", prefix);
end
);



