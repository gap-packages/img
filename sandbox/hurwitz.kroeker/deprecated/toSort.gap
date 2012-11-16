#SetPackageInfo( rec(
#  PackageName := "rationalMapFinder",
#  Version := "1.0beta",
#  AvailabilityTest := ReturnTrue,
#  Autoload := false,
#  BannerString := Concatenation( [
#      "#I  loading the GAP package ", ~.PackageName," in version ",
#      ~.Version, "\n" ] ),
#  PackageDoc := rec(
#      BookName  := "rationalMapFinder",
#      SixFile   := "doc/manual.six",
#      Autoload  := true ) ) );


get21212121ProblemFkt := function()
	local shape,shapeList,scalingRelationList,rfsProblem;
	shape := createShape([2,1]);
	shapeList := [ shape,shape,shape,shape ];
	scalingRelationList:=[[0,-1]];
	rfsProblem := rec(shapeList:=shapeList ,  scalingRelationList:= scalingRelationList);
	Assert(true,IsRFSProblem(rfsProblem));
	return rfsProblem;
end;
BindGlobal("get21212121Problem",get21212121ProblemFkt);

getRFSTestProblemFkt := function()
	local shape,shapeList,scalingRelationList,rfsProblem;
	shape := createShape([2,1]);
	shapeList := [ shape,shape,shape,shape ];
	scalingRelationList:=[ ];
	rfsProblem:=rec(shapeList:=shapeList ,  scalingRelationList:= scalingRelationList);
	Assert(true,IsRFSProblem(rfsProblem));
	return rfsProblem;
end;
BindGlobal("getRFSTestProblem",getRFSTestProblemFkt);



testM2RationalMapSearch := function()

	local m2SearchPathList, rfsObj, partitionList,fieldChar, polTupleList, polTuple, liftOptions, liftedPolTuple,ind,expectedResult,preRationalMap,rationalMap,strictNormalization;
	
	# Read("rfsWrapper.gap");
	
	LoadPackage("float");
	SetFloats(MPC);;
	############# init ###################################
	
	# create an object with the interface 'computeFiniteFieldSolutions(..)' and 'approximateComplexSolutions(..)':
	#rfsObj := createM2RationalMapFinder( "M2", [ Directory("/home/kroeker/bin/") ] );  	
	rfsObj := createM2RationalMapFinder( "M2", [   ] );  
	rfsObj.verbose:=true;
	
	############# 1. (smart) brute force search for polynomialSets ([A,B,C],...) over given finite field Fp where 
	#             polynomials A,B and C matches given shapes, gcd(A,B)=gcd(A,C)=gcd(B,C)=1 and  B-lambda*A-C = 0 for some lambda in Fp
	
	fieldChar := 7;	partitionList := [ [3], [2,1], [1,2] ]; 
	strictNormalization := true;
	polTupleList := rfsObj.computeFiniteFieldSolutions( partitionList, strictNormalization, [], fieldChar );
	
	############# 2. try to compute for a result from step (1) a lift to a polynomial ring over extension of Q and a complex approximation 
	
	polTuple :=  polTupleList[1];
		
	# design decision: liftOptions is an extra object => interface keeps the same while  'liftOptions' is extensible!
	# todo: createLiftOptions.
	liftOptions := rec ( decimalPrecision := 16 ); # all following parameters are optional, 
							# but if   no 'maxLiftDepth' given or was chosen too big,  the computation may run forever and/or consume all memory.
	liftOptions.maxLiftDepth  := 10;  # lift up to    mod fieldChar^(2^maxLiftDepth) ; 
	liftOptions.maxLatticeDim := 100; # if no  'maxLatticeDim'  given or was chosen too big, the computation may run forever and/or consume all memory.
		
	liftedPolTuple := rfsObj.approximateComplexSolutions( polTuple, liftOptions );
	
	## the test part:
	
	preRationalMap := preRationalMapFromRootDataElem( liftedPolTuple.rootData[1], liftedPolTuple.liftedPolynomialRing);
	ind := IndeterminatesOfPolynomialRing( liftedPolTuple.liftedPolynomialRing);
	expectedResult := 3.0_c*(ind[1]^2)+(-2.0_c*(ind[1]^3));
	rationalMap := (preRationalMap[1].numerator)/( preRationalMap[1].denominator);
	Assert(true, rationalMap = expectedResult );
	
	
	# todo: output in case the lift fails.
	
	# liftedPolTuple.rootData containts a list  with the preimages of (infty, zero and 1) respectively
        # the last element of 'rootData' contains the scaling factors [lambda,mue,...] : A-lambda*B=C;  A-mue*B=D; ...
	
	############################################
	
end;

rationalMapSearch21212121Example := function()

	local m2SearchPathList, rfsObj, partitionList,fieldChar, polTupleList, polTuple, liftOptions, liftedPolTuple,branchValueApproxList,strictNormalization,preRationalMap;
	
	# Read("rfsWrapper.gap");
	
	LoadPackage("float");
	SetFloats(MPC);;
	############# init ###################################
	
	# create an object with the interface 'computeFiniteFieldSolutions(..)' and 'approximateComplexSolutions(..)':
	#rfsObj := createM2RationalMapFinder( "M2", [ Directory("/home/kroeker/bin/") ] );
	rfsObj := createM2RationalMapFinder(  "M2",[   ] );  	  	
	
	############# 1. (smart) brute force search for polynomialSets ([A,B,C,D],...) over given finite field Fp where 
	#             polynomials A,B C and D matches given shapes, and for all pairs from [A,B,C,D]  gcd  is one ,
	#	      B - lambda*A = C for some lambda in Fp and
	#	      B - mue*A    = D with mue = -i*lambda  ( -i is determinated via 'branchValueApproxList'-parameter:
	#                                                   branchValueApproxList[1][1] := RealPart(-i), branchValueApproxList[1][2] := ImaginaryPart[-i] ).
 	# 						    First three branch values are omitted and assumed as normalized to [infinity, 0, 1 ]. 
	
	fieldChar := 13;	partitionList := [ [1,2], [2,1], [2,1], [2,1] ]; 
	branchValueApproxList := [ [0/1, -1/2] ]	; # first three branch values ommitted and are assumed [infinity, 0, 1 ].
	strictNormalization :=false; # if false, the algorithm decides, which factors will be normalized to [infinity, 0, 1 ],
				     # otherwise first entries of the first three partitions in 'partitionList' determine which factors to normalize.
	polTupleList := rfsObj.computeFiniteFieldSolutions( partitionList, strictNormalization, branchValueApproxList, fieldChar );
	
	############# 2. try to compute for a result from step (1) a lift to a polynomial ring over extension of Q and a complex approximation 
	
	polTuple :=  polTupleList[1];
		
	# design decision: liftOptions is an extra object => interface keeps the same while  'liftOptions' is extensible!
	# todo: createLiftOptions.
	liftOptions := rec ( decimalPrecision := 16 ); # all following parameters are optional, 
							# but if   no 'maxLiftDepth' given,  the computation may run forever and/or consume all memory.
	liftOptions.maxLiftDepth  := 10;  # lift up to    mod fieldChar^(2^maxLiftDepth) ; 
	liftOptions.maxLatticeDim := 100; # if no  'maxLatticeDim'  given, the computation may run forever and/or consume all memory.
		
	liftedPolTuple := rfsObj.approximateComplexSolutions( polTuple, liftOptions );
	
	preRationalMap := preRationalMapFromRootDataElem( liftedPolTuple.rootData[1], liftedPolTuple.liftedPolynomialRing);
	
	# todo: output in case the lift fails.
	
	# liftedPolTuple.rootData containts the preimages of (infty, zero and 1) respectively
        # the last array contains the scaling factors [lambda,mue,...] : A-lambda*B=C;  A-mue*B=D; ...
	
	############################################
	
end;


checkProblem := function()
	local rfsProblem, searchOptions,resultFileName,rfs,resultVariableName,
	outputOptions,macaulay2Path,M2SearchPathList;
	Read("rfsWrapper.gap");
	rfsProblem := get21212121Problem();
	rfsProblem := 	getRFSTestProblem();
	searchOptions := createDefaultRFSOptionList();
	checkRFSOptionListConsistency(searchOptions);
	
	outputOptions:=rec();
	outputOptions.resultFileName := "43222RFSsearchResult.gap";
	
	outputOptions.resultVariableName := "rfsObj";
	
	#macaulay2Path :=  Null ;
	macaulay2Path :=  "/home/kroeker/bin/" ;
	
	M2SearchPathList:=[ Directory(macaulay2Path) ];
	
	performRFSSearchOverFiniteFields( rfsProblem,searchOptions,outputOptions,M2SearchPathList);
	Read("43222RFSsearchResult.gap"); 	#sets variable $resultVariableName  (here "rfsObj") !  
	


end;


	
# den Namen der Funktion bzw der Ergebisvariabe die Macaulay in die Datei schreiben solk könnte GAP an Macaulay als Parameter übergben
checkLift :=function()	
	local f, rfsProblem, liftOptions,resultFileName,polSet,liftedPolSet,delta,
		outputFileName, resultVariableName,outputOptions,macaulay2Path,
		W1,W2,W3,lambda,indeterminates,homogenVariable,solutionIdx;
	Read("rfsWrapper.gap");
	macaulay2Path :=  Null ;
	macaulay2Path :=  "/home/kroeker/bin/" ;
	
	checkMacaulayRFSPackagePresence(macaulay2Path);
	
	polSet := get43222ExamplePolSet();
	
	liftOptions := rec();
	
	liftOptions.decimalPrecision := 16;


	outputOptions:=rec();
	outputOptions.resultFileName:="43222LiftResult.gap";	
	outputOptions.resultVariableName := "liftedPolSet";
	# additional Parameter 'output style'

	
	performRFSPolSetLift(polSet,liftOptions,outputOptions,macaulay2Path);
	
	Read(outputOptions.resultFileName);# initialises variable "$resultVariableName", here "liftedPolSet"
	
	# liftedPolSet:=varName;
	
	indeterminates := IndeterminatesOfPolynomialRing (liftedPolSet.liftedPolynomialRing);
	homogenVariable:= indeterminates[2];

	for solutionIdx in [1..Size(liftedPolSet.polSetFunctionLists) ]do

		W1 := liftedPolSet.polSetFunctionLists[ solutionIdx ][1];	 
		W2 := liftedPolSet.polSetFunctionLists[ solutionIdx ][2] ;
		W3 := liftedPolSet.polSetFunctionLists[ solutionIdx ][3];
		lambda := liftedPolSet.scalingFactorLists[ solutionIdx  ][3-2];
	
		W1 := Value(W1,[ homogenVariable ],[1.0]);
		W2 := Value(W2,[ homogenVariable ],[1.0]);
		W3 := Value(W3,[ homogenVariable ],[1.0]);
	
		f := W2 - lambda*W1 - W3;

		delta := Sqrt(Sum(List(CoefficientsOfUnivariatePolynomial(f),Norm)));
		Print(delta);
		Print("\n");
	od;

	#f:=  (-V1*W1)/ W2 ;
	#m := IMGMachine(f);
	#perms := List([1..3],i->PermList(Output(m,i)));
	#goal := [ (1,5,11,6)(2,3)(4,10)(7,12)(8,13,9),
	#	 (1,7,13,2)(3,8)(4,5)(6,12)(9,11,10),
	#	 (1,3,9,4)(2,8)(5,10)(6,7)(11,13,12) ];;
	#	 
	#change := RepresentativeAction(SymmetricGroup(13),perms,goal,OnTuples);


end;


	
sandbox:=function()
	
	#fam:=RationalFunctionsFamily(FamilyObj(1));
	#PolynomialByExtRepNC(fam,[[1,1],2]);
end;




