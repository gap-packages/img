
#####################################################################################################################


## todo: newOperation nur fuer die Schnittstellen, nicht fuer private Funktionen...

checkMacaulayPresence@ := function( M2Command )
	Assert(0, IsString( M2Command ) or IsBool( M2Command ) );
  	if M2Command=fail then 
  		Error ("could not find Macaulay2 binary. Please install Macaulay2 on your system and/or setup the environment appropriately to use the 'rationalFunctionSearch' package");
  	else
  		Info( InfoWarning, 1,"M2 was found!");
	fi;
end;


#rationalMapFinder@.checkMacaulayPresence := checkMacaulayPresence@;
rationalMapFinder@.checkMacaulayPresence := NewOperation("checkMacaulayPresence",[IsObject]);
InstallMethod( rationalMapFinder@.checkMacaulayPresence , "checkMacaulayPresence test", 	[IsObject], checkMacaulayPresence@);
Unbind( checkMacaulayPresence@ );
	

checkMacaulayRFSPackagePresence@ := function(  M2Command   )
	local M2InputStringStream, tmpdir, str, outStrStream, strerror, strSuccess, commandString;
	Assert(0, IsString( M2Command ) or IsBool( M2Command ) );
	#rationalMapFinder@.checkMacaulayPresence( M2Command ); 
 	 
	commandString :=        "try (loadPackage(\"rationalFunctionSearch\")\n ) then (print (\"loadsuc\" | \"cess\"))"; 
	Append (commandString, "else (print (\"loader\" | \"ror\") ); ");
 	M2InputStringStream := InputTextString(commandString);

 	tmpdir := DirectoryTemporary();;
 	str := "";; 
 	outStrStream := OutputTextString( str, false );;

	Process( tmpdir, M2Command, M2InputStringStream , outStrStream, [ "--no-prompts"  ] );;
	strerror := ReplacedString(str,"loaderror","");
	strSuccess := ReplacedString(str,"loadsuccess","");
	if ( Size(strerror)<Size(str) ) then 
		Error(" Macaulay2: load Macaulay2 package 'rationalFunctionSearch' failed. Please install the package.");
		# Check if you are in the directory of the package or the Macaulay2 environment is set (correct path in '~/.Macaulay2/init.m2' 
	else
		if (Size(strSuccess)<Size(str)) then 
			Info( InfoWarning, 1,"Macaulay2 package 'rationalFunctionSearch' available!\n");
		else
			Error("checkMacaulayRFSPackagePresence failed...");
		fi;
	fi;
end;

#rationalMapFinder@.checkMacaulayRFSPackagePresence := checkMacaulayRFSPackagePresence@;
rationalMapFinder@.checkMacaulayRFSPackagePresence := NewOperation( "checkMacaulayRFSPackagePresence", [IsObject] );
InstallMethod( rationalMapFinder@.checkMacaulayRFSPackagePresence , "checkMacaulayRFSPackagePresence test",
		[IsObject], checkMacaulayRFSPackagePresence@);
Unbind( checkMacaulayRFSPackagePresence@ );



createM2RFSSearchString@ := 	function  ( rfsProblem, searchOptions, resultDestFileName)
	local   commandString, key;
	
	Assert(0, rationalMapFinder@.IsRFSProblem(rfsProblem) );
	Assert(0, rationalMapFinder@.IsRFSOptionList(searchOptions) );
			
	commandString:= ""; 
	#str := OutputTextString( commandString, true );
	Append (commandString, "loadPackage(\"rationalFunctionSearch\")\n");
	Append (commandString, "shapeArray:= ");
	Append (commandString,String(rfsProblem.shapeList));
	Append (commandString, "\n");
	Append (commandString,"shapeList:= new RFSShapeList from (apply(shapeArray,el->new Shape from el))\n");

	Append (commandString, "scalingArray:= ");
	Append (commandString,String(rfsProblem.scalingRelationList));
	Append (commandString, "\n");
	Append (commandString,"scalingRelationList:= new List from  (apply(scalingArray,el->new List from el))\n");
	# todo: ohne if then else auskomen...?
	if Size(rfsProblem.scalingRelationList)>0 then 
		Append (commandString, "scalingRelationObj := createScalingRelations( scalingRelationList,true  );\n");
	else
		Append (commandString, "scalingRelationObj := createScalingRelations( null );\n");
	fi;

	if (IsBound(rfsProblem.normalizedFactorDegrees)) then 
		Append (commandString, "strictNormRuleSet := createStrictNormalizationRuleSet(    ");
		Append( commandString,String(rfsProblem.normalizedFactorDegrees ) );
		Append (commandString, ");\n");
		Append (commandString, "rfsProblem := createRFSProblem( shapeList, strictNormRuleSet, scalingRelationObj  );\n");
	else
		Append (commandString, "rfsProblem := createRFSProblem( shapeList, scalingRelationObj  );\n");
	
	fi;
	Append (commandString," constructOptions := createPolynomialConstructOptions(null \n");
	for key in RecNames(searchOptions) do 
		if not rationalMapFinder@.IsNull(searchOptions.(key)) then 
			Append (commandString,",\n \"");
		 
			Append(commandString,String(key));
			Append (commandString,"\" => ");
			 
			Append(commandString,String(searchOptions.(key)));
			Append (commandString," \n ");
		fi;
	
	od;
	Append (commandString,")\n");
	Append (commandString,"sortedPolSetList := findFiniteFieldRFSExamples(rfsProblem, constructOptions );\n");
	Append(commandString,"saveFFRFSResultInGAPfile(rfsProblem,constructOptions, sortedPolSetList, ");
	Append(commandString,"\"");
	Append(commandString, resultDestFileName );
	Append(commandString,"\" ");
	Append(commandString,");\n");
	return commandString;
end;
	
#rationalMapFinder@.createM2RFSSearchString := createM2RFSSearchString@;
rationalMapFinder@.createM2RFSSearchString := NewOperation( "createM2RFSSearchString", [ IsRecord, IsRecord, IsString ] );
InstallMethod( rationalMapFinder@.createM2RFSSearchString , "createM2RFSSearchString ",	[ IsRecord, IsRecord, IsString ], createM2RFSSearchString@);
Unbind( createM2RFSSearchString@ );
	

performM2RFSSearchOverFiniteFields@ :=	function ( M2HurwitzMapFinder, rfsProblem, searchOptions )
	local  runDir,  commandString, outStr, outStrStream,  output, outputFileName, inputStream, rfsObj, searchOutputOptions;


	searchOutputOptions := M2HurwitzMapFinder.getSearchOutputOptions();
	commandString :=  M2HurwitzMapFinder.private.createM2RFSSearchString( rfsProblem, searchOptions, searchOutputOptions.resultFileName );
     
	runDir:= searchOutputOptions.outputDirectory;
 	outStr := ""; outStrStream := OutputTextString(outStr,false);;
 	 
	Info( InfoWarning, 2, commandString); 
 	
 	outputFileName := Concatenation( searchOutputOptions.resultFileName, ".M2input" );
 	output := OutputTextFile( outputFileName , true );;
 	WriteAll( output, commandString );
 		 
	Info( InfoWarning, 2,commandString);
	CloseStream(output);
	
 	#inputStream := InputTextFile(outputFileName);
 	inputStream := InputTextString( "" );
 	Process(runDir, M2HurwitzMapFinder.getM2Command(), inputStream , outStrStream, [  String(outputFileName), "--no-prompts"  ] );;
	rfsObj := ReadAsFunction( searchOutputOptions.resultFileName )(); 	#sets variable $resultVariableName  (here "rfsObj") !  
	if (not M2HurwitzMapFinder.debug) then 
		RemoveFile(outputFileName);
		RemoveFile( searchOutputOptions.resultFileName );
	fi;
	
	return rfsObj;
end;

#rationalMapFinder@.performM2RFSSearchOverFiniteFields := performM2RFSSearchOverFiniteFields@;
rationalMapFinder@.performM2RFSSearchOverFiniteFields := NewOperation( "performM2RFSSearchOverFiniteFields", [ IsRecord, IsRecord, IsRecord ] );
InstallMethod( rationalMapFinder@.performM2RFSSearchOverFiniteFields , "performM2RFSSearchOverFiniteFields ",
					[ IsRecord, IsRecord, IsRecord ], performM2RFSSearchOverFiniteFields@);
Unbind( performM2RFSSearchOverFiniteFields@  );


createM2PolSetLiftCommandString@ := function( polSet, liftOptions, resultFileName)
	local rfsProblem, commandString, key, ind ;
	
	commandString:= ""; 
	rfsProblem := polSet.rfsProblem;

	#str := OutputTextString(commandString,true);
	Append (commandString, "loadPackage(\"padicLift\")\n");
	Append (commandString, "loadPackage(\"rationalFunctionSearch\")\n");
	Append (commandString, "shapeArray:= ");
	Append(commandString,String(rfsProblem.shapeList));
	Append (commandString, "\n");
	Append (commandString,"shapeList:= new RFSShapeList from (apply(shapeArray,el->new Shape from el))\n");

	Append (commandString, "scalingArray:= ");
	Append (commandString,String(rfsProblem.scalingRelationList));
	Append (commandString, "\n");
	Append (commandString,"scalingRelationList:= new List from  (apply(scalingArray,el->new List from el))\n");
	# todo: ohne if then else auskomen...?
	if Size(rfsProblem.scalingRelationList)>0 then 
		Append (commandString, "scalingRelationObj := createScalingRelations( scalingRelationList,true  );\n");
	else
		Append (commandString, "scalingRelationObj := createScalingRelations( null );\n");
	fi;
	#todo: get rid if duplicate code...
	if ( IsBound(rfsProblem.normalizedFactorDegrees) ) then 
		Append (commandString, "strictNormRuleSet := createStrictNormalizationRuleSet(    ");
		Append ( commandString,String(rfsProblem.normalizedFactorDegrees ) );
		Append (commandString, ");\n");
		Append (commandString, "rfsProblem := createRFSProblem( shapeList, strictNormRuleSet, scalingRelationObj  );\n");
	else
		Append (commandString, "rfsProblem := createRFSProblem( shapeList, scalingRelationObj  );\n");
	
	fi;

	Append (commandString, "  rng:=createRFSRing( ");
		Append(commandString, String(Characteristic(polSet.polynomialRing)) );
	Append (commandString, "  ); \n");
	ind := IndeterminatesOfPolynomialRing(polSet.polynomialRing);
	Assert(0, Size(ind)=2 );
	Append( commandString, String(ind[1]) );
	Append (commandString, ":=  commonVariable rng; \n");
	Append( commandString, String(ind[2]) );
	Append (commandString, ":= homogenVariable rng; \n");
	#todo: check in Macaulay if polynomials are homogenized. (when creating createRFSPolSet)

	Append (commandString, "polArray := ");
	Append (commandString, String(polSet.polynomialTupleZZLift));
	Append (commandString, "; \n ");

	Append (commandString, "polSet := createRFSPolynomialSet( rfsProblem, polArray);\n ");

	Append (commandString," liftOptions := createLiftOptions(null ");
	for key in RecNames(liftOptions) do 
		if not rationalMapFinder@.IsNull(liftOptions.(key)) then 
			Append (commandString,",\n \"");
			Append(commandString,String(key));
			Append (commandString,"\" => ");
			Append(commandString,String(liftOptions.(key)));
			Append (commandString,"  ");
		fi;
	
	od;
	Append (commandString," );\n");

	Append (commandString,"  tryLiftAndLLLAndPairPolSet( polSet, \"liftAndLLLOptions\"=>liftOptions ); \n");
	     Append(commandString,"saveRFSLiftResultInGAPfile(polSet,liftOptions, ");
	Append(commandString,"\"");
	Append( commandString, String(resultFileName) );
	Append(commandString,"\"");
	Append(commandString,");\n");	
	commandString:=ReplacedString(commandString,"\\",""); # hab schon vergessen, wozu diese Ersetzung notwendig ist.
	return   commandString;	
end;

#rationalMapFinder@.createM2PolSetLiftCommandString := createM2PolSetLiftCommandString@;
rationalMapFinder@.createM2PolSetLiftCommandString := NewOperation( "createM2PolSetLiftCommandString", [ IsRecord, IsRecord, IsString ] );
InstallMethod( rationalMapFinder@.createM2PolSetLiftCommandString , "createM2PolSetLiftCommandString ",
		 [ IsRecord, IsRecord, IsString ], createM2PolSetLiftCommandString@);
Unbind( createM2PolSetLiftCommandString@  );


performM2RFSPolSetLift@ := function ( M2HurwitzMapFinder, polSet, liftOptions  )
	local runDir, commandString,outStr, outStrStream, output, outputFileName, inputStream, liftedPolSet,liftOutputOptions ;

	liftOutputOptions := M2HurwitzMapFinder.getLiftOutputOptions();

	commandString := M2HurwitzMapFinder.private.createM2PolSetLiftCommandString(polSet, liftOptions, liftOutputOptions.resultFileName);
	#Print(commandString); # -- debug
	 
	runDir:= liftOutputOptions.outputDirectory;
 	outStr := ""; outStrStream := OutputTextString( outStr,false );; 	
 	
 	# reason for saving and loading: Macaulay has problems with "\ cr" in the InputTextString(commandString).
 	# Writing the commandString to a file and then creating  the input stream via InputTextFile works.
 	outputFileName := Concatenation( liftOutputOptions.resultFileName, ".M2input" );
 	Print("outputFileName:");
 	Print(outputFileName);Print("\n");
 	output := OutputTextFile( outputFileName , true );;
 	 
 	WriteAll( output, commandString );
	CloseStream(output);
 	#inputStream := InputTextFile(outputFileName);
	inputStream := InputTextString( "" );
 
	Info( InfoWarning, 2,commandString);
	
 	Process(runDir, M2HurwitzMapFinder.getM2Command(), inputStream , outStrStream, [ String(outputFileName), "--no-prompts"  ] );;
 	
 	if (not M2HurwitzMapFinder.debug) then 
		RemoveFile(outputFileName);
	fi;

 	liftedPolSet := ReadAsFunction( liftOutputOptions.resultFileName)(); 	#sets variable $resultVariableName  (here "rfsObj") !  
 	
 	if (not M2HurwitzMapFinder.debug) then 
		RemoveFile( liftOutputOptions.resultFileName);
	fi;

	return liftedPolSet;
end;

rationalMapFinder@.performM2RFSPolSetLift := performM2RFSPolSetLift@;	
rationalMapFinder@.performM2RFSPolSetLift := NewOperation( "performM2RFSPolSetLift", [ IsRecord, IsRecord, IsRecord ] );
InstallMethod( rationalMapFinder@.performM2RFSPolSetLift , "performM2RFSPolSetLift  ",
		 [ IsRecord, IsRecord, IsRecord ],performM2RFSPolSetLift@);	
Unbind( performM2RFSPolSetLift@  );

###########################################################################


performM2RFSSearchOverFiniteFieldsFI@ := function( M2HurwitzMapFinder, integerPartitionList, strictNormalization, scalingRelationList, finiteFieldChar )
	
	local key, shapeList, rfsProblem, searchOptions, rfsObj, polSetList;

	Assert( 0, IsPrime( finiteFieldChar ) );
	shapeList  := rationalMapFinder@.createShapeList( integerPartitionList );
	rfsProblem := rec( shapeList:=shapeList ,  scalingRelationList := scalingRelationList);
	if (strictNormalization) then 
		Assert(0, Size(integerPartitionList)>2 );
		rfsProblem.normalizedFactorDegrees := [ integerPartitionList[1][1], integerPartitionList[2][1], integerPartitionList[3][1] ];
	fi;

	searchOptions := rationalMapFinder@.createRFSOptionListByChar(finiteFieldChar);	
	rationalMapFinder@.RFSOptionListCheckConsistency(searchOptions);		

	rfsObj := M2HurwitzMapFinder.computeFiniteFieldSolutions( rfsProblem, searchOptions );

	polSetList := [];
	# put all example candidates in one list
	for key in RecNames(rfsObj.polynomialSetTable) do
		Append(polSetList,rfsObj.polynomialSetTable.(key) );
	od;

	return polSetList;
end;

rationalMapFinder@.performM2RFSSearchOverFiniteFieldsFI := performM2RFSSearchOverFiniteFieldsFI@;
rationalMapFinder@.performM2RFSSearchOverFiniteFieldsFI := NewOperation( "performM2RFSSearchOverFiniteFieldsFI", [ IsRecord, IsList, IsBool, IsList, IsPosInt ] );
InstallMethod( rationalMapFinder@.performM2RFSSearchOverFiniteFieldsFI , "performM2RFSSearchOverFiniteFieldsFI test",
		 [ IsRecord, IsList, IsBool, IsList, IsPosInt ], performM2RFSSearchOverFiniteFieldsFI@);
Unbind( performM2RFSSearchOverFiniteFieldsFI@  );


findBinary@ := function (binaryName, searchPathList)
	local  localSearchPathList;
	localSearchPathList := searchPathList;
	if rationalMapFinder@.IsNull( localSearchPathList ) then
		 localSearchPathList :=  DirectoriesSystemPrograms();
	fi;
	return Filename( localSearchPathList, binaryName );
end;

rationalMapFinder@.findBinary := findBinary@;
Unbind( findBinary@ );


createMacaulay2HurwitzMapFinder@ := function ( M2BinaryFileName )
	local M2HurwitzMapFinder,private   ;
		
	M2HurwitzMapFinder:= rec();
	
	
	#################### private : ####################################
	private := rec();
	private.checkMacaulayPresence := rationalMapFinder@.checkMacaulayPresence;
	#Unbind( rationalFunctionSearch.checkMacaulayPresence );
	private.checkMacaulayRFSPackagePresence := rationalMapFinder@.checkMacaulayRFSPackagePresence;
	#Unbind( rationalFunctionSearch.checkMacaulayRFSPackagePresence );
	
	private.createM2RFSSearchString := rationalMapFinder@.createM2RFSSearchString;
	#Unbind( rationalFunctionSearch.createM2RFSSearchString );
	
	private.createM2PolSetLiftCommandString := rationalMapFinder@.createM2PolSetLiftCommandString;
	#Unbind( rationalFunctionSearch.createM2PolSetLiftCommandStrin );
	 
		
	private.performM2RFSSearchOverFiniteFieldsFI := rationalMapFinder@.performM2RFSSearchOverFiniteFieldsFI;
	#Unbind( rationalFunctionSearch.performM2RFSSearchOverFiniteFieldsFI );
	
	private.performM2RFSSearchOverFiniteFields := rationalMapFinder@.performM2RFSSearchOverFiniteFields;
	#Unbind( rationalFunctionSearch.performM2RFSSearchOverFiniteFields );
	
	private.performM2RFSPolSetLift := rationalMapFinder@.performM2RFSPolSetLift;
	
	M2HurwitzMapFinder.private := Immutable(private);
	################### public #####################################	
		
	M2HurwitzMapFinder.computeFiniteFieldSolutionsFI := function( integerPartitionList, strictNormalization, scalingRelationList, finiteFieldChar  )
	
		return  M2HurwitzMapFinder.private.performM2RFSSearchOverFiniteFieldsFI( M2HurwitzMapFinder , 
												 integerPartitionList,
												 strictNormalization, 
												 scalingRelationList, 
												 finiteFieldChar);
	end;
	
	
	M2HurwitzMapFinder.computeFiniteFieldSolutions := function( rfsProblem, searchOptions  )
		return  M2HurwitzMapFinder.private.performM2RFSSearchOverFiniteFields( M2HurwitzMapFinder ,  rfsProblem, searchOptions );
	end;
	

	M2HurwitzMapFinder.approximateComplexSolutions			:= function(  polSet, liftOptions)
		return M2HurwitzMapFinder.private.performM2RFSPolSetLift( M2HurwitzMapFinder,  polSet, liftOptions);
	end;
	 

	M2HurwitzMapFinder.getM2Command := function () 
		return M2BinaryFileName;
	end;
	
	
	M2HurwitzMapFinder.checkMacaulayPresence:= function() 
		private.checkMacaulayPresence( M2HurwitzMapFinder.getM2Command() );
	end;
	
	M2HurwitzMapFinder.checkMacaulayRFSPackagePresence :=   function() 
		private.checkMacaulayRFSPackagePresence( M2HurwitzMapFinder.getM2Command() );
	end;
	
	
	
	##########################################
	M2HurwitzMapFinder.getSearchOutputOptions := function()
		local searchOutputOptions,randInt;
		searchOutputOptions := rec();
		randInt:= Random(10000000,100000000000);
		searchOutputOptions.resultFileName := Concatenation("RFSsearchResult",String(randInt),".gap");
	
		searchOutputOptions.resultVariableName := Concatenation("rfsObj",String(randInt)); 
	
		searchOutputOptions.outputDirectory := DirectoryTemporary();;
	
		searchOutputOptions.resultFileName := Filename( searchOutputOptions.outputDirectory , searchOutputOptions.resultFileName );
		return Immutable(searchOutputOptions);
	end;
	##########################################
	 
	M2HurwitzMapFinder.getLiftOutputOptions := function()
		local liftOutputOptions, randInt;
		liftOutputOptions := rec();
		randInt:= Random(10000000,100000000000);
		liftOutputOptions.resultFileName := Concatenation("RFSLiftResult",String(randInt),".gap");
	
		liftOutputOptions.resultVariableName := Concatenation("liftedPolSet",String(randInt)); 
		liftOutputOptions.outputDirectory := DirectoryTemporary();;	
		liftOutputOptions.resultFileName := Filename( liftOutputOptions.outputDirectory , liftOutputOptions.resultFileName );
		return Immutable(liftOutputOptions);
	end;
	##########################################
	
	########### checks#######################
	M2HurwitzMapFinder.checkMacaulayPresence();
	M2HurwitzMapFinder.checkMacaulayRFSPackagePresence();
	
	M2HurwitzMapFinder.debug := false;
	
	return Immutable(M2HurwitzMapFinder);
end;

rationalMapFinder@.createMacaulay2HurwitzMapFinder := createMacaulay2HurwitzMapFinder@ ;
rationalMapFinder@.createMacaulay2HurwitzMapFinder := NewOperation("createMacaulay2HurwitzMapFinder" ,[IsString] );
InstallMethod(rationalMapFinder@.createMacaulay2HurwitzMapFinder, " creates a HurwitzMap finder using Macaulay2 routines", [IsString],createMacaulay2HurwitzMapFinder@);

 



