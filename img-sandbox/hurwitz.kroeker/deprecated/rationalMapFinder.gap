
# "rationalMapFinder.gap" Hurwitz rational map finder.
# Documentation: coming soon; see for starters the EXAMPLE file.

# TODO: 
#  -clean code (partly done)
#  -documentation 


LoadPackage ("guava"); #DivisorsMultivariatePolynomial, DegreeMultivariatePolynomial
#DivisorsMultivariatePolynomial:  factor multivariate Polynomial!
LoadPackage ("float"); #Root finding and complex numbers

# TODO: where 'setFloats' should be called?
SetFloats(MPC);;

if not IsBound(Null) then
	BindGlobal("Null", MakeImmutable([]) );
fi;

rationalMapFinder@ := rec();

rationalMapFinder@.IsConsistent := NewProperty("IsConsistent", IsObject );
rationalMapFinder@.checkConsistency := NewOperation("checkConsistency", [IsObject] );


##########################################################################################################################
IsIntegerPartition@ := function(partition)
	local shapeCopy,entry;
	if not IsList(partition) then
		return false;
	fi;
	if Size(partition)=0 then
		return false;
	fi;
	for entry in partition do
		if not IsInt(entry) or  entry<1 then 
			return false;
		fi;
	od;
	return true;
end;

rationalMapFinder@.IsIntegerPartition := IsIntegerPartition@;
#rationalMapFinder@.IsIntegerPartition := NewProperty("IsIntegerPartition", IsList);
#InstallMethod( rationalMapFinder@.IsIntegerPartition , "IsIntegerPartition test",[IsList], IsIntegerPartition@ );
Unbind( IsIntegerPartition@ );

# rationalMapFinder@.IsIntegerPartition := IsIntegerPartitionFilter; # dangerous, because IsIntegerPartitionFilter is not immutable.



IsShape@ := function(shape)
	local shapeCopy,entry;
	if not rationalMapFinder@.IsIntegerPartition(shape) then
		return false;
	fi;
	shapeCopy := ShallowCopy(shape);
	Sort(shapeCopy);
	return (shape=Reversed(shapeCopy));
end;

rationalMapFinder@.IsShape := IsShape@;
#rationalMapFinder@.IsShape := NewProperty("IsShape",IsObject);
#InstallMethod( rationalMapFinder@.IsShape , "IsShape test",[IsObject], IsShape@ );
Unbind( IsShape@ );


# create a shape from a unsorted integer list.# Todo: allow only IntegerPartitions as parameter.
createShape@ := function(list)
	local shapeCopy,entry,preShape;
	if not rationalMapFinder@.IsIntegerPartition(list) then
		Error("parameter has to be a integer partition (a integer list with entries>1");
	fi;
	preShape :=ShallowCopy(list);
	Sort(preShape);
	return Reversed(preShape);
end;


rationalMapFinder@.createShape := createShape@;
#rationalMapFinder@.createShape := NewOperation("createShape", [IsObject] );
#InstallMethod( rationalMapFinder@.createShape , "createShape test",[IsObject], createShape@ );
Unbind( createShape@ );

##########################################################################################################################

IsShapeList@ := function (shapeList  )
	local entry,degree;
	if not IsList(shapeList) then
		Info( InfoWarning, 2, "parameter has to be a list" );
		return false;
	fi;
	if Size(shapeList)>0 then
	
		for entry in shapeList do
			if not rationalMapFinder@.IsShape(entry) then
				Info( InfoWarning, 2,"Shape list entries  has to be shapes");
				return false;
			fi;
		od;
		degree := Sum(shapeList[1]);
		for entry in shapeList do
			if not Sum(entry)=degree then
				Info( InfoWarning, 2,"IsShapeList: shapes expected to have same degree ");
				return false;
			fi;
		od;
	fi;
	return true;	
end;

rationalMapFinder@.IsShapeList :=  IsShapeList@;
#rationalMapFinder@.IsShapeList := NewProperty("IsShapeList", IsObject );
#InstallMethod( rationalMapFinder@.IsShapeList , "IsShapeList filter",[IsObject], IsShapeList@ );
Unbind( IsShapeList@ );


createShapeList@ := function(list)
	local shape,entry,shapeList;
	shapeList:=[];
	for entry in list do
		shape := rationalMapFinder@.createShape(entry);
		Append(shapeList,[shape]);
	od;
	Assert(0, rationalMapFinder@.IsShapeList(shapeList));
	return shapeList;
end;

rationalMapFinder@.createShapeList := createShapeList@;
#rationalMapFinder@.createShapeList := NewOperation("createShapeList", [IsList] );
#InstallMethod( rationalMapFinder@.createShapeList , "createShapeList  ",[IsList],createShapeList@);
Unbind( createShapeList@ );


#########################################################################################################################
IsScalingRelationList@ := function(scalingRelation)
	local elem,defaultError,part;
	if not IsList(scalingRelation) then 
		 Info( InfoWarning, 2,"rfsProblem is not a record");
		return false;
	fi;

	defaultError := "Error: scalingRelation entries must be pairs of real and imaginary complex number part ";	
	for elem in scalingRelation do
		if not  IsList(elem) then
			Info( InfoWarning, 2,defaultError);		 
			return false;
		fi;

		if not  Size(elem)=2 then
				Info( InfoWarning, 2,defaultError);
			return false;
		fi;
		for part in elem do	
			if  not  part in Rationals then
				Info( InfoWarning, 2,defaultError);	 
				return false;
			fi;
		od;
	od;
 	return true;
end;


rationalMapFinder@.IsScalingRelationList :=  IsScalingRelationList@;
#rationalMapFinder@.IsScalingRelationList := NewProperty("IsScalingRelationList", IsObject );
#InstallMethod( rationalMapFinder@.IsScalingRelationList , "IsScalingRelationList  ",[IsObject], IsScalingRelationList@);
Unbind( IsScalingRelationList@ );

##################################################################################################################

IsRFSProblem@ := function (rfsProblem)
	local keyList,key,rnames ;
	if not IsRecord(rfsProblem) then 
		Info( InfoWarning, 2,"IsRFSProblemFkt: rfsProblem is not a record");
		return false;
	fi;
	 keyList:=["shapeList","scalingRelationList"];
	    rnames:=    RecNames(rfsProblem);
	
	    for key in 	keyList do 
		if    Position(rnames,key)=false then
			Info( InfoWarning, 2,"IsRFSProblemFkt: key shapeList or scalingRelations missing\n");
		    	return false;
		fi;
	    od;
	    if not Length(keyList)=Length(rnames) then 
	    	return false;
	    fi;
	    if not rationalMapFinder@.IsShapeList(rfsProblem.shapeList) then
	    	return false;
	    fi;
	    if not rationalMapFinder@.IsScalingRelationList(rfsProblem.scalingRelationList) then
	    	return false;
	    fi;
	    if (Size(rfsProblem.scalingRelationList)+3 <> Size(rfsProblem.shapeList) ) then
			Info( InfoWarning, 2,"IsRFSProblemFkt: to many or too less scaling relations.\n");
	    	    	return false;
	   fi;
	   return true;
	    
end;

rationalMapFinder@.IsRFSProblem := IsRFSProblem@;
#rationalMapFinder@.IsRFSProblem := NewProperty("IsRFSProblem", IsObject );
#InstallMethod( rationalMapFinder@.IsRFSProblem , "IsRFSProblem  ",[IsObject], IsRFSProblem@ );
Unbind( IsRFSProblem@ );

###########################################################################################################

IsNull@ := function (object)
	return (IsList(object) and IsEmpty(object));
end;

rationalMapFinder@.IsNull := IsNull@;
#rationalMapFinder@.IsNull := NewProperty("IsNull", IsObject );
#InstallMethod( rationalMapFinder@.IsNull , "IsNull  ",[IsObject], IsNull@ );
Unbind( IsNull@ );


IsIntOrNull@ := function(object)
	return (rationalMapFinder@.IsNull(object) or IsInt(object));
end;

rationalMapFinder@.IsIntOrNull := IsIntOrNull@ ;
#rationalMapFinder@.IsIntOrNull := NewProperty("IsIntOrNull", IsObject );
#InstallMethod( rationalMapFinder@.IsIntOrNull , "IsIntOrNull  ",[IsObject],IsIntOrNull@);
Unbind( IsIntOrNull@ );

############################################################################################################

IsRFSOptionList@ := function (rfsOptionList)
	local keyList,key,rnames,errormsg ;
	if not IsRecord(rfsOptionList) then 
		Info( InfoWarning, 2,"rfsOptionList is not a record\n");
		return false;
	fi;
	# keyList:=["minChar","maxChar","softExampleLimit","parallelize", "resultFileName"];
	 keyList:=["minChar","maxChar","softExampleLimit","parallelize"];
	    rnames:=    RecNames(rfsOptionList);
	     for key in 	keyList do 
		if    Position(rnames,key)=false or Position(rnames,key)=fail then
			errormsg:= Concatenation("key " , key ," is missing \n");
			Info( InfoWarning, 2,errormsg);
		    	return false;
		fi;
	    od;
    	#if not IsString(rfsOptionList.resultFileName) then 
    	#	# todo: check if resultFileName can be created/deleted
	#	Info( InfoWarning, 2,"rfsOptionList.resultFileName must be an string!\n");
	#	return false;
	#fi;
	if not rationalMapFinder@.IsIntOrNull(rfsOptionList.minChar)   then 
		Info( InfoWarning, 2,"rfsOptionList.minChar must be an int or Null!\n");
		return false;
	fi;
	# TODO: eliminate dublicate code...
	if not rationalMapFinder@.IsIntOrNull(rfsOptionList.maxChar)   then 
		Info( InfoWarning, 2,"rfsOptionList.maxChar must be an int or Null!\n");
		return false;
	fi;
	
	if not rationalMapFinder@.IsIntOrNull(rfsOptionList.softExampleLimit)  then 
		Info( InfoWarning, 2,"rfsOptionList.softExampleLimit must be an int!\n");
		return false;
	fi; 
	

	if not IsBool(rfsOptionList.parallelize) then 
		Info( InfoWarning, 2,"rfsOptionList.parallelize must be a boolean!\n");
		return false;
	fi;
	return true;
end;

rationalMapFinder@.IsRFSOptionList :=  IsRFSOptionList@;
#rationalMapFinder@.IsRFSOptionList := NewProperty("IsRFSOptionList", IsObject );
#InstallMethod( rationalMapFinder@.IsRFSOptionList , "IsRFSOptionList  ",[IsObject],IsRFSOptionList@);
Unbind( IsRFSOptionList@ );


RFSOptionListIsConsistent@ := function(rfsOptionList)
	Assert(0,rationalMapFinder@.IsRFSOptionList(rfsOptionList));
	if rationalMapFinder@.IsNull(rfsOptionList.minChar) then 
		return false;
	fi;
	if not  rfsOptionList.minChar >1 then 
		Info( InfoWarning, 2,"rfsOptionList.minChar must be >1!");
		return false;
	fi;	
	if   not rationalMapFinder@.IsNull(rfsOptionList.maxChar)  then 
		if   (rfsOptionList.maxChar-rfsOptionList.minChar) <0 then
			Info( InfoWarning, 2,"rfsOptionList.maxChar < minChar");
			return false;
		fi;	
	else    # maxChar is null
		if rationalMapFinder@.IsNull(rfsOptionList.softExampleLimit) or not rfsOptionList.softExampleLimit>0 then 
			Info( InfoWarning, 2,"in case rfsOptionList.maxChar<rfsOptionList.minChar condition rfsOptionList.softExampleLimit>0 required!");
			return false;			
		fi; 
		
	fi;
	
	if  not rationalMapFinder@.IsNull(rfsOptionList.softExampleLimit) then 
		if (rfsOptionList.softExampleLimit<1) then 
			Info( InfoWarning, 2,"softExampleLimit has to be >=1 or Null (=[]) !");
			return false;
		fi;
	fi;
	
	return true;
end;

rationalMapFinder@.RFSOptionListIsConsistent := RFSOptionListIsConsistent@;
#InstallMethod( rationalMapFinder@.IsConsistent , "IsConsistent   ",[rationalMapFinder@.IsRFSOptionList], RFSOptionListIsConsistent@);
#rationalMapFinder@.RFSOptionListIsConsistent := NewOperation("RFSOptionListIsConsistent", [IsObject] );
#InstallMethod( rationalMapFinder@.RFSOptionListIsConsistent , "RFSOptionListIsConsistent   ",[IsObject], RFSOptionListIsConsistent@);
Unbind( RFSOptionListIsConsistent@ );


RFSOptionListCheckConsistency@ := function(rfsOptionList)
	Assert(0, rationalMapFinder@.RFSOptionListIsConsistent (rfsOptionList));
end;

rationalMapFinder@.RFSOptionListCheckConsistency := RFSOptionListCheckConsistency@;
#rationalMapFinder@.RFSOptionListCheckConsistency := NewOperation("RFSOptionListCheckConsistency", [IsObject] );
#InstallMethod( rationalMapFinder@.RFSOptionListCheckConsistency , "RFSOptionListCheckConsistency   ",[IsObject], RFSOptionListCheckConsistency@);

Unbind( RFSOptionListCheckConsistency@ );


createDefaultRFSOptionList@ := function()
	local rfsOptionList;
	rfsOptionList := rec (minChar:=2, maxChar:=[], softExampleLimit:=1,   parallelize:=false);
	Assert(0,rationalMapFinder@.IsRFSOptionList(rfsOptionList));
	return rfsOptionList;
end;

rationalMapFinder@.createDefaultRFSOptionList := createDefaultRFSOptionList@;
#rationalMapFinder@.createDefaultRFSOptionList := NewOperation("createDefaultRFSOptionList", [] );
#InstallMethod( rationalMapFinder@.createDefaultRFSOptionList , "createDefaultRFSOptionList   ",[], createDefaultRFSOptionList@);
Unbind( createDefaultRFSOptionList@ );

createRFSOptionListByCharFkt@ := function(characteristic)
	local rfsOptionList;
	Assert(0, rationalMapFinder@.IsPrime(characteristic));
	rfsOptionList := rec (minChar:=characteristic, maxChar:=characteristic, softExampleLimit:=[],  parallelize:=false);
	Assert(0, rationalMapFinder@.IsRFSOptionList(rfsOptionList));
	return rfsOptionList;
end;

rationalMapFinder@.IsPrime := NewProperty("IsPrime", IsInt );

InstallMethod( rationalMapFinder@.IsPrime, "IsPrime", [IsInt], 
function(primecandidate)
	 return IsPrime(primecandidate);
end );

rationalMapFinder@.createRFSOptionListByChar := NewOperation("createRFSOptionListByChar", [IsPosInt] );
InstallMethod( rationalMapFinder@.createRFSOptionListByChar , "createRFSOptionListByChar   ",[IsPosInt],createRFSOptionListByCharFkt@);
Unbind( createRFSOptionListByCharFkt@ );

##################################################################################################################################



IsRFSPolynomialSetFkt@ := function(polSet)
    local indeterminates, rng, rnames, key, keyList, pol,smallDegreeCount;
    if not  IsRecord(polSet) then
       return false;
    fi;
    # check if ring and polynomialTuple are elements of RecNames
    keyList:=["polynomialRing","degree","polynomialTuple"];
    rnames:=    RecNames(polSet);
    if not Length(keyList)=Length(rnames) then 
    	return false;
    fi;
    for key in 	keyList do 
        if    Position(rnames,key)=false then
	    	return false;
	fi;
    od;
    if not (IsPolynomialRing(polSet.polynomialRing) ) then
    	return false;
    fi;
    
    indeterminates := IndeterminatesOfPolynomialRing(polSet.polynomialRing);
     if   Size(indeterminates)<>2 then
     	return false;
     fi;	
    if   Length(polSet.polynomialTuple)<3 then 
    	return false;
    fi;
    if not IsPosInt(polSet.degree) then 
    	return false; 
    fi;
    smallDegreeCount:=0;
    for pol in polSet.polynomialTuple do
       	if not IN(pol,polSet.polynomialRing) then 
       		return false;
       	fi;
    	#if not IsUnivariatePolynomial(pol) then 
    	#	return false;
    	#fi;
    	#TODO: DegreeMultivariatePolynomial temporarily removed
    	#if DegreeMultivariatePolynomial(pol,polSet.polynomialRing)<>polSet.degree then 
    	#	return false;
    	#fi;
	#if  DegreeMultivariatePolynomial(pol,polSet.polynomialRing)<>polSet.degree then 
	#	return false;
    	#	smallDegreeCount:=smallDegreeCount+1;
    	#fi;
    od;
    if smallDegreeCount>1 then 
    	return false;
    fi;
    #todo: check GCD = 1 ?
    return true;
end;



rationalMapFinder@.IsRFSPolynomialSet := NewProperty("IsRFSPolynomialSet", IsObject );
InstallMethod( rationalMapFinder@.IsRFSPolynomialSet , "IsRFSPolynomialSet   ",[IsObject],IsRFSPolynomialSetFkt@);
Unbind( IsRFSPolynomialSetFkt@ );
 
 
IsPolynomialListFkt@ := function(polList)
	local entry;
	for entry in polList do
		if not IsPolynomial(entry) then
			return false;
		fi;
	od;
	return true;
end;

 
rationalMapFinder@.IsPolynomialList := NewProperty("IsPolynomialList", IsList );
InstallMethod( rationalMapFinder@.IsPolynomialList , "IsPolynomialList   ",[IsList],IsPolynomialListFkt@);
Unbind( IsPolynomialListFkt@ );

 
createRFSPolynomialSetFkt@ := function(parPolynomialList,rng)
	local polSetDegree,polSet;
	Assert(0, rationalMapFinder@.IsPolynomialList(parPolynomialList) );
	if Size(parPolynomialList)<3 then
		Error("polynomialTuple shorter than 3");
	fi;
	
	#TODO temporarily removed polSetDegree:= DegreeMultivariatePolynomial(parPolynomialList[1],polSet.polynomialRing);
	
	polSetDegree:=[];
	polSet := rec( polynomialRing:=rng, degree:= polSetDegree, polynomialTuple := parPolynomialList );
	#SetInfoLevel(InfoWarning,2); 
	if not rationalMapFinder@.IsRFSPolynomialSet(polSet)  then

		Error("createRFSPolynomialSet failed, use higher InfoWarning level for more info.");
	fi;
	return polSet;
end;

rationalMapFinder@.createRFSPolynomialSet := NewOperation("createRFSPolynomialSet", [ IsObject,  IsPolynomialRing] );
InstallMethod( rationalMapFinder@.createRFSPolynomialSet , "createRFSPolynomialSet ",[ IsObject, IsPolynomialRing ], createRFSPolynomialSetFkt@);
Unbind( createRFSPolynomialSetFkt@ );


#####################################################################################################################
## TODO: move examples to documentation!

get43222RFSProblemFkt := function()
	local shape, shapeList,scalingRelationList,rfsProblem;
	shape := rationalMapFinder@.createShape([4,3,2,2,2]);
	shapeList := [ shape,shape,shape ];
	scalingRelationList := [];
	rfsProblem := rec(shapeList:=shapeList ,  scalingRelationList:= scalingRelationList);
	Assert(0,rationalMapFinder@.IsRFSProblem(rfsProblem));
	return rfsProblem;
end;



get43222Char11ExamplePolSetFkt := function()
	local shape,shapeList,scalingRelationList,rfsProblem,polRing,ind,polSet,s,t,nicePolRing;
	shape := rationalMapFinder@.createShape([4,3,2,2,2]);
	shapeList := [ shape,shape,shape ];
	scalingRelationList:=[];
	rfsProblem := rec(shapeList:=shapeList ,  scalingRelationList:= scalingRelationList);
	Assert(0,rationalMapFinder@.IsRFSProblem(rfsProblem));
	
	polSet:=rec();
	 polRing := PolynomialRing(  Field( Z(11)) ,["t","s"] : new );
	
	  ind:=IndeterminatesOfPolynomialRing(polRing);
	 t:=ind[1];
	 s:=ind[2];
	 polSet:=rec();
	 polSet.rfsProblem := rfsProblem; 
	 polSet.polynomialRing := polRing; 
	 polSet.polynomialTuple:= [	 (s)^4  *( t -5*s )^3*( t^3 +3*t^2*s +2*t*s^2 +3*s^3 )^2 ,
					 (t)^4  *( t +3*s )^3*( t^3 -3*t*s^2 -5*s^3 )^2,
					 (t-s)^4*( t -3*s )^3*( t^3 -2*t*s^2 -3*s^3 )^2];

	 nicePolRing := PolynomialRing(  Rationals ,["t","s"] : new);
	 polSet.ZZPolRing:=nicePolRing;					 
         ind := IndeterminatesOfPolynomialRing(nicePolRing);
	 t := ind[1];
	 s := ind[2];

	  polSet.polynomialTupleZZLift:= [ (s)^4  *( t -5*s )^3*( t^3 +3*t^2*s +2*t*s^2 +3*s^3 )^2 ,
					 (t)^4  *( t +3*s )^3*( t^3 -3*t*s^2 -5*s^3 )^2,
					 (t-s)^4*( t -3*s )^3*( t^3 -2*t*s^2 -3*s^3 )^2];

	Assert(0, rationalMapFinder@.IsRFSPolynomialSet(polSet));
	return polSet;
end;
 


# creates polynomials [A,B,C,...] from single rootData with  B-lambdaA = C, B-mueA = D, etc.
createPolynomialListFromRootDataElemFkt@ := function ( preimageList, polynomialRing )
	local currentPolynomial,polynomialList, pos,ind, preimageData ;
	polynomialList := [];
	 
	ind := IndeterminatesOfPolynomialRing(polynomialRing);
	for pos in [1..Size(preimageList)] do
		currentPolynomial:=1.0;
		for preimageData in preimageList[pos] do
			if (preimageData[1]<>infinity) then 
				currentPolynomial := currentPolynomial*( ( ind[1] - preimageData[1] )^preimageData[2] );
			#else
				#currentPolynomial=currentPolynomial 
			fi;
			
		od;
		Append(polynomialList,[currentPolynomial]);
	od;
	return polynomialList; 
end;

rationalMapFinder@.createPolynomialListFromRootDataElem := NewOperation("createPolynomialListFromRootDataElem", [ IsObject,IsPolynomialRing] );
InstallMethod( rationalMapFinder@.createPolynomialListFromRootDataElem , "createPolynomialListFromRootDataElem ",[ IsObject,IsPolynomialRing ], createPolynomialListFromRootDataElemFkt@);
Unbind( createPolynomialListFromRootDataElemFkt@ );

preRationalMapFromRootDataElemFkt@ := function ( rootDataElement, polynomialRing )
	local polynomialList,rationalMapList, ind,currPos,scalingFactor, num , denom ;
	polynomialList := rationalMapFinder@.createPolynomialListFromRootDataElem( rootDataElement.preimageLists, polynomialRing ) ;
	rationalMapList :=[];
	ind := IndeterminatesOfPolynomialRing(polynomialRing);
	currPos := 3;
	for scalingFactor in rootDataElement.scalingFactorList do
		num := polynomialList[2];
		denom := polynomialList[1]*scalingFactor;
		Append(rationalMapList,[ rec( numerator:=num , denominator:=denom ) ] );
		currPos := currPos+1;
	od;
	return rationalMapList;
end;

rationalMapFinder@.preRationalMapFromRootDataElem := NewOperation("preRationalMapFromRootDataElem", [ IsObject,IsPolynomialRing] );
InstallMethod( rationalMapFinder@.preRationalMapFromRootDataElem , "preRationalMapFromRootDataElem ",[ IsObject,IsPolynomialRing ], preRationalMapFromRootDataElemFkt@);
Unbind( preRationalMapFromRootDataElemFkt@ );


M2EXEC@ := Filename( DirectoriesSystemPrograms() , "M2");
	
ReadPackage("fr","hurwitz/MacaulayMapFinderWrapper.gap");
ReadPackage("fr","hurwitz/FunctionalStyleInterface.gap");
#ReadPackage("fr","hurwitz/AllOperations.gap");


rationalMapFinder@ := Immutable(rationalMapFinder@);
	




