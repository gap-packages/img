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

# dependencies: package 'float'!


# TODO: GF(11) vs GF(121) vs ZmodnZ(121) : what are the type of its elements  and how to check for it?


# todo: update documentation...
# adds a function to a global record structure described by a list 'recordstructure' of nesting record names.
# Example: 
# if recordstructure[1] is a global record structure, 
# installs recordstructure[1].recordstructure[2]. ... recordstructure[i].functionName := functionRef
# in case functionName is a new record entry.

 # todo: get rid of duplicate code (see also InstallGlobalRecordFunctionEx) 
 InstallGlobalFunction( InstallGlobalRecordFunctionOrMethod@FR,
function( recordstructure, functionName, comments, params, functionRef, installFkt, reinstall)

    local headRec, currentRec, name, i,fullName, globalName;
    
    Assert(0, Size( recordstructure)>0);
      
    Assert(0, IsBoundGlobal( recordstructure[1]) );

    fullName := Concatenation(recordstructure[1],"\.");
    headRec := ShallowCopy( ValueGlobal( recordstructure[1] ) );
    currentRec := headRec;
    Assert(0, IsRecord(currentRec) );
    for i in [2..Size(recordstructure)] do
        name := recordstructure[i];
        Assert(0, name in RecNames( currentRec ) );
        currentRec.(name) := ShallowCopy(currentRec.(name));
        currentRec := currentRec.(name);
         Assert(0, IsRecord(currentRec) );
        fullName := Concatenation(fullName,name,"\.");
    od;
   
    if functionName in RecNames(currentRec) then 
        if reinstall then
        
             if (installFkt=InstallMethod or installFkt=InstallOtherMethod) then 
                Error("cannot reinstall method for an operation");
             else
                currentRec.(functionName) := functionRef;
             fi;
         
        else # not reinstall
            if not IsOperation( currentRec.(functionName) ) then 
               Error(Concatenation( "function '", functionName, "' is already installed!" ));
            else       #  is operation 

              if not  (installFkt=InstallMethod or installFkt=InstallOtherMethod) then
                Error("trying overwrite an operation");
              else
                installFkt( currentRec.(functionName), comments, params, functionRef );
              fi;
            fi;
        fi;
    fi;
     if not functionName in RecNames(currentRec) then
       if (installFkt=InstallMethod or installFkt=InstallOtherMethod) then
            Error(Concatenation( "operation '", functionName, "' is not installed!" ));
        else
            currentRec.(functionName) := functionRef;
        fi;
    fi;
    
    fullName := Concatenation( fullName, functionName );
    #if IS_READ_ONLY_GLOBAL(fullName) then 
    #    MakeReadWriteGlobal( fullName );
    #fi;
    #if IsBoundGlobal(fullName) then 
    #    UnbindGlobal(  fullName ) ;    
    #fi;
    #
    #BindGlobal( fullName, currentRec.(functionName) ) ;
    
    #MakeReadOnlyGlobal(fullName);
    globalName := Concatenation( functionName,  recordstructure[1] );
    Assert(0, Size( recordstructure[1])>0);
    if not recordstructure[1][1]='@' then 
        globalName := Concatenation( functionName, "@", recordstructure[1] );
    fi;
    if IS_READ_ONLY_GLOBAL(globalName) then 
      MakeReadWriteGlobal( globalName );
    fi;
    if IsBoundGlobal(globalName) then 
        UnbindGlobal(  globalName ) ;    
    fi;
    BindGlobal( globalName ,  currentRec.(functionName) ) ;
    #MakeReadOnlyGlobal(globalName);
    SetName(currentRec.(functionName),fullName );
   
    MakeReadWriteGlobal( recordstructure[1] );
    UnbindGlobal( recordstructure[1] );
    MakeImmutable( headRec );
    BindGlobal( recordstructure[1], headRec);
   return;
end
);


 InstallGlobalFunction( InstallGlobalRecordMethod@FR,
 function( recordstructure, functionName, comments, params, functionRef)
     InstallGlobalRecordFunctionOrMethod@FR( recordstructure, functionName, comments, params, functionRef, InstallMethod, false);
 end
);


 InstallGlobalFunction( InstallGlobalRecordOtherMethod@FR,
 function( recordstructure, functionName, comments, params, functionRef)
     InstallGlobalRecordFunctionOrMethod@FR( recordstructure, functionName, comments, params, functionRef, InstallOtherMethod,false);
 end
);


InstallGlobalFunction( InstallGlobalRecordFunctionEx@FR,
function( recordstructure, functionName, functionRef, reinstall)
    local installFkt;
    installFkt := function( variable, comments, params, functionRef, reinstall )
        variable := functionRef;
    end;
    InstallGlobalRecordFunctionOrMethod@FR( recordstructure, functionName, "",[], functionRef, installFkt, reinstall);
end
);



InstallGlobalFunction( InstallGlobalRecordFunction@FR,
function( recordstructure, functionName, functionRef)
     Assert(0, not IsOperation(functionRef) );
     InstallGlobalRecordFunctionEx@FR( recordstructure, functionName,  functionRef,false);
end
);

InstallGlobalFunction( InstallGlobalRecordOperation@FR,
function( recordstructure, functionName, parameterTypes)
     InstallGlobalRecordFunctionEx@FR( recordstructure, functionName,  NewOperation(functionName, parameterTypes), false);
end
);

InstallGlobalFunction( ReInstallGlobalRecordFunction@FR,
function( recordstructure, functionName, functionRef)
    Assert(0, not IsOperation(functionRef) );
    InstallGlobalRecordFunctionEx@FR( recordstructure, functionName,  functionRef, true);
end
);


InstallGlobalFunction( ReInstallGlobalRecordOperation@FR,
function( recordstructure, functionName, parameterTypes)
    InstallGlobalRecordFunctionEx@FR( recordstructure, functionName,  NewOperation(functionName, parameterTypes),true);
end
);



 
  InstallGlobalRecordOperation@FR ( ["@FR@Utils"], "Degree", 
 [IsPolynomial] 
 );
 
#InstallOtherMethod( Degree , 
 InstallGlobalRecordMethod@FR ( ["@FR@Utils"], "Degree",
"get degree of a multivariate polynomial", [IsPolynomial], 
 function( polynomial )
    
    local   coeffData, monomData, degree,monomialDegrees, pos, i;

    coeffData := ExtRepPolynomialRatFun(polynomial);
    monomialDegrees := List([1..Size(coeffData)/2]);
    for pos in [1..Size(coeffData)/2] do
        degree := 0;
	    if not IsZero(coeffData[pos*2]) then 
	            monomData := coeffData[pos*2-1];
	            for i in [1..Size(monomData)/2] do
	                degree:=degree+monomData[i*2];
	            od;
	    fi;
        monomialDegrees[pos]:=degree;
    od;
   return Maximum(monomialDegrees);
end
);




###############################################################################################################


InstallMethod( FlattenList@FR, "remove the top level nesting ", [IsList], 
function(list)
    local result, entry;
    
    Assert(0, IsList(list));
    
    result := [];
    for entry in list do
        if IsList(entry) then
        Append(result,entry);
        else
            Append( result, [entry] );
        fi;
    od;
    list := result;
    return list;
end
);

InstallMethod( FirstElement@FR, "get first list element ", [IsList], 
function(list)
    if Size(list)=0 then 
        return fail;
    fi;
    return list[1];
end
);

InstallMethod( LastElement@FR, "get last list element ", [IsList], 
function(list)
    if Size(list)=0 then 
        return fail;
    fi;
    return list[Size(list)];
end
);


InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_FLATTEN_LIST", 
function()
	Assert(0, [] = FlattenList@FR( [] ));
	Assert(0, [1,2,1] = FlattenList@FR( [1,[2,1]] ));
	Assert(0, [1,2,[1]] = FlattenList@FR( [1,[2,[1]]] ));
	Assert(0, [[1],1] = FlattenList@FR( [[],[[1]],1] ));
end
);


#################################### GET/SET COEFFICIENTS ################################################################


InstallMethod( IsMonomial@FR, "", [IsObject], 
function (monomial)
	local  monomData;
	if not IsPolynomial (monomial) then 
		return false;
	fi;
	monomData := ExtRepPolynomialRatFun(monomial);
	if Size(monomData) <>2 or not IsOne(monomData[2]) then
		return false;
	fi;
	return true;
end
);


InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_IS_MONOMIAL", 
function()
    local rng, indet, x, y;
    rng := PolynomialRing( ZmodnZ(11)  ,["x","y"] );
    indet  := IndeterminatesOfPolynomialRing(rng);
    x := indet[1];
    y := indet[2];
    Assert(0, IsMonomial@FR(x) );
    Assert(0, IsMonomial@FR(x*y) );
    Assert(0, not IsMonomial@FR(2*x*y) );
    Assert(0, not IsMonomial@FR(x+y) );
    Assert(0, not IsMonomial@FR(3) );
    Assert(0, not IsMonomial@FR(rng) );
end);


InstallMethod( MonomialCoefficient@FR , 
"get coefficient for a given monomial of an polynomial ", [ IsPolynomial, IsPolynomial ],
 function( polynomial, monomial )
    
    local  monomData, coeffData, pos;
 
    if not IsMonomial@FR ( monomial ) then 
    	Error( "MonomialCoefficient: second parameter is not a monomial !" );
    fi;
    
    monomData := ExtRepPolynomialRatFun(monomial);

    coeffData := ExtRepPolynomialRatFun(polynomial);
    for pos in [1..Size(coeffData)/2] do
        if coeffData[pos*2-1]=monomData[1] then
            return coeffData[pos*2];
        fi;
    od;    
   
    return Zero( CoefficientsFamily(FamilyObj(polynomial)) ) ;
end
);


# get coefficients of specific monomials.
# todo:  implementation is not efficient
InstallOtherMethod( Coefficients@FR, 
" get coefficients of specified monomials", [IsPolynomial, IsList],
 function( polynomial, monomials )
     local monomial, monomialCoefficient, coefficients;
   
     # checking:
      for monomial in monomials do     
      	 if not IsMonomial@FR(monomial) then
      	 	Error( "getCoefficients: second parameter has to be a monomial list !" );
      	 fi;
      od;    
   
    coefficients := [];
    for monomial in monomials do
        monomialCoefficient := MonomialCoefficient@FR ( polynomial, monomial );
        Append( coefficients, [ monomialCoefficient ] );
    od;
    return coefficients;
end
);

InstallMethod( CoefficientsEx@FR, 
" get coefficients of specified monomials", [IsPolynomial, IsList],
 function( polynomial, monomials )
    return [ Coefficients@FR( polynomial, monomials), monomials];
 end
 );


InstallOtherMethod( CoefficientsEx@FR , 
"get nonzero coefficients and corresponding monomial list for a polynomial", [IsPolynomial], 
 function( polynomial )
    
    local  coeffList, coeffData, pos, monomialList, idCoeff;
    coeffList := [];
    monomialList := [];

    idCoeff := One( CoefficientsFamily( FamilyObj(polynomial) ) );
    coeffData := ExtRepPolynomialRatFun(polynomial);
    for pos in [1..Size(coeffData)/2] do
	if not IsZero(coeffData[pos*2]) then 
	        Append( coeffList, [ coeffData[pos*2] ]);
	        Append( monomialList, [ PolynomialByExtRep( FamilyObj( polynomial), [ coeffData[pos*2-1] , idCoeff ]   ) ]    );
	fi;
    od;
   return [ coeffList, monomialList ];
end
);




InstallOtherMethod( Coefficients@FR, 
" get coefficients of specified monomials", [IsPolynomial],
 function( polynomial )
    return CoefficientsEx@FR( polynomial)[1] ;
 end
 );


#InstallGlobalFunction( TEST_MONOMIAL_COEFFICIENT@FR ,
InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_MONOMIAL_COEFFICIENT", 
 function()
    local rng, indeterminates,x,y,polynomial;
    rng := PolynomialRing( ZmodnZ(11)  ,["x","y"] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
    y := indeterminates[2];
    polynomial := (x^4-4)^3*(4*y^2+2);
    Assert(0, Z(11)^4 = MonomialCoefficient@FR(polynomial, x^4*y^2));
    Assert(0, Zero(Z(11)) = MonomialCoefficient@FR(polynomial, x^42*y^2));
    Assert(0, Z(11)^2 = MonomialCoefficient@FR(polynomial, x^0*y^0));
end
);


#InstallGlobalFunction( TEST_COEFFICIENTS@FR ,
InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_COEFFICIENTS", 
function()
    local rng, indeterminates,x,y,polynomial;
     rng := PolynomialRing( ZmodnZ(11)  ,["x","y"] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
    y := indeterminates[2];
    polynomial := (x^4-4)^3*(4*y^2+2);
    Assert(0, [Z(11)^4, Zero(Z(11))] =Coefficients@FR(polynomial, [x^4*y^2, x^42*y^2]));    

end
);


#################################### POLYNOMIAL DIFFERENTIATION ################################################################

# Jacobian: compute ( d[fktlist_i] / d[indeterminants]_j )
InstallGlobalFunction( Jacobian@FR, 
function( fktlist, indeterminants )
    local cols, mat, row, col, fkt;
    if not  IsList(fktlist) then
	    Error("Jacobian: first parameter has to be a list of polynomials! \n");
    fi;
    if not  IsList(indeterminants) then
	    Error("Jacobian: second parameter has to be a list of indeterminates! \n");
    fi;
    
    for fkt in fktlist do
    	if not IsPolynomial(fkt) then
    		Error("Jacobian: first parameter has to be a list of polynomials! \n");
    	fi;
    od;
    
    mat:=   List( [1..Size(fktlist)], n->
                                        List( [1..Size(indeterminants)], l->0) 
                );
    for row in [1..Size(fktlist)] do 
        for col in  [1..Size( indeterminants)] do 
            mat[row][col] := Derivative( fktlist[row], indeterminants[col] ) ;
        od;
    od;
    return mat;
end
);


#InstallGlobalFunction(TEST_JACOBIAN@FR, 
InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_JACOBIAN", 
function() 
    local rng, ind,x,y,scalar,pol,jacobian;
    rng := PolynomialRing(Rationals,2);
    ind := IndeterminatesOfPolynomialRing(rng);
    x := ind[1];
    y := ind[2];
    scalar:=5/3;
    pol := scalar*x;
    
    jacobian := Jacobian@FR( [pol,y^2], ind);
    Assert(0, jacobian = [ [Derivative(pol,x), Derivative(pol,y)], [Derivative(y^2,x),Derivative(y^2,y)] ] );
end
);


#################################### COERCE POLYNOMIALS AND SCALARS #########################################################

InstallGlobalFunction( CoerceScalar@FR ,
function(scalar, dstRing)
    local intVal, coercedVal;

    coercedVal :=   scalar ;
    if Int(scalar)* One(scalar)=scalar then 
        coercedVal := Int(scalar) ;
         if Characteristic( scalar )>0 and Int(coercedVal)>Characteristic( scalar )/2 then
            coercedVal := coercedVal-Characteristic( scalar );
         fi;
        coercedVal := coercedVal* One(dstRing);
    fi;
    coercedVal := coercedVal * One(dstRing);
 
    if   IsRing(dstRing) then 
        Assert(0, coercedVal in dstRing);
    fi;
    if  IsFamily(dstRing) then
        Assert(0, FamilyObj(coercedVal) = dstRing );
    fi;
    return coercedVal;
end
);


#InstallGlobalFunction( TEST_COERCE_SCALAR@FR ,
InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_COERCE_SCALAR", 
function()
    local scalar, dstRing;
    scalar := 1/3;
    dstRing := Integers;

    # CoerceScalar@FR( scalar,dstRing ); # TODO: fails and should fail!, but how it can be used in a test ?

    dstRing := GF(11);
    CoerceScalar@FR( scalar,dstRing ); #ok.

    dstRing := ZmodnZ(11);
    Assert(0, One(dstRing)*scalar=  CoerceScalar@FR( scalar,dstRing )); #ok.

    scalar := 23;
    dstRing := Integers;

     Assert(0, One(dstRing)*scalar=  CoerceScalar@FR( scalar,dstRing )); #ok.

    dstRing := GF(11);
     Assert(0, One(dstRing)*scalar=  CoerceScalar@FR( scalar,dstRing )); #ok.

    dstRing := ZmodnZ(121);
   Assert(0, One(dstRing)*scalar=  CoerceScalar@FR( scalar,dstRing )); #ok.
    
end
);


# note: will not work for Galois fields and floats!
# has some problems for iterative dstRings...
InstallGlobalFunction( CoercePolynomial@FR ,
function( polynomial, dstRing )
   local pos, intVal,fam,  polRep, polRepCopy, coercedPol, coercedVal, scalar;

  fam := ElementsFamily(FamilyObj( dstRing ));

   polRep := ExtRepPolynomialRatFun( polynomial );
    polRepCopy := ShallowCopy(polRep);
    for pos in [1..Size(polRep)/2] do
        if IsPolynomial( polRepCopy[2*pos]) and IsPolynomialRing(dstRing) then
            polRepCopy[2*pos] := CoercePolynomial@FR(  polRep[2*pos],  CoefficientsRing(dstRing) )  ;    
        else
	         polRepCopy[2*pos] := CoerceScalar@FR(  polRep[2*pos],  CoefficientsRing(dstRing) )  ;      
        fi;
    od;        
    coercedPol :=  PolynomialByExtRep(fam, polRepCopy);
    return coercedPol;
end
);




#InstallGlobalFunction( TEST_COERCE_POLYNOMIAL@FR ,
InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_COERCE_POLYNOMIAL", 
function()
      local rng, ind, x, pol, dstRng, baseField, scalar, dstRing, coercedPol, dstInd, expectedResult;

    rng := PolynomialRing(Rationals,1);
    ind:=IndeterminatesOfPolynomialRing(rng);
    x := ind[1];
    scalar:=5/3;
    pol := scalar*x;

    baseField := ZmodnZ(11);
    dstRng := PolynomialRing(baseField ,1);
    coercedPol := CoercePolynomial@FR(pol, dstRng);
    dstInd := IndeterminatesOfPolynomialRing(dstRng);
    expectedResult := dstInd[1]*Z(11)^6;
    Assert(0, coercedPol=expectedResult);

    # CoerceScalar@FR( scalar,dstRng ); #ok.

    baseField := ZmodnZ(121);
    dstRng := PolynomialRing(baseField ,1);
    coercedPol := CoercePolynomial@FR(pol, dstRng);
    dstInd := IndeterminatesOfPolynomialRing(dstRng);
    expectedResult := dstInd[1]*ZmodnZObj(42,121);
    Assert(0, coercedPol=expectedResult);
     
    CoerceScalar@FR( scalar,dstRng ); #ok.

   # baseField := Integers;
   # dstRng := PolynomialRing(baseField ,1);
   # CoercePolynomial@FR(pol, dstRng);  # fails and probably should fail, but how to use in a test?
end
);


# coerce polynomial or scalar elements in vec to elements in dstRing.
# works only for prime fields...
InstallGlobalFunction( CoerceTensor@FR ,
function( tensor,  dstRing )
    local coercedTensor,polRepCopy, coercedPol, coordinate, pos, polRep, fam, intVal;
 
    if not IsList(tensor) then
        return  CoerceTensor@FR( [tensor], dstRing)[1];
    fi;
    coercedTensor := List( [1..Size(tensor)], n->0 ) ;
        for coordinate in  [1..Size( tensor)] do 
            if IsList( tensor[coordinate] ) then
                 coercedTensor[coordinate] := CoerceTensor@FR( tensor[coordinate], dstRing );
            else
                if IsPolynomial( tensor[coordinate] ) then
                    coercedTensor[coordinate] :=  CoercePolynomial@FR( tensor[coordinate], dstRing) ;
                else
                    coercedTensor[coordinate] :=  CoerceScalar@FR( tensor[coordinate] , dstRing) ;
                fi;
            fi;
        od;
    return coercedTensor;
end
);





#################################### EVALUATE POLYNOMIALS ################################################################



# EvalPolynomialTensor: substitute all indeterminates in the tensor by corresponding values.
# parameters:  ( tensor, indeterminates , values )
# precondition: tensor elements are polynomials over indeterminates and indeterminates belong to the same ring.
# postconditon:  indeterminates[i] is replaced by values[i];
# note: why was/ is the multiplication with One(values[1]) necessary?  
#        => assume we want to evaluate an polynomial 1/3*x^0 where x is ZmodnZObj( 1, 121 ). 
#       The result will be 1/3 while  1/3* ZmodnZObj( 1, 121 ) will give ZmodnZObj( 81, 121 ).
# end note
# todo: check if all values belongs to the same Ring.
#
InstallGlobalFunction( EvalPolynomialTensor@FR ,
 function( tensor, indeterminates , values )
    local pos,  evaluatedTensor;

    if not  Size(indeterminates) = Size(values) then 
    	Error("EvalPolynomialTensor: number of indeterminates and values must be the same");
    fi;
    
    for pos in [1..Size(values)] do 
    	if not One( values[1] )=One(values[pos]) then
    		Error("EvalPolynomialTensor: values must belong to the same ring ");
    	fi;
    od;
     

    if not IsList(tensor) then
        return  EvalPolynomialTensor@FR( [tensor], indeterminates , values )[1];
    fi;

     evaluatedTensor :=  List( [1..Size(tensor)], n->0);
    for pos in [1..Size(tensor)] do 
          if IsList( tensor[pos] ) then
              evaluatedTensor[pos]:= EvalPolynomialTensor@FR (tensor[pos], indeterminates, values);
          else
              evaluatedTensor[pos]:= One(values[1])* Value( tensor[pos], indeterminates, values );
               # evaluatedTensor[pos]:=Value( tensor[pos], indeterminates, values );
          fi;
    od;
  
    return evaluatedTensor;
end
);


# for some reason  EvalPolynomialTensor with implicit coertion (multiplication with One) was bad, and therefore I wrote EvalPolynomialTensorWeak
# but I forgot the caused problem
#
#InstallGlobalFunction( EvalPolynomialTensorWeak ,
# function( tensor, indeterminates , values )
#    local pos,  evaluatedTensor;
#
#    if not  Size(indeterminates) = Size(values) then 
#    	Error("EvalPolynomialTensorWeak: number of indeterminates and values must be the same");
#    fi;
#    
#    for pos in [1..Size(values)] do 
#    	if not One( values[1] )=values[pos] then
#    		Error("EvalPolynomialTensorWeak: values must belong to the same ring ");
#    	fi;
#    od;
#
#    if not IsList(tensor) then
#        return  EvalPolynomialTensorWeak( [tensor], indeterminates , values )[1];
#    fi;
#
#     evaluatedTensor :=  List( [1..Size(tensor)], n->0);
#    for pos in [1..Size(tensor)] do 
#          if IsList( tensor[pos] ) then
#              evaluatedTensor[pos] := EvalPolynomialTensorWeak (tensor[pos], indeterminates, values);
#          else
#               evaluatedTensor[pos] := Value( tensor[pos], indeterminates, values );
#          fi;
#    od;
#  
#    return evaluatedTensor;
#end
#);


#InstallGlobalFunction( TEST_EVAL_POLYNOMIAL_TENSOR@FR ,
InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_EVAL_POLYNOMIAL_TENSOR", 
 function()
    local rng, ind, x,y,mat,dstRng, evaluatedTensor, weakEvaluatedTensor;
    rng := PolynomialRing(Rationals,2);
    ind := IndeterminatesOfPolynomialRing(rng);
    x:=ind[1];
    y:=ind[2];
    mat:=[[1/3*x^0, x^0,x+y]];
    dstRng := ZmodnZ(121);
    
    evaluatedTensor := EvalPolynomialTensor@FR(mat, [x,y],[ZmodnZObj( 1, 121 ),ZmodnZObj( 2, 121 )]);

    #weakEvaluatedTensor := EvalPolynomialTensorWeak@FR(mat, [x,y],[ZmodnZObj( 1, 121 ),ZmodnZObj( 2, 121 )]);
    #CoerceTensor@FR( weakEvaluatedTensor, dstRng); # liefert_e_ nicht das was man erwartet 
 
    EvalPolynomialTensor@FR(mat, [x,y],[ZmodnZObj( 1, 121 ),ZmodnZObj( 2, 121 )]);
end
);


InstallGlobalFunction( SUBSTITUTE_POLYNOMIAL_COEFFICIENTS@FR,
 function( vec, ind , solution, dstFam  )
    local pos,  evaluatedVec, fam, polRep,coeffPos, polRepCopy,  coeffVal , coercedPol ;

    if not  Size(ind) = Size(solution) then 
    	Error("SubstitutePolynomialCoefficients: number of indeterminates and values must be the same");
    fi;
    
    for pos in [1..Size(solution)] do 
    	if not One( solution[1] )=One(solution[pos]) then
    		Error("SubstitutePolynomialCoefficients: because of impicit coercion solution elements expected belong to the same ring ");
    	fi;
    od;
    
    if not IsList(vec) then
        return  SUBSTITUTE_POLYNOMIAL_COEFFICIENTS@FR( [vec], ind , solution, dstFam )[1];
    fi;

   
    evaluatedVec :=  List( [1..Size(vec)], n->0);
    for pos in [1..Size(vec)] do 
          if IsList( vec[pos] ) then
              evaluatedVec[pos] := SUBSTITUTE_POLYNOMIAL_COEFFICIENTS@FR (vec[pos], ind, solution,dstFam);
          else
            polRep := ExtRepPolynomialRatFun(  vec[pos] );
            polRepCopy := ShallowCopy(polRep);
            for coeffPos in [1..Size(polRep)/2] do
                coeffVal := Value( polRep[2*coeffPos],ind,  solution );
                polRepCopy[2*coeffPos] := coeffVal*One( solution[1] );
            od;      
            coercedPol :=  PolynomialByExtRep(dstFam, polRepCopy);
            evaluatedVec[pos] := coercedPol;
          fi;
    od;
    return evaluatedVec;
end
);


#InstallGlobalFunction( TEST_SUBSTITUTE_POLYNOMIAL_COEFFICIENTS@FR,
InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_SUBSTITUTE_POLYNOMIAL_COEFFICIENTS", 
function()
    local rng, ind, a, b, iterRng,iterInd, x, y, pol, dstRng, dstFam, result,
    PREV_ITER_POLY_WARN;
    
    rng := PolynomialRing(Rationals,3);
    ind := IndeterminatesOfPolynomialRing(rng);
    a := ind[1];
    b := ind[2];
    PREV_ITER_POLY_WARN := ITER_POLY_WARN;
    ITER_POLY_WARN := false;
    iterRng := PolynomialRing(rng,2);
    ITER_POLY_WARN := PREV_ITER_POLY_WARN;
    
    iterInd := IndeterminatesOfPolynomialRing(iterRng);
    x := iterInd[1];
    y := iterInd[2];
    
    pol := a*b*x+b*y;
    dstRng := PolynomialRing( Rationals,2 );
    dstFam := FamilyObj( One(dstRng) );
    result := SUBSTITUTE_POLYNOMIAL_COEFFICIENTS@FR( pol , ind , [2,1,0], dstFam  );
    Assert(0,  CoercePolynomial@FR(result,iterRng) = CoercePolynomial@FR(2*x+y, iterRng ) );
    
    #Assert(0,  result = CoercePolynomial@FR(2*x+y, dstRng ) );  
    
end
);

#################################### POLYNOMIAL PROPERTIES ################################################################

InstallMethod( CountPolynomialVariables@FR ,
"count variables appeared in a polynomial. Zero monomials are ignored. " , [IsPolynomial],
function(polynomial)
    local coeffData, variableIdx, pos, variableIndices, indIdx;

    if not IsPolynomial(polynomial) then
        Error("parameter is not a polynomial!");
    fi;
    
    variableIndices := [];
    coeffData := ExtRepPolynomialRatFun(polynomial);
    for pos in [1..Size(coeffData)/2] do
        for indIdx in [1..Size( coeffData[2*pos-1])/2] do
              Append( variableIndices, [ coeffData[2*pos-1][2*indIdx-1] ] );
        od;
    od;
  
    return Size( Set(variableIndices) );
end
);


#InstallGlobalFunction( TEST_COUNT_POLYNOMIAL_VARIABLES@FR ,
InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_COUNT_POLYNOMIAL_VARIABLES", 
function()
    local rng, indeterminates,x,y;
    rng := PolynomialRing( ZmodnZ(11)  ,["x","y"] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
    y := indeterminates[2];

    Assert(0, CountPolynomialVariables@FR(y)=1);
    Assert(0, CountPolynomialVariables@FR(x*y)=2);
    Assert(0, CountPolynomialVariables@FR(x+y)=2);
end
);


InstallGlobalRecordFunction@FR ( ["@FR@Utils",], "IsMonic", 
function(pol)
    Assert(0, IsUnivariatePolynomial(pol));
    return IsOne(LeadingCoefficient(pol));
end
);


InstallGlobalRecordFunction@FR ( ["@FR@Utils",], "TEST_IS_MONIC", 
function()
    local rng, indeterminates, x;
    rng := PolynomialRing( ZmodnZ(11)  ,["x"] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
    Assert(0, IsMonic@FR@Utils(x));  
    Assert(0, not IsMonic@FR@Utils(2*x));
end
);


InstallGlobalRecordFunction@FR ( ["@FR@Utils",], "IsIndeterminate", 
function(variable)
    return IsUnivariatePolynomial(variable) and Degree(variable)=1 and IsOne(LeadingCoefficient(variable)) and Size(ExtRepPolynomialRatFun(variable))=2;
end
);

InstallGlobalRecordFunction@FR ( ["@FR@Utils",], "TEST_IS_INDETERMINATE", 
function()
    local rng, indeterminates, x,y;
    rng := PolynomialRing( ZmodnZ(11)  ,["x","y"] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
    y := indeterminates[2];
    Assert(0, IsIndeterminate@FR@Utils(x));
    Assert(0, IsIndeterminate@FR@Utils(y));
    Assert(0, not IsIndeterminate@FR@Utils(y+x));
    Assert(0, not IsIndeterminate@FR@Utils(1+x));
    Assert(0, not IsIndeterminate@FR@Utils(2*x));
end
);


 InstallGlobalRecordOperation@FR ( ["@FR@Utils"], "HomogenizedPolynomial", 
  [IsPolynomial, IsObject,IsObject] 
 );
 

InstallGlobalRecordMethod@FR ( ["@FR@Utils"], "HomogenizedPolynomial","homogenize polynomial", [IsPolynomial, IsObject, IsObject],
function(pol, homogenVariable, degree)

    local coeffData, coeffs, monomials,newPol;
    if IsZero(pol) then
        return pol;
    fi;
    
    Assert(0,IsIndeterminate@FR@Utils(homogenVariable) );
    Assert(0, degree>= Degree@FR@Utils(pol) );
    coeffData := CoefficientsEx@FR(pol);
    coeffs := coeffData[1];
    monomials := coeffData[2];
    monomials:= List(monomials, monomial-> monomial*homogenVariable^(degree- Degree@FR@Utils(monomial)) );
    newPol:= coeffs*monomials;
    return newPol;
end
);


InstallGlobalRecordOtherMethod@FR ( ["@FR@Utils"], "HomogenizedPolynomial","homogenize polynomial", [IsPolynomial, IsObject],
function(pol, homogenVariable)

     return @FR@Utils.HomogenizedPolynomial(pol, homogenVariable,  Degree@FR@Utils(pol) );
end
);


 ReInstallGlobalRecordOperation@FR ( ["@FR@Utils"], "IndeterminatesOfPolynomial", 
  [IsPolynomial] 
 );
 
# todo: had a bug. Missing test...
InstallGlobalRecordMethod@FR ( ["@FR@Utils"], "IndeterminatesOfPolynomial","get polynomial indeterminates", [IsPolynomial],
function(polynomial)
    local result, coeffData, monomialData, monomialIdx, indeterminateIdPos, variableIds, variableList, one;
    
    result := [];
    coeffData := ExtRepPolynomialRatFun(polynomial);
    
    variableIds := [];
    if Size(coeffData)<1 then
        return [];
    fi;
    for monomialIdx in [1..Size(coeffData)/2] do
    
        monomialData := coeffData[monomialIdx*2-1];
        for indeterminateIdPos in [1..Size(monomialData)/2] do
            Append(variableIds, [ monomialData[indeterminateIdPos*2-1] ]);
        od;
    od;
    variableIds := Set(variableIds);
    one := One( CoefficientsFamily(FamilyObj(polynomial)));
   
    variableList := List(variableIds, id->  PolynomialByExtRep( FamilyObj( polynomial), [ [ id, 1], one ] ) );  
    
    return variableList;
end
);



 InstallGlobalRecordOperation@FR ( ["@FR@Utils"], "IsHomogenized", 
  [IsPolynomial] 
 );
 

InstallGlobalRecordMethod@FR ( ["@FR@Utils"], "IsHomogenized","check if  polynomial is homogenized", [IsPolynomial],
function(pol)

    local coeffData, coeffs, monomials,monomialDegrees;
    if IsZero(pol) then
        return true;
    fi;
    
    coeffData := CoefficientsEx@FR(pol);
    coeffs := coeffData[1];
    monomials := coeffData[2];
    monomialDegrees := List( monomials, monomial-> Degree@FR@Utils(monomial) );
    monomialDegrees := Set(monomialDegrees);
    return Size(monomialDegrees)=1;
end
);

 

 InstallGlobalRecordOperation@FR ( ["@FR@Utils"], "DehomogenizedPolynomial", 
  [IsPolynomial, IsObject] 
 );
 
 
InstallGlobalRecordMethod@FR ( ["@FR@Utils"], "DehomogenizedPolynomial","dehomogenize polynomial", [IsPolynomial,IsPolynomial],
function(pol, homogenVariable)
    if not IsHomogenized@FR@Utils(pol)  then
        Error("polynomial is not homogenized");
    fi;
    
     if not IsIndeterminate@FR@Utils(homogenVariable) then
        Error ("parameter homogenVariable is not a variable");
     fi;
     if not homogenVariable in IndeterminatesOfPolynomial@FR@Utils(pol) then
        Error ("parameter homogenVariable is not an inteterminate of the polynomial");
     fi;
     return Value( pol, [ homogenVariable ], [ One(homogenVariable)] );
end
);


InstallGlobalRecordOtherMethod@FR ( ["@FR@Utils"], "DehomogenizedPolynomial","dehomogenize polynomial", [IsPolynomial],
function(pol)

    local coeffData, coeffs, monomials,newPol,indeterminates;

    if not IsHomogenized@FR@Utils(pol)  then
        Error("polynomial is not homogenized");
    fi;
    if IsZero(pol) then
        return pol;
    fi;
    indeterminates := IndeterminatesOfPolynomial@FR@Utils(pol);
    if Size(indeterminates)<=1 then
        return pol;
    fi;
    return DehomogenizedPolynomial@FR@Utils( pol, indeterminates[ Size(indeterminates)] );
end
);



InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_DEHOMOGENIZE_POLYNOMIAL", 
function()
    local rng, ind, x,y, pol, hpol, dhpol;

	 rng := PolynomialRing(Rationals,["x","y"]);
	 #rng := PolynomialRing(GF(121),["x","y"]);
	 ind := IndeterminatesOfPolynomialRing(rng);
	 x:=ind[1];
	 y:=ind[2];
	 
	 pol := 2*(2*x^2-3)^2*(x-4);
	 hpol:=HomogenizedPolynomial@FR@Utils(pol,y);
	 dhpol:=DehomogenizedPolynomial@FR@Utils(hpol,y);
	 Assert(0, dhpol=pol);
	 dhpol:=DehomogenizedPolynomial@FR@Utils(hpol);
	 Assert(0, dhpol=pol);
end
);





InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_HOMOGENIZE_POLYNOMIAL", 
function()
    local rng, ind, x,y, pol, hpol, coeffData, monom, monomials;

	 rng := PolynomialRing(Rationals,["x","y"]);
	 #rng := PolynomialRing(GF(121),["x","y"]);
	 ind := IndeterminatesOfPolynomialRing(rng);
	 x:=ind[1];
	 y:=ind[2];
	 
	 pol := 2*(2*x^2-3)^2*(x-4);
	 Assert(0, not IsHomogenized@FR@Utils(pol) );
	 hpol:= @FR@Utils.HomogenizedPolynomial(pol,y,6);
	 Assert(0,   IsHomogenized@FR@Utils(hpol) );
	 coeffData := CoefficientsEx@FR(hpol);
	 monomials := coeffData[2];
	 for monom in monomials do
	    Assert(0, Degree@FR@Utils(monom)=6 );
	 od;	  
	 
	 hpol:= @FR@Utils.HomogenizedPolynomial(pol,y);
     Assert(0,   IsHomogenized@FR@Utils(hpol) );
	 coeffData := CoefficientsEx@FR(hpol);
	 monomials := coeffData[2];
	 for monom in monomials do
	    Assert(0, Degree@FR@Utils(monom)=5 );
	 od;	 
end
);


# check if polynomiyl is constant: use  IsConstantRationalFunction

#################################### POLYNOMIAL FACTORS AND PRODUCTS ################################################################

InstallMethod( IsPower@FR, "check if data structure is a power data structure", [IsList],
function(power)
	if Size(power)<>2  then
		return false;
	fi;
	if fail=ApplicableMethod(\^,[power[1],power[2]]) then
		return false;
	fi;
	return true;
end
);


InstallGlobalFunction( CreatePower@FR,
function(base,exponent)
  if IsPower@FR( [base, exponent ] ) then
  	return [ base, exponent ];
  fi;
  Assert(0, false);
end
 );

# return factors of a polynomial with the property that for each pair of the factors their Gcd is always at most a constant.
# also the unique factors do not contain scalars (only factors of degree>0!)
 InstallMethod( DistinctMonicFactors@FR , 
 " # return distinct monic factors (no constants) of an univariate polynomial. ", [ IsPolynomial ],
 function(polynomial)
 	local factors, factor1, factor2 ;
 	if not IsUnivariatePolynomial(polynomial) and not (IsHomogenized@FR@Utils(polynomial) and IndeterminateNumber@FR(polynomial)=2 ) then
 		Error("DistinctMonicFactors@FR: parameter is not an univatiate or homogenized polynomial");
 	fi;
 	factors :=  ShallowCopy (Factors( polynomial) );
     	factors[1] := StandardAssociate( factors[1] );
  
  	if (  Degree@FR@Utils(factors[1]) ) =0 then 
  		Remove(factors,1);
  	fi;
	factors:= Set ( factors );
	
	return factors;
end 
);


#InstallGlobalFunction( TEST_DISTINCT_MONIC_FACTORS@FR,
InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_DISTINCT_MONIC_FACTORS", 
function()
    local rng, indeterminates, x, pol, result;
    rng := PolynomialRing( ZmodnZ(11)  ,["x" ] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
    pol := 4*(x-3)^10*(3*x-2)^3;
    result := DistinctMonicFactors@FR(pol);
    Assert(0, result = [(x-3),(x-8) ]);
    pol := 4*x^0;
    result := DistinctMonicFactors@FR(pol);
    Assert(0, Size(result) = 0 );
end
);


# computes the value of a product (a product is a list of powers. A power is a pair [ base, exponent ] ).
InstallGlobalFunction( PRODUCT_VALUE@FR ,
function( product ) 
    local value, power;
   
    value := 1 ;
    for power in product do
        value := value* power[1]^power[2];
    od;
    return value;
end );


 
InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_PRODUCT_VALUE", 
function()

    local rng, indeterminates, x, product;
    
    rng := PolynomialRing( ZmodnZ(11)  ,["x" ] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
    
    product := [[2,3]];
    Assert(0, 2^3= PRODUCT_VALUE@FR( product));
    
    product := [ [x-3,3] ];
    Assert(0, (x-3)^3 = PRODUCT_VALUE@FR( product));
    
    product := [ [x-3,3], [x,2] ] ;
    Assert(0, (x-3)^3*x^2 = PRODUCT_VALUE@FR( product));
    
     product := [  ] ;
    Assert(0, 1 = PRODUCT_VALUE@FR( product));
end
);


InstallMethod( UNIQUE_PRODUCT_HOMOGENIZED@FR,
"factors an univariate polynomial over rationals or finite fields into  power factors ", [ IsPolynomial ],
function( polynomial )
    local localPolynomial, degree, factors, unit , uniqueFactors, uniqueProduct, uniqueProductPart, n, l,ind;
    
    if  not IsHomogenized@FR@Utils(polynomial) or IndeterminateNumber@FR(polynomial)<>2 then
 		Error("UNIQUE_PRODUCT_HOMOGENIZED: first parameter is not a  homogenized polynomial");
    fi;
    degree:=Degree@FR@Utils(polynomial);
    ind := IndeterminatesOfPolynomial@FR@Utils(polynomial);
    Assert(0, Size(ind)=2 );
    localPolynomial := DehomogenizedPolynomial@FR@Utils(polynomial,ind[2]);
    
    factors := ShallowCopy(Factors(localPolynomial));
    unit    := factors[1]/StandardAssociate( factors[1] );
    factors[1] := StandardAssociate( factors[1] );

    uniqueFactors := Set(factors);
        
    uniqueProduct := [  ];
    if not IsOne(unit) then
        Append( uniqueProduct, [[unit,1] ] );
    fi;
    uniqueProductPart := List( [1..Size(uniqueFactors)], 
                               n->[ HomogenizedPolynomial@FR@Utils(uniqueFactors[n],ind[2]), Number( factors , function(l) return l=uniqueFactors[n]; end ) ] );
    if Degree@FR@Utils(localPolynomial)< Degree@FR@Utils(polynomial) then
        Append(uniqueProductPart, [ [ ind[2],  Degree@FR@Utils(polynomial)-Degree@FR@Utils(localPolynomial)] ] );
    fi;
                                   
    Append( uniqueProduct,  uniqueProductPart  );

    return uniqueProduct;
    
end
);



# UNIQUE_PRODUCT: factors an univariate polynomial over rationals or finite fields into  power factors 
# ( a power is a pair of [base,exponent] ) with distinct bases
# see also 'Factors'
InstallMethod( UNIQUE_PRODUCT@FR,
"factors an univariate polynomial over rationals or finite fields into  power factors ", [ IsPolynomial ],
function( polynomial )
    local factors, unit , uniqueFactors, uniqueProduct, uniqueProductPart, n, l;
    
    if not IsUnivariatePolynomial(polynomial) and 
       not (IsHomogenized@FR@Utils(polynomial) and IndeterminateNumber@FR(polynomial)=2  ) then
 		Error("UNIQUE_PRODUCT: first parameter is not a univariate or homogenized polynomial");
    fi;
    if IsHomogenized@FR@Utils(polynomial) and IndeterminateNumber@FR(polynomial)=2 then
        return UNIQUE_PRODUCT_HOMOGENIZED@FR(polynomial);
    fi;
    factors := ShallowCopy(Factors(polynomial));
    unit    := factors[1]/StandardAssociate( factors[1] );
    factors[1] := StandardAssociate( factors[1] );

    uniqueFactors := Set(factors);
        
    uniqueProduct := [  ];
    if not IsOne(unit) then
        Append( uniqueProduct, [[unit,1] ] );
    fi;
    uniqueProductPart := List( [1..Size(uniqueFactors)], 
                               n->[ uniqueFactors[n], Number( factors , function(l) return l=uniqueFactors[n]; end ) ] );
    Append( uniqueProduct,  uniqueProductPart  );
    return uniqueProduct;
    
end
);
    

# just an alternative implamentation for UNIQUE_PRODUCT@FR
InstallMethod( UNIQUE_PRODUCT_1@FR, 
"factors an univariate polynomial over rationals or finite fields into  power factors ", [ IsPolynomial ],
function( polynomial )
   local uniqueProduct, factors, factor, multiplicity, tmp, scalarFactor , value ;
   if not IsUnivariatePolynomial(polynomial) and not IsHomogenized@FR@Utils(polynomial) then
 		Error("UNIQUE_PRODUCT_1: first parameter is not a univariate or homogenized polynomial");
   fi;
 	
    uniqueProduct := [];
    factors := DistinctMonicFactors@FR( polynomial) ;
    Degree@FR@Utils(polynomial);
    for factor in factors do
        tmp:=polynomial;
        multiplicity := 0;
        tmp:=tmp/factor;
        while  Degree@FR@Utils( DenominatorOfRationalFunction(tmp) )<=0 do
            tmp := tmp/factor;
            multiplicity:=multiplicity+1;
        od;
        Append( uniqueProduct, [ [factor, multiplicity] ] );
    od;

    value := PRODUCT_VALUE@FR(uniqueProduct);
    scalarFactor := polynomial/value;
    if not IsOne(scalarFactor) then
        Append( uniqueProduct, [ [scalarFactor, 1] ] );
    fi;

    return uniqueProduct;
end 
);



InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_UNIQUE_PRODUCT", 
function()
    local rng, indeterminates, x, y, expectedProduct, pol, result;
    rng := PolynomialRing( ZmodnZ(11)  ,["x","y" ] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
    y := indeterminates[2];
        
    pol := (x-3)^3 ;
    
    result := UNIQUE_PRODUCT@FR(pol);
    
    Assert(0, result = [ [ x-3, 3 ] ]);
    
    pol := 3*(x-3)^3 ;    
    result := UNIQUE_PRODUCT@FR(pol);
    Assert(0, result = [  [ x-3, 3 ],[ One(rng)*3,1]  ] );
    
    pol := (x-3)^3*x^2;
    result := UNIQUE_PRODUCT@FR(pol);
    expectedProduct := [  [x,2], [x-3,3] ] ;
    Assert(0, expectedProduct = result );
    
    pol := (x-3)^3*x^2;
    pol:=HomogenizedPolynomial@FR@Utils(pol,y,6);
    result := UNIQUE_PRODUCT@FR(pol);
    expectedProduct := [  [x,2], [x-3,3],[y,1] ] ;

     pol := x^0;
     result := UNIQUE_PRODUCT@FR(pol);
     expectedProduct := [  ] ;
        Assert(0, expectedProduct = result );

     pol := 5*x^0;
     result := UNIQUE_PRODUCT@FR(pol);
     expectedProduct := [ [One(rng)*5,1 ] ] ;
    Assert(0, expectedProduct = result );
    
end
);


# removes constant factors from a list of polynomial powers.
# (a polynomial constant factor is a power data [ base, exponent ] where Degree(base)=0 ) ;
  InstallGlobalFunction( REMOVE_CONSTANT_FACTORS@FR ,
function( powers )
    local factor, result;
    result:= [];
    for factor in  powers do
        Assert( 0, Size(factor)=2 );
        Assert( 0, factor[2] in Integers );
        if IsPolynomial( factor[1] ) and  Degree@FR@Utils(factor[1])>0 then
            Append(result, [ factor ]);
        fi;
    od;
    return result;
end );



# sort a list of powers ( a power is a pair of [base,exponent] ) by exponent 
  InstallGlobalFunction(  SORT_POWERS_BY_EXPONENT@FR , 
  function( factors )
    local  result, factor,    tmpFactors, currentExponent, currentFactorList;
    result := []; 
    tmpFactors := [];
    for factor in factors do
        Append(tmpFactors, [ [ factor[2],factor[1] ] ]) ;
    od;
    Sort(tmpFactors);
     
    currentExponent:=Null@FR;
    currentFactorList:=[];
    while Size(tmpFactors)>0 do
        if currentExponent=Null@FR or tmpFactors[1][1]>currentExponent then
            if Size(currentFactorList)>0 then
                Append(result, [ currentFactorList ]);
            fi;
            currentExponent:=tmpFactors[1][1];
            currentFactorList := [];
        fi;
        Append(currentFactorList, [ [ tmpFactors[1][2],tmpFactors[1][1] ] ]);
        Remove(tmpFactors,1);
    od;
    if Size(currentFactorList)>0 then
         Append(result, [ currentFactorList ]);
    fi;
    return result;
end 
);


InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_SORT_POWERS_BY_EXPONENT", 
function()
    local factors,
    sortedFactors, expectedResult;
    factors := [ [ 3,2 ], [ 3,1 ], [ 4,2 ], [ 3,3 ] ];
    sortedFactors := SORT_POWERS_BY_EXPONENT@FR( factors );
    
    expectedResult := [ [ [ 3,1 ]], [ [3,2], [ 4,2 ] ], [ [ 3,3 ] ] ];
   
    Assert(0,  expectedResult=sortedFactors );
    
    factors := [   ];
    sortedFactors := SORT_POWERS_BY_EXPONENT@FR( factors );
    
    expectedResult := [  ];
    Assert(0,  expectedResult=sortedFactors );
end
);


InstallGlobalRecordFunction@FR (["@FR@Utils","Internal"], "RemoveLineByLeadingString",
function( lines, leadingString, separators, last)
   local localRow;
   Assert(0, IsList(lines));
   Assert(0, IsList(separators));   
   Assert(0, last=true or last=false);
   Assert(0, IsString(leadingString) );
   
   if Size(lines)=0 then
        return lines;
   fi;
   
   localRow := SplitString(lines[1],separators);
   if last then
    localRow := SplitString(lines[Size(lines)],separators);
   fi;
   
   while not fail=Position(localRow,"") do Remove (localRow, Position(localRow,"")); od;
   
  
   if leadingString in localRow then
  
      Assert(0, localRow[1]=leadingString);
      
      if last then 
          lines := List([1..Size(lines)-1], i->lines[i]);
       else       
          lines := List([2..Size(lines)], i->lines[i]);
      fi;
   fi;
   return lines;   
end
);

# wenn local in der zweiten zeile, dann ... entferne die erste UND zweite Zeile, sonst nur die erste

# UNSTABLE! do not use extensively! Would either require a system function ' Function.body.toString()' 
#  ( function body without variable declarations ) or a begin and end body marker variable. 

# e.g. 
#
# testfkt:=
# function()
#   local var1, BEGIN_MARKER, END_MARKER;
#   BEGIN_MARKER;
#   var1:=3;
#   END_MARKER;
# end;

# String(@FR@Utils.Tests.TEST_SORT_POWERS_BY_EXPONENT);
InstallGlobalRecordFunction@FR (["@FR@Utils","Internal"], "CreateTestString",
function( testRecordVariableString, prefix)
    local testRecord, prefixString, str,strs, fullStr, name,strNew, pos, localRow, separators, line,strsCopy;

    testRecord := EvalString(testRecordVariableString);
    fullStr:="";
    str := "\n";
    
    separators := [' ',',',';' ];
            
   for name in RecNames(testRecord) do
    
       
        fullStr := Concatenation(fullStr, "#\n#\n");
        
        fullStr := Concatenation(fullStr, "# ",testRecordVariableString, ".", name, " : \n" );
        Assert(0, IsFunction(testRecord.(name)) );
     
        str := StringPrint(testRecord.(name));
       
        strs := SplitString(str,['\n']);
        
        strs := @FR@Utils.Internal.RemoveLineByLeadingString(strs, "function", separators, false);
      
        localRow := SplitString(strs[1],separators);
        while not fail=Position(localRow,"") do Remove (localRow, Position(localRow,"")); od;
        
               
        if (localRow[1]="local") then
         
            strs := JoinStringsWithSeparator(strs,"\n");
            strs := List([Position(strs,';')+1..Size(strs)], n->  strs[n] );
            strs := SplitString(strs,['\n']);
        fi;
        
        #strs:= @FR@Utils.Internal.RemoveLineByLeadingString(strs, "local", separators, false);
        strs:= @FR@Utils.Internal.RemoveLineByLeadingString(strs, "end", separators, true);
       
        
        strs:= @FR@Utils.Internal.RemoveLineByLeadingString(strs, "return", separators, true);
                
        

        str:="";
        for line in strs do 
            if Size(line)>0 and line[Size(line)]=';' then 
                line:=Concatenation(line,";");
            fi;          
            if prefix then
                line:=Concatenation("gap> ",line);
            fi;
            line:=Concatenation(line,"\n");
            str:= Concatenation(str,line);
        od;
     
        fullStr:= Concatenation(fullStr,"#\n" ,str);
    od;    
    return fullStr;
end
);



 InstallGlobalRecordOperation@FR ( ["@FR@Utils"], "LinearFactors", 
  [IsPolynomial,IsObject] 
 );
 
 # get distinct linear factors with given multiplicity
  InstallGlobalRecordMethod@FR ( ["@FR@Utils"], "LinearFactors", "get linear factors" , 
  [IsPolynomial, IsObject] ,
  function(polynomial, multiplicity)
    
        local power,   powerList, result;
        if not multiplicity=Null@FR and  not IsPosInt(multiplicity) then
            Error("expected second parameter to be multiplicity or Null@FR");
        fi;
        
        result := [];
        
         powerList := UNIQUE_PRODUCT@FR(polynomial);
         for power in powerList do
             if Degree@FR@Utils(power[1])=1 and 
                (multiplicity=Null@FR or 
                power[2] = multiplicity) then
                     Append( result, [power[1]] );
             fi;
         od;
         return result;
  end
 );
 
 
  # get distinct linear factors with arbitrary multiplicity
  InstallGlobalRecordOtherMethod@FR ( ["@FR@Utils"], "LinearFactors", "get linear factors" , 
  [IsPolynomial] ,
  function(polynomial)
       return @FR@Utils.LinearFactors(polynomial, Null@FR);
  end
 );
 

InstallGlobalRecordFunction@FR ( ["@FR@Utils","Tests"], "TEST_LINEAR_FACTORS", 
function()
    local rng, indeterminates, x, expectedProduct, pol, factors;
    rng := PolynomialRing( ZmodnZ(11)  ,["x" ] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
        
    pol := (x-3)^3*x*(x^2-2) ;
    
    factors := @FR@Utils.LinearFactors(pol,3);
    Assert(0, Size(factors)=1);
    
   factors := @FR@Utils.LinearFactors(pol);
    Assert(0, Size(factors)=2); 
end
);



 
InstallGlobalRecordFunction@FR (["@FR@Utils"], "CreateTestString",
function(prefix)
    return @FR@Utils.Internal.CreateTestString(" @FR@Utils.Tests", prefix);
end
);



#E hurwitzUtils.gi . . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
