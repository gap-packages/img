# This interface is explicitly not thread-safe by design!!
# please use the object-oriented one.



rationalMapFinder@.SEARCH_HURWITZ_MAP_OVER_FINITE_FIELD   := NewOperation("SEARCH_HURWITZ_MAP_OVER_FINITE_FIELD",[ IsList, IsBool, IsObject, IsPosInt] );


#############################################################################
##
#F SEARCH_HURWITZ_MAP_OVER_FINITE_FIELD( partitionList, strictNormalization, branchValueApproxList, fieldChar )
##
## <#GAPDoc Label="SEARCH_HURWITZ_MAP_OVER_FINITE_FIELD">
## <ManSection>
##   <Func Name="SEARCH_HURWITZ_MAP_OVER_FINITE_FIELD" Arg=" partitionList, strictNormalization, branchValueApproxList, fieldChar "/>
##   <Returns>A list of polynomialTuples over the finite field ZZ/fieldchar which satisfies the request</Returns>
##   <Description>
##     This function searches for polynomial sets over finite field ZZ/<C>fieldChar</C>satisfying multiplicity structure given by <C>partitionList</C>,
##	The number of polynomials of each set is equal to Length(<C>partitionList</C>)
## 	and the polynomials will satisfy the equations
##	polSet[1]-scaling[i]*polSet[2]=polSet[2+i]; i in 1..Length(partitionList)-2
##	
##
##     <P/>  The algorithm normalizes three factors in each polynomial set to [inf,0 1] respectively 
##	     If the parameter <C> strictNormalization</C> is true, 
##		then the algorithm tries to normalize a factor with multiplicity=partitionList[1][1] of first polynomial  in polynomial set to infinity ,
##		a factor with multiplicity=partitionList[2][1] of second polynomial  to zero  
##		and  factor with multiplicity=partitionList[3][1] of third polynomial to one
##		This is not always possible and thus less solutions will be found.
##
##     <P/>  To find corresponding compex approximations of the 
##
##
##     <P/> The following example ...
##
## <Example><![CDATA[
## gap> LoadPackage("fr");
## gap>	ReadPackage("fr","hurwitz/rationalMapFinder.gap");
## gap> finiteFieldSolutionList := rationalMapFinder@FR.SEARCH_HURWITZ_MAP_OVER_FINITE_FIELD( [ [1,2], [2,1], [2,1], [2,1] ], false, [ [0/1, -1/2] ] , 13 );
## gap> finiteFieldSolutionList;
##[ rec( polynomialRing := GF(13)[t,s], polynomialTuple := [ t*s^2+Z(13)^9*s^3, t^3+Z(13)*t^2*s, t^3+Z(13)*t^2*s+Z(13)^5*t*s^2+Z(13)^2*s^3, t^3+Z(13)*t^2*s+Z(13)^7*t*s^2+Z(13^4*s^3 ], 
##      polynomialTupleZZLift := [ t*s^2+5*s^3, t^3+2*t^2*s, t^3+2*t^2*s+6*t*s^2+4*s^3, t^3+2*t^2*s-2*t*s^2+3*s^3 ], rfsProblem := rec( scalingRelationList := [ [ 0, -1/2 ] ], ##shapeList := [ [ 2, 1 ], [ 2, 1 ], [ 2, 1 ], [ 2, 1 ] ] ) ), 
##  rec( polynomialRing := GF(13)[t,s], polynomialTuple := [ t*s^2+Z(13)^4*s^3, t^3+Z(13)^2*t^2*s, t^3+Z(13)^2*t^2*s+Z(13)*t*s^2+Z(13)^5*s^3, t^3+Z(13)^2*t^2*s+Z(13)^9*t*s^2+Z(13*s^3 ], 
##      polynomialTupleZZLift := [ t*s^2+3*s^3, t^3+4*t^2*s, t^3+4*t^2*s+2*t*s^2+6*s^3, t^3+4*t^2*s+5*t*s^2+2*s^3 ], rfsProblem := rec( scalingRelationList := [ [ 0, -1/2 ] ], shapeList := [ [ 2, 1 ], [ 2, 1 ], [ 2, 1 ], [ 2, 1 ] ] ) ), 
##  rec( polynomialRing := GF(13)[t,s], polynomialTuple := [ t*s^2+Z(13)^4*s^3, t^3+Z(13)^2*t^2*s, t^3+Z(13)^2*t^2*s+Z(13)*t*s^2+Z(13)^5*s^3, t^3+Z(13)^2*t^2*s+Z(13)^9*t*s^2+Z(13)*s^3 ], 
 ##     polynomialTupleZZLift := [ t*s^2+3*s^3, t^3+4*t^2*s, t^3+4*t^2*s+2*t*s^2+6*s^3, t^3+4*t^2*s+5*t*s^2+2*s^3 ], rfsProblem := rec( scalingRelationList := [ [ 0, -1/2 ] ], shapeList := [ [ 2, 1 ], [ 2, 1 ], [ 2, 1 ], [ 2, 1 ] ] ) ), 
##  rec( polynomialRing := GF(13)[t,s], polynomialTuple := [ t*s^2+Z(13)^4*s^3, t^3+Z(13)^2*t^2*s, t^3+Z(13)^2*t^2*s+Z(13)*t*s^2+Z(13)^5*s^3, t^3+Z(13)^2*t^2*s+Z(13)^9*t*s^2+Z(13)*s^3 ], 
##      polynomialTupleZZLift := [ t*s^2+3*s^3, t^3+4*t^2*s, t^3+4*t^2*s+2*t*s^2+6*s^3, t^3+4*t^2*s+5*t*s^2+2*s^3 ], rfsProblem := rec( scalingRelationList := [ [ 0, -1/2 ] ], shapeList := [ [ 2, 1 ], [ 2, 1 ], [ 2, 1 ], [ 2, 1 ] ] ) ) ]
## ]]></Example>
##   </Description>
## </ManSection>
## <#/GAPDoc>
##

#if not IsBoundGlobal ("HurwitzMapFinder@FR") then
#HurwitzMapFinder@FR:=fail;
#fi;
 

InstallMethod(
rationalMapFinder@.SEARCH_HURWITZ_MAP_OVER_FINITE_FIELD,
	"search for Hurwitz map candidates over a finite field. The solutions has to be lifted and approximated, see APPROX_HURWITZ_MAP_CANDIDATES",
	[ IsList, IsBool, IsObject, IsPosInt ],
	function( partitionList, strictNormalization, branchValueApproxList, fieldChar )
		if (HurwitzMapFinder@FR=fail) then
			Unbind(HurwitzMapFinder@FR);
		fi;
		if not IsBoundGlobal ("HurwitzMapFinder@FR") then
			BindGlobal("HurwitzMapFinder@FR", rationalMapFinder@.createMacaulay2HurwitzMapFinder( M2EXEC@ ) );
		fi;
		return HurwitzMapFinder@FR.computeFiniteFieldSolutionsFI( partitionList, strictNormalization, branchValueApproxList, fieldChar );
	end
);


#############################################################################
##
#F SEARCH_HURWITZ_MAP_OVER_FINITE_FIELD( partitionList, strictNormalization, branchValueApproxList, fieldChar )
##
## <#GAPDoc Label="APPROX_HURWITZ_MAP_CANDIDATES">
## <ManSection>
##   <Func Name="APPROX_HURWITZ_MAP_CANDIDATES" Arg=" finiteFieldSolution, liftOptions "/>
##   <Returns>A data structure from which it is possible to construct approximations of rational maps which are lifts 
##		of the input from a finite field Hurwitz map candidate search </Returns>
##   <Description>
##     This function lifts a result from finite field Hurwitz map candidate search (see  <Ref Label="SEARCH_HURWITZ_MAP_OVER_FINITE_FIELD" >) 
##     to an extension field of rational numbers and approximates the map over complex numbers.
##	One finite field search result is lifted to a set of rational maps which are represented by root data:
##	 liftedPolTuple<C>.rootData</C> is a list where each element represents an Huritz map candidate.
##	 A Hurwitz map representation is given by a list with the preimages of the branch values of the map and a scaling factor list.
##	  The first three preimage lists consists of the preimages of (infinity, zero and 1) respectively.
##        Each preimage list consists if pairs of a branch value preimage and its multiplicity.       
##	  The polynomial W_1 is the polynomial with roots as preimages of infinity and W_2 the polynomial with roots as preimages of zero.
##        (_link to preprint_)
##	  The polynomials W_i are constructed via W_2 - scalingFactor[i-2]*W1 for i>2.
##	<P/>
##
##   </Description>
## </ManSection>
## <#/GAPDoc>
rationalMapFinder@.APPROX_HURWITZ_MAP_CANDIDATES   := NewOperation("APPROX_HURWITZ_MAP_CANDIDATES", [IsRecord, IsRecord] );

InstallMethod(
       rationalMapFinder@.APPROX_HURWITZ_MAP_CANDIDATES,
       "approximate Hurwitz map candidates from data obtained by a finite field search (SEARCH_HURWITZ_MAP_OVER_FINITE_FIELD) ", 
       [IsRecord, IsRecord],
	function(  finiteFieldSolution, liftOptions )
		#if (HurwitzMapFinder@FR=fail) then
		#	Unbind(HurwitzMapFinder@FR);
		#fi;
		if not IsBoundGlobal ("HurwitzMapFinder@FR") then
			BindGlobal("HurwitzMapFinder@FR", rationalMapFinder@.createMacaulay2HurwitzMapFinder( M2EXEC@ ) );
		fi;
		return HurwitzMapFinder@FR.approximateComplexSolutions( finiteFieldSolution, liftOptions );
	end
);

# todo : 
# todo : implement a black box rational map where the nominator and denominator are not expanded but represented by products of factors
# this would probably increase precision. 

rationalMapFinder@.CREATE_PRE_RATIONAL_MAP := rationalMapFinder@FR.preRationalMapFromRootDataElem;

