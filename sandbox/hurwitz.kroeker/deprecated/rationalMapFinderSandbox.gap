
# merge RFS  - not for random search and for incomplete characteristics. 
## in case with different precision, take roots with higher precision.
# rfsProblem: shapeList, scalingRelations (alternatively , not yet, minimal Polynomials for relations? )
# RFS: rfsProblem, Options - no!, (Statistics), rfsSetList  ( HTcreate, HTadd; package "orb" )  - at beginning rfsSetTable is empty.
# rfsSolutionList: polIdealList not emply, equationList not empty, rootList not empty.
# rfsSet: link to rfsProblem, ( polynomialRing: link to Ring in [s,t], polynomialList), ( idealRing: link to IdealRing, idealList, equationList, rootList)
# polIdealList - identical for same shape list
# equationList(final) : depends on shape list and scalingRelations
# STREAmse: http://www.gap-system.org/Manuals/doc/htm/ref/CHAP010.htmhttp://www.gap-system.org/Manuals/doc/htm/ref/CHAP010.htm
# http://www.gap-system.org/Manuals/doc/htm/ref/CHAP010.htm

# RECORDS http://www.gap-system.org/Manuals/doc/htm/ref/CHAP027.htm
# Lists http://www.gap-system.org/Gap3/Manual3/C027S000.htm

# types of Objects : http://www.gap-system.org/Manuals/doc/htm/ref/CHAP013.htm#SECT005
# Object example: http://www.cs.st-andrews.ac.uk/~alexk/circle/chap2.html

# creating new Objects: http://www.gap-system.org/Manuals/doc/htm/prg/CHAP003.htm#SSEC017.2
# objects and Elements: http://www.gap-system.org/Manuals/doc/htm/ref/CHAP012.htm#SECT006

# Doku Übersicht: http://www.gap-system.org/Manuals/doc/htm/ref/chapters.htm

# Field: 


# - M2-Aufruf für jeden RFSSet einzeln?

LoadPackage ("orb"); # Hash tables.

# TODO:
# In Future: RFS Object has property searchType with values 'full', 'random'
# Fehlermeldungen aus Macaulay an GAP weiterreichen...


# Structure: RFS; RFS.Problem; RFS.Problem.ShapeList,RFS.Problem.scalingRelations, RFS.PolSetList, 
# Functions: Create RFS (from RFSProblem)
# 

# Questions for gap:
#-------------------
# how to redirect the stderr output from a Process?
# how to delete a file and check for file existence?
# gibtes sowas wie referenzen in GAP?
# howto define own types?  - partly OK
# howto define 
# NULL-value ?
# corresponding ring of a polynomial?
# HashTable in gap? -   HTcreate, package "orb"
# howto define multiple functions in a string and eval them? - InputTextString + Read(input).
# howto get degree of a Multivariate Polynomial  -DegreeMultivariatePolynomial(), package "guave"
# how to get all available Methods for a given ObjectType?
# try catch
# wie funktionier overloading, falls überhaupt
# records: kann man vermeiden dass die Einträge alphabetisch geordnet werden?



# nun ist mir immer noch nicht klar wie minExamples in Zusammenhang mit minChar und maxChar wirken soll
# hört man in der charakteristik auf, in der das letzte benötigte minExample gefunden wurde, oder 
#maxchar oder inexamples muss gesetzt sein - erfüllt wenn minexamples >0.
# glaube es war so: minexamples erreicht=>höre auf;maxchar erreicht =>höre auf.
# wieso nicht maxexamples: wenn maxchar gesetzt  höre bei maxexamples auf; wenn maxChar gesetzt, wird minChar ignoriert oder darf nicht gesetzt sein.
# wenn maxchar nicht gesetzt, höre bei minexamples auf.
#todo: decimalPrecision option!
# isBound
# bin zum Schluss gekommen, dass man nur maxExamples braucht.
 

#rfsProblem: shapeList, scalingRelations
# example:
shape := createShape([4,3,2,2,2]);
IsShape(shape);
shapeList := [shape,shape,shape];
IsShapeList(shapeList);
scalingRelationList:=[];
IsScalingRelationList(scalingRelationList);
rfsProblem:=rec(shapeList:=shapeList ,  scalingRelationList:= scalingRelationList);
IsRFProblem(rfsProblem);



# todo: translate back scaled Factors in solutionPoints when returning result.


coeffRing:=Field( Z(7) );
rnf:=PolynomialRing(coeffRing,["t","s"]); 
tensorRing:=PolynomialRing(rnf,["a1","a2"]); 
ind:=IndeterminatesOfPolynomialRing(rnf);
 t:=ind[1];
 s:=ind[2];
 
 coeffRing2:=Field( Z(7) );
rnf2:=PolynomialRing(coeffRing2,["x","y"]); 
ind:=IndeterminatesOfPolynomialRing(rnf2);
x:=ind[1];
y:=ind[2];
 


LoadPackage("guava"); --DivisorsMultivariatePolynomial # factor multivariate Polynomial!

# x:= commonVariable(polSet);
# y:= homogenVariable(polSet);

DeclareCategory( "IsPolynomialSet", IsMultiplicativeElementWithInverse );



-- howto get a ring when polynomial is given?

IsPolynomialSet := function(polSet)
    local indeterminates;
    local rng;
    if   Length(polSet)<>2 then
    	return false;
    
    rng:=polSet[1];
    
    polList:=polSet[2];
    
    
    
     end;

getSamplePolynomial := function( )
 local indeterminates, coeffRing,rng,pol,t;
 coeffRing:=Field( Z(7) );
 rng:=PolynomialRing(coeffRing,["t" ]); 
 indeterminates :=IndeterminatesOfPolynomialRing(rng);
 t:=indeterminates[1];
 pol:=t+t;
 return pol;
end;

# not correct for the IdealRing...
commonVariable := function(rng)
 local indeterminates;
 indeterminates :=IndeterminatesOfPolynomialRing(rng);
 Assert(2,Length(indeterminates));
 return indeterminates[1];
end;


homogenVariable := function(rng)
 local indeterminates;
 indeterminates:=IndeterminatesOfPolynomialRing(rng);
 Assert(2,Length(indeterminates));
 return indeterminates[2];
end;


--questions: how to get the Ring of a ringElement


KnownAttributesOfObject(pol);
KnownPropertiesOfObject(pol);
-- liefert z.B. NICHT LeadingCoefficient , CoefficientsRing.
CoefficientsRing(rnf);

--get polynomial summands;
-- get polynomial 

--terms: only for MonoidPolynomials...WTF
-- CoefficientsOfUnivariatePolynomial? multivariate pendant

 - package mvp not present or does not work
 
 
coeffRing := Field( Z(7) );
rng:=PolynomialRing( Field(Z(7)) ,["t"] );
ind:=IndeterminatesOfPolynomialRing(rng);
t:=ind[1];
pol:=t+8;
getRing(pol);

Ideal(rnf,[pol,pol]);
GeneratorsOfTwoSidedIdeal

Value(s+t,[s],[1]);

String(pol);
String(rnf);

 EvalString
 polString:=" t:= \"t\"; rng := PolynomialRing( Field(Z(7)) ,[t] ); ";
 
-- funktioniert;
s := "";; str := OutputTextString(s,false);;
PrintTo(str,"rnf:=");
PrintTo(str,rnf);
PrintTo(str,";");
PrintTo(str,bla);

PrintTo(str,"commonVariable1 := ");
PrintTo(str,commonVariable);
PrintTo(str,";\n");
PrintTo(str,"homogenVariable1 :=");
PrintTo(str,homogenVariable);
PrintTo(str,";\n");
input:=InputTextString(s);
Read(input);


GroebnerBasis([pol],MonomialLexOrdering()); 


-
date:= rec(year:=1992, month:="Jan", day:=13);
date.year; 

rnames := RecNames(date);

date.( rnames[1] );

polSet:=rec( polynomialRing:=rnf, degree:=1, polynomialList:=[t,t,t] );


 
    #todo: check if the 	
    
coeffRing:=Field( Z(7) );
rnf:=PolynomialRing(coeffRing,["t","s"]); 
ind:=IndeterminatesOfPolynomialRing(rnf);
 t:=ind[1];
 s:=ind[2];
 
polSet:=Immutable (rec( polynomialRing:=rnf, degree:=1, polynomialList:=[t,t,t] ) );


polSet1:=rec( polynomialRing:=rnf, degree:=1, polynomialList:=[t,t,t],bla:="bla" );
polSet2:=rec( polynomialRing:=rnf, degree:=1  );

DeclareOperation("IsRFSPolynomialSet",[IsObject]);


DeclareCategory( "RFSPolynomialSet", IsRFSPolynomialSet ); 


IN(pol, rng); # das war aber nicht gesucht...


DeclareCategory("IsRFSPolynomialSet",
IsComponentObjectRep);
DeclareOperation("IsShort",[IsRFSPolynomialSet]);
DeclareOperation("NrLetters",[IsBlubb]);


DeclareOperation("IsRFSPolynomialSet",[IsObject]);
InstallMethod(IsRFSPolynomialSet,"for Objects",
[IsObject],IsRFSPolynomialSet1
);




BindGlobal("BlubbsFamily",NewFamily("BlubbsFamily"));

DeclareCategory("IsBlubb",IsComponentObjectRep and IsAttributeStoringRep);


# DeclareCategory("IsBlubb",  IsRFSPolynomialSet); geht nicht...




DeclareRepresentation("IsBlubbDenseRep",
IsBlubb,["polynomialRing","degree","polynomialList"]);

DeclareRepresentation("IsBlubbDenseRep",
IsRFSPolynomialSet,["polynomialRing","degree","polynomialList"]);



BindGlobal("BlubbDenseType",
NewType(BlubbsFamily,IsBlubbDenseRep));


oPolSet:=Objectify(BlubbDenseType,polSet);
oPolSet1:=Objectify(BlubbDenseType,polSet1);
oPolSet2:=Objectify(BlubbDenseType,polSet2);


IsBlubbDenseRep(oPolSet) ; => true
IsBlubbDenseRep(oPolSet1) ; => true
IsBlubbDenseRep(oPolSet2) ; => true (bad);

oPolSet!.polynomialRing;

DeclareOperation("polynomialRing",[IsBlubb]);
DeclareAttribute("polynomialRing1",IsBlubb);
DeclareAttribute("degree",IsBlubb);
DeclareAttribute("tmp",IsBlubb);


InstallMethod(polynomialRing,"for dense Blubbs",
[IsBlubbDenseRep],
function(bl)
return bl!.polynomialRing;
end);

InstallMethod(polynomialRing1,"for dense Blubbs",
[IsBlubbDenseRep],
function(bl)
return bl!.polynomialRing;
end);

InstallMethod(degree,"for dense Blubbs",
[IsBlubbDenseRep],
function(bl)
return bl!.degree;
end);




# knownProperties: Attributes with values true/false.

# 1. Problem:  Setxxx wirkt nur bei erster Anwendung...
# 2. Problem: setzt man die eigenschaft via !.xxx:= , sieht man keine Änderung bei KnownAttributesOfObject(). Und: auf die Eigenschaft kann auch noch nicht zugegriffen werden!

 


