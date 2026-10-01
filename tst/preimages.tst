#############################################################################
##
##  Preimages under maps whose image GAP cannot cheaply test membership in.
##  With checked PreImagesRepresentative these run into coset enumerations
##  or return fail for elements that do have preimages.
##
#############################################################################

gap> START_TEST("preimages.tst");
gap> SetFloats(IEEE754FLOAT);

# SubFRMachine returns fail when a decomposition leaves the subgroup
gap> m := PolynomialSphereMachine(2,[1/7],[]);;
gap> G := StateSet(m);;
gap> H := SphereGroup(4);;
gap> f := GroupHomomorphismByImages(H,G,GeneratorsOfGroup(H),[G.1^2,G.2,G.3,(G.3*G.2*G.1^2)^-1]);;
gap> SubFRMachine(m,f);
fail

# EpimorphismToOut is not surjective as GAP sees it
gap> a := AutomorphismSphereMachine(m);;
gap> Length(GeneratorsOfFRMachine(a));
6

# SubFRMachine and lifting on a 6-generator sphere group
gap> m := Mating(PolynomialSphereMachine(2,[3/7],[]),PolynomialSphereMachine(2,[],[1/6]));;
gap> gens := GeneratorsOfGroup(StateSet(m));;
gap> i := SphereGroup([1,4,3,2,5],[0,0,0,0,0],FreeGroup("f1","f2","f3","g1","x"));;
gap> tm := ChangeFRMachineBasis(m,[gens[1]^-1*gens[5],One(StateSet(m))]);;
gap> inj := GroupHomomorphismByImages(i,StateSet(m),GeneratorsOfGroup(i),[gens[1]^gens[5],gens[2],gens[3],gens[4],gens[1]*gens[6]/gens[1]*gens[5]]);;
gap> m2 := SubFRMachine(tm,inj);
<sphere machine with alphabet [ 1 .. 2 ] on Group( [ f1, f2, f3, g1, x ] ) / [ f1*g1*f3*f2*x ]>
gap> DegreeOfP1Map(P1MapBySphereMachine(m2));
2

#
gap> STOP_TEST("preimages.tst", 1);
