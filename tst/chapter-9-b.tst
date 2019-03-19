#############################################################################
##
#W  chapter-9-b.tst                 FR Package              Laurent Bartholdi
##
#H  @(#)$Id$
##
#Y  Copyright (C) 2011,  Laurent Bartholdi
##
#############################################################################
##
##  This file tests the functions explained in chapter 9 of the manual,
##  that use the DLL
##
#############################################################################

gap> START_TEST("fr:chapter 9 (2/2)");
gap> 
gap> Info(InfoFR,1,"12.8 P1 points");
#I  12.8 P1 points
gap> P1Distance(P1infinity,P1infinity);
0.
gap> 
gap> Info(InfoFR,1,"Shishikura-Tan Lei matings");
#I  Shishikura-Tan Lei matings
gap> SetFloats(IEEE754FLOAT);
gap> z := P1z;
<z>
gap> a := RootsFloat((z-1)*(3*z^2-2*z^3)+1);
[ 0.598631+0.565259i, -0.426536, 0.598631-0.565259i, 1.72927 ]
gap> c := RootsFloat((z^3+z)^3+z);
[ 0., 0.557573+0.540347i, -0.557573+0.540347i, -0.557573-0.540347i,
  0.557573-0.540347i, 0.264425+1.26049i, -0.264425+1.26049i,
  -0.264425-1.26049i, 0.264425-1.26049i ]
gap> am := List(a,a->SphereMachine((a-1)*(3*z^2-2*z^3)+1));
[ <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Grou\
p( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Grou\
p( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Grou\
p( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Grou\
p( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]> ]
gap> cm := List(c,c->SphereMachine(z^3+c));
[ <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f2) on Grou\
p( [ f1, f2 ] )/[ f2*f1 ]>,
  <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Grou\
p( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Grou\
p( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Grou\
p( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Grou\
p( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Grou\
p( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Grou\
p( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Grou\
p( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>,
  <sphere machine with alphabet [ 1 .. 3 ] and adder FRElement(...,f4) on Grou\
p( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]> ]
gap> m := ListX(am,cm,Mating);;
gap> P1MapBySphereMachine(am[2]);
<z^3+(-2.1398032727526402_z)*z+(-0.36023134081173841_z)>
gap> P1MapBySphereMachine(m[9+2]);
<((1.4686858846983573-1.0290716144155767i_z)*z^3+(1.9850056687104112-1.7893204\
8860004i_z)*z^2+(-1.8069445331386629-0.83371844877356505i_z)*z+(-1.66797918021\
99766+0.67701449505226674i_z))/((1.3981038147796185-1.4150160285253968i_z)*z^3\
+(1.8103784467122102-2.3473103431921412i_z)*z^2+(-2.1376242871497211-0.5505179\
5026414962i_z)*z+1._z)>
gap> 
gap> Info(InfoFR,1,"An obstructed mating");
#I  An obstructed mating
gap> P1MapBySphereMachine(m[9+6]);
rec(
  machine := <sphere machine with alphabet [ 1, 2, 3 ] on Group( [ f1, f2, f3, g1, g2, g3 ] ) / [ f\
3*f2*f1*g3*g2*g1 ]>, matrix := [ [ 1/2, 1 ], [ 1/2, 0 ] ],
  multicurve := [ f1*f3*f2*f1*f3*f2*f1*f3*f2*f1*g1*g3*g2*g1*g3*g2*g1*g3*g2*g1^G,
      f1^-1*f2^-1*f1*f3*f2*f1*f3*f2*f1*f3*f2*f1*g2^-1*g3*g2*g1*g3*g2*g1*g3*g2*g1^G ] )
gap> 
gap> Info(InfoFR,1,"Testing Triangulations");
#I  Testing Triangulations
gap> if IsBound(MacFloat) then Float := MacFloat; fi;
gap> oct := List([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[-1.,0.,0.],[0.,-1.,0.],[0.,0.,-1.]],P1Sphere);;
gap> s := Sqrt(Float(1/3));;
gap> cube := List([[s,s,s],[s,s,-s],[s,-s,s],[-s,s,s],[s,-s,-s],[-s,s,-s],[-s,-s,s],[-s,-s,-s]],P1Sphere);;
gap> DelaunayTriangulation(cube);
<triangulation with 11 vertices, 54 edges and 18 faces>
gap> DelaunayTriangulation(cube{[1,5]});
<triangulation with 6 vertices, 24 edges and 8 faces>
gap> p := List([[0.,0.,1.],[0.,0.,-1.],SphereP1(P1Point(1.e-4)),SphereP1(P1Point(0.,1.e4))],P1Sphere);
[ <0+0i>, <P1infinity>, <0.0001+0i>, <0+10000i> ]
gap> DelaunayTriangulation(p,100.);
<triangulation with 32 vertices, 180 edges and 60 faces>
gap> 
gap> Info(InfoFR,1,"Testing RationalFunction");
#I  Testing RationalFunction
gap> f := P1MapBySphereMachine(PolynomialSphereMachine(2,[],[7/16]):param_unicritical);
<z^2+(-1.771257023356245+0.066161509079844388i_z)>
gap> 
gap> Info(InfoFR,1,"Testing Pilgrim's obstructed blowup of the torus");
#I  Testing Pilgrim's obstructed blowup of the torus
gap> F := SphereGroup(4,[0,0,0,0],FreeGroup("a","b","c","d"));
<sphere group on generators [ a, b, c, d ], ordering d*c*b*a>
gap> a := F.1;; b := F.2;; c := F.3;; d := F.4;; o := One(F);;
gap> M := FRMachine(F,[[c^-1,o,o,o,c],[o,o,o,d,d^-1],[a,o,o,a^-1,o],[b,o,d,a,c]],
>                   [(1,5)(2,4,3),(1,2)(4,5),(1,4)(2,3,5),()]);
<sphere machine with alphabet [ 1 .. 5 ] on Group( [ a, b, c, d ] ) / [ d*c*b*\
a ]>
gap> P1MapBySphereMachine(M);
rec(
  machine := <sphere machine with alphabet [ 1, 2, 3, 4, 5 ] on Group( [ a, b,\
 c, d ] ) / [ d*c*b*a ]>, matrix := [ [ 1 ] ], multicurve := [ a*c^G ] )
gap> 
gap> Info(InfoFR,1,"Testing mating of airplane with z^2+i");
#I  Testing mating of airplane with z^2+i
gap> m := Mating(PolynomialSphereMachine(2,[3/7],[]),PolynomialSphereMachine(2,[],[1/6]));
<sphere machine with alphabet [ 1 .. 2 ] on Group( [ f1, f2, f3, g1, g2, g3 ] \
) / [ f3*f2*f1*g3*g2*g1 ]>
gap> Unbind(f1); Unbind(f2); Unbind(f3); Unbind(g1); Unbind(g2); Unbind(g3); 
gap> AssignGeneratorVariables(StateSet(m));
#I  Assigned the global variables [ f1, f2, f3, g1, g2, g3 ]
gap> i := FreeGroup("f1","f2","f3","g1","x");
<free group on the generators [ f1, f2, f3, g1, x ]>
gap> i := SphereGroup([1,4,3,2,5],[0,0,0,0,0],i);
<sphere group on generators [ f1, f2, f3, g1, x ], ordering f1*g1*f3*f2*x>
gap> tm := ChangeFRMachineBasis(m,[f1^-1*g2,One(StateSet(m))]);;
gap> inj := GroupHomomorphismByImages(i,StateSet(m),GeneratorsOfGroup(i),[f1^g2,f2,f3,g1,f1*g3/f1*g2]);
[ f1, f2, f3, g1, x ] -> [ g2^-1*f1*g2, f2, f3, g1, f1*g3*f1^-1*g2 ]
gap> m2 := SubFRMachine(tm,inj);
<sphere machine with alphabet [ 1 .. 2 ] on Group( [ f1, f2, f3, g1, x ] ) / [\
 f1*g1*f3*f2*x ]>
gap> P1MapBySphereMachine(m2);
<((-1.1898473846113464+0.24902829777699695i_z)*z^2+(1.9648756079548171+1.99403\
77105934426i_z)*z+(-0.65573526775326108-2.4418235808459454i_z))/((3.4179914358\
760768+2.610275400738987i_z)*z+1._z)>
gap> 
gap> STOP_TEST( "chapter-9-b.tst", 10^10 );
fr:chapter 9 (2/2)
msecs: 13604

#E chapter-9-b.tst . . . . . . . . . . . . . . . . . . . . . . . . .ends here
