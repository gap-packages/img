#############################################################################
##
#W machine.gd                                               Laurent Bartholdi
##
#Y Copyright (C) 2013, Laurent Bartholdi
##
#############################################################################
##
##  FR machines with stateset a sphere group
##
#############################################################################

#############################################################################
##
#M SphereMachine
##
## <#GAPDoc Label="SphereMachine">
##
DeclareProperty("IsSphereMachine", IsFRMachine);
## <ManSection>
##   <Filt Name="IsSphereMachine" Arg="m"/>
##   <Filt Name="IsPolynomialSphereMachine" Arg="m"/>
##   <Description>
##     The categories of <E>Sphere</E> and <E>polynomial</E> machines.
##     Sphere machines are group FR machines
##     whose underlying group is a sphere group, see <Ref Oper="SphereGroup"/>.
##
##     <P/> A polynomial machine is a group FR machine with a distinguished
##     state (which must be a generator of the stateset), stored as the
##     attribute <Ref Attr="AddingElement"/>; see
##     <Ref Oper="AsPolynomialSphereMachine"/>. If it is normalized, in the sense
##     that the wreath recursion of the adding element <C>a</C> is
##     <C>[[a,1,...,1],[d,1,...,d-1]]</C>, then the basepoint is assumed
##     to be at <M>+\infty</M>; the element <C>a</C> describes a
##     clockwise loop around infinity; the <M>k</M>th preimage of the basepoint
##     is at <M>\exp(2i\pi(k-1)/d)\infty</M>, for <M>k=1,\dots,d</M>; and
##     there is a direct connection from basepoint <M>k</M> to <M>k+1</M> for
##     all <M>k=1,\dots,d-1</M>.
##
##     <P/> The last category is the intersection of the first two.
##   </Description>
## </ManSection>
##
DeclareAttribute("AsSphereMachine", IsGroupFRMachine);
DeclareOperation("AsSphereMachine", [IsGroupFRMachine, IsWord]);
DeclareOperation("AsSphereMachine", [IsGroupFRMachine, IsSphereGroup]);
## <ManSection>
##   <Oper Name="AsSphereMachine" Arg="m[,w]"/>
##   <Returns>A sphere machine.</Returns>
##   <Description>
##     This function creates a new sphere machine, starting from a group
##     FR machine <A>m</A>. If a state <A>w</A> is specified, and that
##     state defines the trivial FR element, then it is used
##     as relator; if <A>w</A> is a sphere group, then it is used as the
##     new stateset.
##     Finally, if no relator and no group is specified, and the product (in some ordering)
##     of the generators is trivial, then that product is used as
##     relator. In other cases, the method returns <K>fail</K>.
##
##     <P/> A standard FR machine can be recovered from a sphere machine
##     by <Ref Oper="AsGroupFRMachine" BookName="FR"/>,
##     <Ref Oper="AsMonoidFRMachine" BookName="FR"/>,
##     and <Ref Oper="AsSemigroupFRMachine" BookName="FR"/>.
## <Example><![CDATA[
## gap> m := UnderlyingFRMachine(BasilicaGroup);
## <Mealy machine on alphabet [ 1 .. 2 ] with 3 states>
## gap> g := AsGroupFRMachine(m);
## <FR machine with alphabet [ 1 .. 2 ] on Group( [ f1, f2 ] )>
## gap> AsSphereMachine(g,Product(GeneratorsOfFRMachine(g)));
## <FR machine with alphabet [ 1 .. 2 ] on Group( [ f1, f2, t ] )/[ f1*f2*t ]>
## gap> Display(last);
##  G  |              1         2
## ----+-----------------+---------+
##  f1 |          <id>,2      f2,1
##  f2 |          <id>,1      f1,2
##   t | f2^-1*f1*f2*t,2   f1^-1,1
## ----+-----------------+---------+
## Relator: f1*f2*t
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareGlobalFunction("NewSphereMachine");
## <ManSection>
##   <Attr Name="NewSphereMachine" Arg="..."/>
##   <Returns>A new sphere machine, based on string descriptions.</Returns>
##   <Description>
##     This command constructs a new sphere machine, in a format similar to
##     <Ref Func="FRGroup" BookName="FR"/>; namely, the arguments are strings of the form
##     "gen=&lt;word-1,...,word-d&gt;perm"; each <C>word-i</C> is a word in the
##     generators; and <C>perm</C> is a transformation,
##     either written in disjoint cycle or in images notation. The underlying
##     group of the machine is a sphere group.
##
##     <P/><C>word-i</C> is allowed to be the
##     empty string; and the "&lt;...&gt;" may be skipped altogether.
##     Each <C>word-i</C> may also contain inverses.
##
##     <P/>The extra final arguments describe relations in the underlying
##     sphere group; at least one relation is required, the product of the
##     generators in an appropriate order.
##         
##     <P/>The following examples construct realizable foldings of the
##     polynomial <M>z^3+i</M>, following Cui's arguments.         
## <Example><![CDATA[
## gap> fold1 := NewSphereMachine("a=<,,b,,,B>(1,2,3)(4,5,6)","b=<,,b*a/b,,,B*A/B>",
##      "A=<,,b*a,,,B*A>(3,6)","B=(1,6,5,4,3,2)","a*B*A*b");
## gap> <FR machine with alphabet [ 1, 2, 3, 4, 5, 6 ] on Group( [ a, b, A, B ] )/[ a*B*A*b ]>                                
## gap> fold2 := NewSphereMachine("a=<,,b,,,B>(1,2,3)(4,5,6)","b=<,,b*a/b,,,B*A/B>",
##      "A=(1,6)(2,5)(3,4)","B=<B*A,,,b*a,,>(1,4)(2,6)(3,5)","a*B*A*b");;
## gap> P1MapBySphereMachine(fold1); P1MapBySphereMachine(fold2);
## ...
## ]]></Example>
##   </Description>
## </ManSection>
## <#/GAPDoc>
##
#############################################################################

#############################################################################
##
#M PolynomialMachine
##
## <#GAPDoc Label="PolynomialFRMachine">
DeclareProperty("IsKneadingMachine",IsFRMachine);
DeclareProperty("IsPlanarKneadingMachine",IsFRMachine);
InstallTrueMethod(IsBoundedFRMachine,IsKneadingMachine);
InstallTrueMethod(IsLevelTransitive,IsKneadingMachine);
## <ManSection>
##   <Prop Name="IsKneadingMachine" Arg="m"/>
##   <Prop Name="IsPlanarKneadingMachine" Arg="m"/>
##   <Returns>Whether <A>m</A> is a (planar) kneading machine.</Returns>
##   <Description>
##     A <E>kneading machine</E> is a special kind of Mealy machine, used
##     to describe postcritically finite complex polynomials. It is a
##     machine such that its set of permutations is "treelike" (see
##     <Cite Key="MR2162164" Where="§6.7"/>) and such that each non-trivial
##     state occurs exactly once among the outputs.
##
##     <P/> Furthermore, this set of permutations is <E>treelike</E> if
##     there exists an ordering of the states that their product in that
##     order <M>t</M> is an adding machine; i.e. such that <M>t</M>'s
##     activity is a full cycle, and the product of its states along that
##     cycle is conjugate to <M>t</M>. This element <M>t</M> represents the
##     Carathéodory loop around infinity.
## <Example><![CDATA[
## gap> M := BinaryKneadingMachine("0");
## BinaryKneadingMachine("0*")
## gap> Display(M);
##    |  1     2
## ---+-----+-----+
##  a | c,2   b,1
##  b | a,1   c,2
##  c | c,1   c,2
## ---+-----+-----+
## gap> IsPlanarKneadingMachine(M);
## true
## gap> IsPlanarKneadingMachine(GrigorchukMachine);
## false
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareAttribute("AddingElement", IsSphereMachine);
## <ManSection>
##   <Attr Name="AddingElement" Arg="m"/>
##   <Returns>The element generating the adding submachine.</Returns>
##   <Description>
##     This attribute stores the product of generators that is an
##     adding machine.
##     In essence, it records an ordering of the generators whose
##     product corresponds to the Carathéodory loop around infinity.
##
##     <P/> The following example illustrates Wittner's shared mating
##     of the airplane and the rabbit. In the machine <C>m</C>, an
##     airplane is represented by <C>Group(a,b,c)</C> and a rabbit is
##     represented by <C>Group(x,y,z)</C>; in the machine <C>newm</C>,
##     it is the other way round.
## <Example><![CDATA[
## gap> f := FreeGroup("a","b","c","x","y","z");;
## gap> AssignGeneratorVariables(f);
## gap> m := AsSphereMachine(FRMachine(f,[[a^-1,b*a],[One(f),c],[a,One(f)],[z*y*x,
##        x^-1*y^-1],[One(f),x],[One(f),y]],[(1,2),(),(),(1,2),(),()]));;
## gap> Display(m);
##  G |      1             2   
## ---+---------+-------------+
##  a |  a^-1,2         b*a,1  
##  b |  <id>,1           c,2  
##  c |     a,1        <id>,2  
##  x | z*y*x,2   x^-1*y^-1,1  
##  y |  <id>,1           x,2  
##  z |  <id>,1           y,2  
## ---+---------+-------------+
## Relator: z*y*x*c*b*a
## gap> iso := GroupHomomorphismByImages(f,f,[a,b^(y^-1),c^(x^-1*y^-1*a^-1),x^(b*a*z*a^-1),y,z^(a^-1)],[a,b,c,x,y,z]);;
## gap> newm := ChangeFRMachineBasis(m^iso,[a^-1*y^-1,y^-1*a^-1*c^-1]);;
## gap> Display(newm);
##  G |          1         2   
## ---+-------------+---------+
##  a | a^-1*c^-1,2   c*a*b,1  
##  b |      <id>,1       c,2  
##  c |         a,1    <id>,2  
##  x |       z*x,2    x^-1,1  
##  y |      <id>,1       x,2  
##  z |         y,1    <id>,2  
## ---+-------------+---------+
## Relator: c*a*b*y*z*x
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareSynonym("IsPolynomialSphereMachine",IsSphereMachine and HasAddingElement);
DeclareAttribute("AsPolynomialSphereMachine",IsFRMachine);
DeclareOperation("AsPolynomialSphereMachine",[IsFRMachine,IsWord]);
## <ManSection>
##   <Oper Name="AsPolynomialSphereMachine" Arg="m [,adder [,relator]]"/>
##   <Returns>A polynomial sphere machine.</Returns>
##   <Description>
##     The first function creates a new polynomial sphere machine, starting from
##     a group or Mealy machine. A <E>polynomial</E> machine is one that
##     has a distinguished adding element, <Ref Attr="AddingElement"/>.
##
##     <P/> If the argument is a Mealy machine, it must be planar (see
##     <Ref Prop="IsPlanarKneadingMachine"/>). If the argument is a group
##     machine, its permutations must be treelike, and its outputs must be
##     such that, up to conjugation, each non-trivial state appears
##     exactly once as the product along all cycles of all states.
##
##     <P/> If a second argument <A>adder</A> is supplied, it is checked to
##     represent an adding element, and is used as such.
##
## <Example><![CDATA[
## gap> M := PolynomialMealyMachine(2,[1/7],[]);
## <Mealy machine on alphabet [ 1 .. 2 ] with 4 states>
## gap> Mi := AsPolynomialSphereMachine(M);
## !!!
## ]]></Example>
##   </Description>
## </ManSection>
## <#/GAPDoc>
##
#############################################################################

#############################################################################
##
#M Operations
##
## <#GAPDoc Label="SphereOperations">
DeclareOperation("PolynomialMealyMachine",[IsPosInt,IsList,IsList]);
DeclareOperation("PolynomialMealyMachine",[IsPosInt,IsList]);
DeclareOperation("PolynomialSphereMachine",[IsPosInt,IsList,IsList,IsRecord]);
DeclareOperation("PolynomialSphereMachine",[IsPosInt,IsList,IsList]);
DeclareOperation("PolynomialSphereMachine",[IsPosInt,IsList,IsRecord]);
DeclareOperation("PolynomialSphereMachine",[IsPosInt,IsList]);
## <ManSection>
##   <Oper Name="PolynomialSphereMachine" Arg="d,per[,pre][,options]"/>
##   <Oper Name="PolynomialMealyMachine" Arg="d,per[,pre]"/>
##   <Returns>A sphere or Mealy machine.</Returns>
##   <Description>
##     This function creates a sphere or Mealy machine that describes
##     a topological polynomial. The polynomial is described symbolically
##     in the language of <E>external angles</E>. For more details, see
##     <Cite Key="MR762431"/> and <Cite Key="MR812271"/> (in the quadratic
##      case), <Cite Key="MR1149891"/> (in the preperiodic case), and
##     <Cite Key="math.DS/9305207"/> (in the general case).
##
##     <P/> <A>d</A> is the degree of the polynomial. <A>per</A> and
##     <A>pre</A> are lists of angles or preangles. In what follows,
##     angles are rational numbers, considered modulo 1.  Each entry in
##     <A>per</A> or <A>pre</A> is either a rational (interpreted as an
##     angle), or a list of angles <M>[a_1,\ldots,a_i]</M> such that
##     <M>da_1=\ldots=da_i</M>. The angles in <A>per</A> are angles landing
##     at the root of a Fatou component, and the angles in <A>pre</A> land
##     on the Julia set.
##
##     <P/> Note that, for sphere machines, the last generator of the machine
##     produced is an adding machine, representing a loop going
##     counterclockwise around infinity (in the compactification of
##     <M>\mathbb C</M> by a disk, this loop goes <E>clockwise</E> around
##     that disk).
##
##     <P/> In constructing a polynomial sphere machine, one may specify a
##     record <A>options</A>, which may contain the following fields:
##     <C>mealy</C> (boolean, default <K>false</K>) specifies if a formal
##     construction is required; <C>adding</C> specifying that the adding
##     machine should have the most compact representation; and
##     <C>orbispace</C> (boolean, default <K>false</K>) asking the
##     constructed group to have orbispace points of minimal degree.
##
##     <P/> In
##     a <E>formal</E> recursion, distinct angles give distinct generators;
##     while in a non-formal recursion, distinct angles, which land at the
##     same point in the Julia set, give a single generator. The simplest
##     example where this occurs is angle <M>5/12</M> in the quadratic
##     family, in which angles <M>1/3</M> and <M>2/3</M> land at the same
##     point -- see the example below.
##
##     <P/> The attribute <C>Correspondence(m)</C> records the angles
##     landing on the generators: <C>Correspondence(m)[i]</C> is a list
##     <C>[a,s]</C> where <M>a</M> is an angle landing on generator <C>i</C>
##     and <M>s</M> is <K>"Julia"</K> or <K>"Fatou"</K>.
##
##     <P/> If only one list of angles is supplied, then <Package>IMG</Package>
##     guesses that all angles with denominator coprime to <A>n</A> are
##     Fatou, and all the others are Julia.
##
##     <P/> The inverse operation, reconstructing the angles from the sphere
##     machine, is <Ref Oper="SupportingRays"/>.
## <Example><![CDATA[
## gap> PolynomialSphereMachine(2,[0],[]); # the adding machine
## <FR machine with alphabet [ 1 .. 2 ] on Group( [ f1, f2 ] )/[ f2*f1 ]>
## gap> Display(last);
##  G  |     1        2
## ----+--------+--------+
##  f1 | <id>,2     f1,1
##  f2 |   f2,2   <id>,1
## ----+--------+--------+
## Relator: f2*f1
## gap> Display(PolynomialSphereMachine(2,[1/3],[])); # the Basilica
##  G  |      1         2
## ----+---------+---------+
##  f1 | f1^-1,2   f2*f1,1
##  f2 |    f1,1    <id>,2
##  f3 |    f3,2    <id>,1
## ----+---------+---------+
## Relator: f3*f2*f1
## gap> Display(PolynomialSphereMachine(2,[],[1/6])); # z^2+I
##  G  |            1         2
## ----+---------------+---------+
##  f1 | f1^-1*f2^-1,2   f2*f1,1
##  f2 |          f1,1      f3,2
##  f3 |          f2,1    <id>,2
##  f4 |          f4,2    <id>,1
## ----+---------------+---------+
## Relator: f4*f3*f2*f1
## gap> PolynomialSphereMachine(2,[],[5/12]);
## <FR machine with alphabet [ 1, 2 ] and adder f4 on Group( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>
## gap> Correspondence(last);
## [ [ [ 1/3, 2/3 ], "Julia" ], [ [ 5/12 ], "Julia" ], [ [ 5/6 ], "Julia" ] ]
## gap> PolynomialSphereMachine(2,[],[5/12],rec(formal:=true));
## <FR machine with alphabet [ 1, 2 ] and adder f5 on Group( [ f1, f2, f3, f4, f5 ] )/[ f5*f4*f3*f2*f1 ]>
## gap> Correspondence(last);
## [ [ 1/3, "Julia" ], [ 5/12, "Julia" ], [ 2/3, "Julia" ], [ 5/6, "Julia" ] ]
## ]]></Example>
##     The following construct the examples in Poirier's paper:
## <Listing><![CDATA[
## PoirierExamples := function(arg)
##     if arg=[1] then
##         return PolynomialSphereMachine(2,[1/7],[]);
##     elif arg=[2] then
##         return PolynomialSphereMachine(2,[],[1/2]);
##     elif arg=[3,1] then
##         return PolynomialSphereMachine(2,[],[5/12]);
##     elif arg=[3,2] then
##         return PolynomialSphereMachine(2,[],[7/12]);
##     elif arg=[4,1] then
##         return PolynomialSphereMachine(3,[[3/4,1/12],[1/4,7/12]],[]);
##     elif arg=[4,2] then
##         return PolynomialSphereMachine(3,[[7/8,5/24],[5/8,7/24]],[]);
##     elif arg=[4,3] then
##         return PolynomialSphereMachine(3,[[1/8,19/24],[3/8,17/24]],[]);
##     elif arg=[5] then
##         return PolynomialSphereMachine(3,[[3/4,1/12],[3/8,17/24]],[]);
##     elif arg=[6,1] then
##         return PolynomialSphereMachine(4,[],[[1/4,3/4],[1/16,13/16],[5/16,9/16]]);
##     elif arg=[6,2] then
##         return PolynomialSphereMachine(4,[],[[1/4,3/4],[3/16,15/16],[7/16,11/16]]);
##     elif arg=[7] then
##         return PolynomialSphereMachine(5,[[0,4/5],[1/5,2/5,3/5]],[[1/5,4/5]]);
##     elif arg=[9,1] then
##         return PolynomialSphereMachine(3,[[0,1/3],[5/9,8/9]],[]);
##     elif arg=[9,2] then
##         return PolynomialSphereMachine(3,[[0,1/3]],[[5/9,8/9]]);
##     else
##         Error("Unknown Poirier example ",arg);
##     fi;
## end;
## ]]></Listing>
##   </Description>
## </ManSection>
##
DeclareAttribute("SupportingRays",IsFRMachine);
## <ManSection>
##   <Attr Name="SupportingRays" Arg="m"/>
##   <Returns>A <C>[degree,fatou,julia]</C> description of <A>m</A>.</Returns>
##   <Description>
##     This operation is the inverse of <Ref Oper="PolynomialSphereMachine"/>:
##     it computes a choice of angles, describing landing rays on Fatou/Julia
##     critical points.
##
##     <P/> If there does not exist a complex realization, namely if the
##     machine is obstructed, then this command returns an obstruction, as
##     a record. The field <K>minimal</K> is set to false, and a proper
##     sub-machine is set as the field <K>submachine</K>. The field
##     <K>homomorphism</K> gives an embedding of the stateset of
##     <K>submachine</K> into the original machine, and <K>relation</K> is
##     the equivalence relation on the set of generators of <A>m</A> that
##     describes the pinching.
## <Example><![CDATA[
## gap> r := PolynomialSphereMachine(2,[1/7],[]);
## <FR machine with alphabet [ 1, 2 ] and adder f4 on Group( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>
## gap> F := StateSet(r);; SetName(F,"F");
## gap> SupportingRays(r);
## [ 2, [ [ 1/7, 9/14 ] ], [  ] ] # actually returns the angle 2/7
## gap> # now CallFuncList(PolynomialSphereMachine,last) would return the machine r
## gap> twist := GroupHomomorphismByImages(F,F,GeneratorsOfGroup(F),[F.1^(F.2*F.1),F.2^F.1,F.3,F.4])^-1;
## [ f1, f2, f3, f4 ] -> [ f1*f2*f1^-1, f2*f1*f2*f1^-1*f2^-1, f3, f4 ]
## gap> List([-5..5],i->2*SupportingRays(r*twist^i)[2][1][1]);
## [ 4/7, 5/7, 4/7, 4/7, 5/7, 2/7, 4/7, 4/7, 2/7, 4/7, 4/7 ]
## gap> r := PolynomialSphereMachine(2,[],[1/6]);;
## gap> F := StateSet(r);;
## gap> twist := GroupHomomorphismByImages(F,F,GeneratorsOfGroup(F),[F.1,F.2^(F.3*F.2),F.3^F.2,F.4]);;
## gap> SupportingRays(r);
## [ 2, [  ], [ [ 1/12, 7/12 ] ] ]
## gap> SupportingRays(r*twist);
## [ 2, [  ], [ [ 5/12, 11/12 ] ] ]
## gap> SupportingRays(r*twist^2);
## rec(
##   transformation := [ [ f1, f2^-1*f3^-1*f2^-1*f3^-1*f2*f3*f2*f3*f2, f2^-1*f3^-1*f2^-1*f3*f2*f3*f2,
##           f4 ] -> [ f1, f2, f3, f4 ],
##       [ f1^-1*f2^-1*f1^-1*f2^-1*f1*f2*f1*f2*f1, f1^-1*f2^-1*f1^-1*f2*f1*f2*f1, f3, f4 ] ->
##         [ f1, f2, f3, f4 ],
##       [ f1^-1*f2^-1*f3^-1*f2*f1*f2^-1*f3*f2*f1, f2, f2*f1^-1*f2^-1*f3*f2*f1*f2^-1, f4 ] ->
##         [ f1, f2, f3, f4 ], [ f1, f3*f2*f3^-1, f3, f4 ] -> [ f1, f2, f3, f4 ],
##       [ f1, f2, f2*f3*f2^-1, f4 ] -> [ f1, f2, f3, f4 ],
##       [ f1, f3*f2*f3^-1, f3, f4 ] -> [ f1, f2, f3, f4 ],
##       [ f1, f2, f2*f3*f2^-1, f4 ] -> [ f1, f2, f3, f4 ],
##       [ f1, f3*f2*f3^-1, f3, f4 ] -> [ f1, f2, f3, f4 ] ], machine := <FR machine with alphabet
##     [ 1, 2 ] and adder f4 on Group( [ f1, f2, f3, f4 ] )/[ f4*f3*f2*f1 ]>, minimal := false,
##   submachine := <FR machine with alphabet [ 1, 2 ] and adder f3 on Group( [ f1, f2, f3 ] )>,
##   homomorphism := [ f1, f2, f3 ] -> [ f1, f2*f3, f4 ],
##   relation := <equivalence relation on <object> >, niter := 8 )
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareAttribute("SphereMachineWithNormalizedAdder",IsSphereMachine);
DeclareOperation("SphereMachineWithNormalizedAdder",[IsSphereMachine,IsAssocWord]);
## <ManSection>
##   <Attr Name="NormalizedAdderSphereMachine" Arg="m [adder]"/>
##   <Returns>A sphere machine.</Returns>
##   <Description>
##     This function returns a new sphere machine, in which the adding element
##     has been put into a standard form <M>t=[t,1,\dots,1]s</M>, where
##     <M>s</M> is the long cycle <M>i\mapsto i-1</M>. If the adding element
##     <A>adder</A> is not specified, then <A>m</A> should be a polynomial
##     sphere machine, and <A>adder</A> is its <C>AddingElement</C>.
##   </Description>
## </ManSection>
##
DeclareAttribute("SimplifiedSphereMachine",IsSphereMachine);
## <ManSection>
##   <Attr Name="SimplifiedSphereMachine" Arg="m"/>
##   <Returns>A simpler sphere machine.</Returns>
##   <Description>
##     This function returns a new sphere machine, with hopefully simpler
##     transitions. The simplified machine is obtained by applying
##     automorphisms to the stateset. The sequence of automorphisms
##     (in increasing order) is stored as a correspondence; namely,
##     if <C>n=SimplifiedSphereMachine(m)</C>, then
##     <C>m^Product(Correspondence(n))=n</C>.
## <Example><![CDATA[
## gap> r := PolynomialSphereMachine(2,[1/7],[]);;
## gap> F := StateSet(r);; SetName(F,"F");
## gap> twist := GroupHomomorphismByImages(F,F,GeneratorsOfGroup(F),[F.1^(F.2*F.1),F.2^F.1,F.3,F.4]);;
## gap> m := r*twist;; Display(m);
##  G  |                     1            2
## ----+------------------------+------------+
##  f1 |          f1^-1*f2^-1,2   f3*f2*f1,1
##  f2 | f1^-1*f2^-1*f1*f2*f1,1       <id>,2
##  f3 |          f1^-1*f2*f1,1       <id>,2
##  f4 |                   f4,2       <id>,1
## ----+------------------------+------------+
## Adding element: f4
## Relator: f4*f3*f2*f1
## gap> n := SimplifiedSphereMachine(m);
## <FR machine with alphabet [ 1, 2 ] and adder f4 on F>
## gap> Display(n);
##  G  |            1            2
## ----+---------------+------------+
##  f1 | f2^-1*f1^-1,2   f1*f2*f3,1
##  f2 |        <id>,1         f1,2
##  f3 |        <id>,1         f2,2
##  f4 |          f4,2       <id>,1
## ----+---------------+------------+
## Adding element: f4
## Relator: f4*f1*f2*f3
## gap> n = m^Product(Correspondence(n));
## true
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareOperation("Mating",[IsPolynomialSphereMachine,IsPolynomialSphereMachine]);
DeclareOperation("Mating",[IsPolynomialSphereMachine,IsPolynomialSphereMachine,IsBool]);
DeclareAttribute("EquatorElement",IsSphereMachine);
DeclareAttribute("EquatorTwist",IsSphereMachine);
## <ManSection>
##   <Oper Name="Mating" Arg="m1,m2 [,formal]"/>
##   <Attr Name="EquatorElement" Arg="m"/>
##   <Attr Name="EquatorTwist" Arg="m"/>
##   <Returns>A sphere machine.</Returns>
##   <Description>
##     This function "mates" two polynomial sphere machines.
##
##     <P/> The mating is defined as follows: one removes a disc around
##     the adding machine in <A>m1</A> and <A>m2</A>; one applies complex
##     conjugation to <A>m2</A>; and one glues the hollowed spheres along
##     their boundary circle.
##
##     <P/> The optional argument <A>formal</A>, which defaults to
##     <K>true</K>, specifies whether a <E>formal</E> mating should be done;
##     in a non-formal mating, generators of <A>m1</A> and <A>m2</A> which
##     have identical angle should be treated as a single generator. A
##     non-formal mating is of course possible only if the machines are
##     realizable -- see <Ref Oper="SupportingRays"/>.
##
##     <P/> The attribute <C>Correspondence</C> is a pair of homomorphisms,
##     from the statesets of <A>m1,m2</A> respectively to the stateset of the
##     mating.
##
##     <P/> The attribute <C>EquatorElement</C> is set, and records the
##     original adding elements of <A>m1,m2</A>, which have become the equator
##     of the mating.
##
##     <P/> Note that there are <M>d-1</M> different matings between
##     polynomials of degree <M>d</M>: each has <M>d-1</M> fixed rays at
##     angles <M>2\pi ik/(d-1)</M>. This command constructs the mating in
##     which rays at angle <M>0</M> are matched to each other. To obtain the
##     other matings, multiply the machine by a power of its
##     <C>EquatorTwist</C>.
## <Example><![CDATA[
## gap> # the Tan-Shishikura examples
## gap> SetP1Points(PMCOMPLEX);
## gap> z := Indeterminate(@IMG.field);;
## gap> a := RootsFloat((z-1)*(3*z^2-2*z^3)+1);;
## gap> c := RootsFloat((z^2+1)^3*z^2+1);;
## gap> am := List(a,a->SphereMachine((a-1)*(3*P1z^2-2*P1z^3)+1));;
## gap> cm := List(c,c->SphereMachine(P1z^3+c));;
## gap> m := ListX(am,cm,Mating);;
## gap> # m[1] is realizable
## gap> P1MapBySphereMachine(m[1]);
## ((1.66408+I*0.668485)*z^3+(-2.59772+I*0.627498)*z^2+(-1.80694-I*0.833718)*z
##   +(1.14397-I*1.38991))/((-1.52357-I*1.27895)*z^3+(2.95502+I*0.234926)*z^2
##   +(1.61715+I*1.50244)*z+1)
## gap> # m[29] is obstructed, and is the original Tan-Shishikura map
## gap> P1MapBySphereMachine(m[29]);
## rec(
##   machine := <sphere machine with alphabet [ 1, 2, 3 ] on Group( [ f1, f2, f3, g1, g2,\
##  g3 ] ) / [ f3*f2*f1*g3*g2*g1 ]>, matrix := [ [ 1, 1/2 ], [ 1, 0 ] ],
##   multicurve := [ (f1*f3*f2)^2*f1*g1^-1*(g3*g2*g1)^3*f2^-1*f3*f2^G,
##       f1^-1*f2^-1*(f3*f2*f1)^4*g2*(g3*g2*g1)^3^G ] )
## gap> but the other mating of the same polynomials is not obstructed:
## gap> P1MapBySphereMatrix(m[29]*EquatorTwist(m[29]));
## <((-1.4495156808145406+0.44648102591936722i_z)*z^3+
## (-1.1286550578708263-0.40162285610021786i_z)*z^2+
## (1.0326873942952213-0.11770300021984977i_z)*z+
## (1.0940864612174037+0.24650956710141259i_z))/
## ((0.85917327990384307-0.8755042485587835i_z)*z^3+
## (0.9573881709899621-0.14875521653926685i_z)*z^2+
## (-0.68923589444039035+0.48120812618585479i_z)*z+1._z)>
## ]]></Example>
##   </Description>
## </ManSection>
## <#/GAPDoc>
##
#############################################################################

#############################################################################
##
#M Automorphisms of sphere machines
##
## <#GAPDoc Label="AutomorphismsFRMachines">
DeclareOperation("AutomorphismVirtualEndomorphism",[IsGroupHomomorphism]);
DeclareOperation("AutomorphismSphereMachine",[IsSphereMachine]);
## <ManSection>
##   <Attr Name="AutomorphismVirtualEndomorphism" Arg="v"/>
##   <Attr Name="AutomorphismSphereMachine" Arg="m"/>
##   <Returns>A description of the pullback map on Teichmüller space.</Returns>
##   <Description>
##     Let <A>m</A> be a sphere machine, thought of as a biset for the
##     fundamental group <M>G</M> of a punctured sphere. Let <M>M</M> denote
##     the automorphism of the surface, seen as a group of outer
##     automorphisms of <M>G</M> that fixes the conjugacy classes of punctures.
##
##     <P/> Choose an alphabet letter <A>a</A>, and consider the
##     virtual endomorphism <M>v:G_a\to G</M>. Let <M>H</M> denote the
##     subgroup of <M>M</M> that fixes all conjugacy classes of <M>G_a</M>.
##     then there is an induced virtual endomorphism <M>\alpha:H\to M</M>,
##     defined by <M>t^\alpha=v^{-1}tv</M>. This is the homomorphism
##     computed by the first command. Its source and range are in fact
##     groups of automorphisms of range of <A>v</A>.
##
##     <P/> The second command constructs an FR machine associated with
##     <A>\alpha</A>. Its stateset is a free group generated by elementary
##     Dehn twists of the generators of <A>G</A>.
##
## <Example><![CDATA[
## gap> SetP1Points(PMCOMPLEX);
## gap> z := Indeterminate(@IMG.field);;
## gap> # a Sierpinski carpet map without multicurves
## gap> m := SphereMachine((z^2-z^-2)/2/COMPLEX_I);
## <FR machine with alphabet [ 1, 2, 3, 4 ] on Group( [ f1, f2, f3, f4 ] )/[ f3*f2*f1*f4 ]>
## gap> AutomorphismSphereMachine(i);
## <FR machine with alphabet [ 1, 2 ] on Group( [ x1, x2, x3, x4, x5, x6 ] )>
## gap> Display(last);
##  G  |     1        2
## ----+--------+--------+
##  x1 | <id>,2   <id>,1  
##  x2 | <id>,1   <id>,2  
##  x3 | <id>,2   <id>,1  
##  x4 | <id>,2   <id>,1  
##  x5 | <id>,1   <id>,2  
##  x6 | <id>,2   <id>,1  
## ----+--------+--------+
## gap> # the original rabbit problem
## gap> m := PolynomialSphereMachine(2,[1/7],[]);;
## gap> v := VirtualEndomorphism(m,1);;
## gap> a := AutomorphismVirtualEndomorphism(v);
## MappingByFunction( <group with 20 generators>, <group with 6 generators>, function( a ) ... end )
## gap> Source(a).1;
## [ f1, f2, f3, f4 ] -> [ f3*f2*f1*f2^-1*f3^-1, f2, f3, f3*f2*f1^-1*f2^-1*f3^-1*f2^-1*f3^-1 ]
## gap> Image(a,last);
## [ f1, f2, f3, f4 ] -> [ f1, f2, f2*f1*f3*f1^-1*f2^-1, f3^-1*f1^-1*f2^-1 ]
## gap> # so last2*m is equivalent to m*last
## ]]></Example>
##   </Description>
## </ManSection>
## <#/GAPDoc>
##
#############################################################################

#############################################################################
##
## <#GAPDoc Label="RotatedSpider">
DeclareOperation("LiftOfConjugacyClass", [IsGroupFRMachine,IsConjugacyClassGroupRep]);
## <ManSection>
##   <Oper Name="LiftOfConjugacyClass" Arg="m c"/>
##   <Returns>A list of conjugacy classes and multiplicities.</Returns>
##   <Description>
##     This command computes the preimage of the conjugacy class <A>c</A>
##     by the sphere machine <A>m</A>, namely, it applies the wreath recursion
##     to a representative of <A>c</A> and collects the products on all
##     cycles. It returns then a list of pairs <C>[cc,len]</C> where <C>cc</C>
##     is the conjugacy class of a product on a cycle, and <C>len</C> is the
##     length of the cycle.
##   </Description>
## </ManSection>
##
#DeclareAttribute("ComplexConjugate", IsFRMachine); # already declared for arithmetic objects
## <ManSection>
##   <Oper Name="ComplexConjugate" Arg="m"/>
##   <Returns>An FR machine with inverted states.</Returns>
##   <Description>
##     This function constructs an FR machine whose generating states are
##     the inverses of the original states. If <A>m</A> came from a complex
##     rational map <M>f(z)</M>, this would construct the machine of the
##     conjugate map <M>\overline{f(\overline z)}</M>.
## <Example><![CDATA[
## gap> a := PolynomialSphereMachine(2,[1/7]);
## <FR machine with alphabet [ 1, 2 ] and adder FRElement(...,f4) on <object>/[ f4*f3*f2*f1 ]>
## gap> Display(a);
##  G  |            1            2
## ----+---------------+------------+
##  f1 | f1^-1*f2^-1,2   f3*f2*f1,1
##  f2 |          f1,1       <id>,2
##  f3 |          f2,1       <id>,2
##  f4 |          f4,2       <id>,1
## ----+---------------+------------+
## Adding element: FRElement(...,f4)
## Relator: f4*f3*f2*f1
## gap> Display(ComplexConjugate(a));
##  G  |            1                     2
## ----+---------------+---------------------+
##  f1 | f1*f2*f3*f4,2   f4^-1*f2^-1*f1^-1,1
##  f2 |          f1,1      <identity ...>,2
##  f3 |          f2,1      <identity ...>,2
##  f4 |          f4,2      <identity ...>,1
## ----+---------------+---------------------+
## Adding element: FRElement(...,f4)
## Relator: f1*f2*f3*f4
## gap> ExternalAngle(a);
## {2/7}
## gap> ExternalAngle(ComplexConjugate(a));
## {6/7}
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareOperation("RotatedSpider", [IsPolynomialSphereMachine]);
DeclareOperation("RotatedSpider", [IsPolynomialSphereMachine, IsInt]);
## <ManSection>
##   <Oper Name="RotatedSpider" Arg="m, [p]"/>
##   <Returns>A polynomial FR machine with rotated spider at infinity.</Returns>
##   <Description>
##     This function constructs an isomorphic polynomial FR machine, but with
##     a different numbering of the spider legs at infinity. This rotation is
##     accomplished by conjugating by <C>adder^p</C>, where <C>adder</C> is the
##     adding element of <A>m</A>, and <A>p</A>, the rotation parameter, is
##     <M>1</M> by default.
## <Example><![CDATA[
## gap> a := PolynomialSphereMachine(3,[1/4]);
## <FR machine with alphabet [ 1, 2, 3 ] and adder FRElement(...,f3) on <object>/[ f3*f2*f1 ]>
## gap> Display(a);
##  G  |      1        2         3
## ----+---------+--------+---------+
##  f1 | f1^-1,2   <id>,3   f2*f1,1
##  f2 |    f1,1   <id>,2    <id>,3
##  f3 |    f3,3   <id>,1    <id>,2
## ----+---------+--------+---------+
## Adding element: FRElement(...,f3)
## Relator: f3*f2*f1
## gap> Display(RotatedSpider(a));
##  G  |     1            2               3
## ----+--------+------------+---------------+
##  f1 | <id>,2   f2*f1*f3,3   f3^-1*f1^-1,1
##  f2 | <id>,1       <id>,2   f3^-1*f1*f3,3
##  f3 |   f3,3       <id>,1          <id>,2
## ----+--------+------------+---------------+
## Adding element: FRElement(...,f3)
## Relator: f3*f2*f1
## gap> ExternalAngle(a);
## {3/8}
## gap> List([1..10],i->ExternalAngle(RotatedSpider(a,i)));
## [ {7/8}, {1/4}, {7/8}, {1/4}, {7/8}, {1/4}, {7/8}, {1/4}, {7/8}, {1/4} ]
## ]]></Example>
##   </Description>
## </ManSection>
## <#/GAPDoc>
##
#############################################################################

#############################################################################
##
#E DBRationalIMGGroup
##
## <#GAPDoc Label="DBRationalIMGGroup">
DeclareGlobalFunction("DBRationalIMGGroup");
## <ManSection>
##   <Func Name="DBRationalIMGGroup" Arg="sequence/map"/>
##   <Returns>An IMG group from Dau's database.</Returns>
##   <Description>
##     This function returns the iterated monodromy group from a database
##     of groups associated to quadratic rational maps. This database has
##     been compiled by Dau Truong Tan <Cite Key="tan:database"/>.
##
##     <P/> When called with no arguments, this command returns the database
##     contents in raw form.
##
##     <P/> The argments can be a sequence; the first integer is the size
##     of the postcritical set, the second argument is an index for the
##     postcritical graph, and sometimes a third argument distinguishes
##     between maps with same post-critical graph.
##
##     <P/> If the argument is a rational map, the command returns the
##     IMG group of that map, assuming its canonical quadratic rational form
##     form exists in the database.
## <Example><![CDATA[
## gap> DBRationalIMGGroup(z^2-1);
## IMG((z-1)^2)
## gap> DBRationalIMGGroup(z^2+1); # not post-critically finite
## fail
## gap> DBRationalIMGGroup(4,1,1);
## IMG((z/h+1)^2|2h^3+2h^2+2h+1=0,h~-0.64)
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareGlobalFunction("PostCriticalMachine");
## <ManSection>
##   <Func Name="PostCriticalMachine" Arg="f"/>
##   <Returns>The Mealy machine of <A>f</A>'s post-critical orbit.</Returns>
##   <Description>
##     This function constructs a Mealy machine <C>P</C> on the alphabet
##     <C>[1]</C>, which describes the post-critical set of <A>f</A>.
##     It is in fact an oriented graph with constant out-degree 1. It is
##     most conveniently passed to <Ref Oper="Draw" BookName="FR"/>.
##
##     <P/> The attribute <C>Correspondence(P)</C> is the list of values
##     associated with the stateset of <C>P</C>.
## <Example><![CDATA[
## gap> z := Indeterminate(Rationals,"z");;
## gap> m := PostCriticalMachine(z^2);
## <Mealy machine on alphabet [ 1 ] with 2 states>
## gap> Display(m);
##    |  1
## ---+-----+
##  a | a,1
##  b | b,1
## ---+-----+
## gap> Correspondence(m);
## [ 0, infinity ]
## gap> m := PostCriticalMachine(z^2-1);; Display(m); Correspondence(m);
##    |  1
## ---+-----+
##  a | c,1
##  b | b,1
##  c | a,1
## ---+-----+
## [ -1, infinity, 0 ]
## ]]></Example>
##   </Description>
## </ManSection>
## <#/GAPDoc>
##
#############################################################################

#############################################################################
##
#E Conversions
##
## <#GAPDoc Label="Conversions">
DeclareAttribute("KneadingSequence", IsRat);
## <ManSection>
##   <Attr Name="KneadingSequence" Arg="angle" Label="angle"/>
##   <Returns>The kneading sequence associated with <A>angle</A>.</Returns>
##   <Description>
##     This function converts a rational angle to a kneading sequence, to
##     describe a quadratic polynomial.
##
##     <P/> If <A>angle</A> is in <M>[1/7,2/7]</M> and the option
##     <C>marked</C> is set, the kneading sequence is decorated with markings
##     in A,B,C.
## <Example><![CDATA[
## gap> KneadingSequence(1/7);
## [ 1, 1 ]
## gap> KneadingSequence(1/5:marked);
## [ "A1", "B1", "B0" ]
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareGlobalFunction("AllInternalAddresses");
## <ManSection>
##   <Attr Name="AllInternalAddresses" Arg="n"/>
##   <Returns>Internal addresses of maps with period up to <A>n</A>.</Returns>
##   <Description>
##     This function returns internal addresses for all periodic points of
##     period up to <A>n</A> under angle doubling. These internal addresses
##     describe the prominent hyperbolic components along the path from the
##     landing point to the main cardioid in the Mandelbrot set; this is a
##     list of length <C>3k</C>, with at position <C>3i+1,3i+2</C> the
##     left and right angles, respectively, and at position <C>3i+3</C> the
##     period of that component. For example,
##     <C>[ 3/7, 4/7, 3, 1/3, 2/3, 2 ]</C> describes the airplane: a 
##     polynomial with landing angles <M>[3/7,4/7]</M> of period <M>3</M>;
##     and such that there is a polynomial with landing angles <M>[1/3,2/3]</M>
##     and period <M>2</M>.
## <Example><![CDATA[
## gap> AllInternalAddresses(3);
## [ [  ], [ [ 1/3, 2/3, 2 ] ], 
## [ [ 1/7, 2/7, 3 ], [ 3/7, 4/7, 3, 1/3, 2/3, 2 ], [ 5/7, 6/7, 3 ] ] ]
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareGlobalFunction("ExternalAnglesRelation");
## <ManSection>
##   <Func Name="ExternalAnglesRelation" Arg="degree, n"/>
##   <Returns>An equivalence relation on the rationals.</Returns>
##   <Description>
##     This function returns the equivalence relation on <C>Rationals</C>
##     identifying all pairs of external angles that land at a
##     common point of period up to <A>n</A> under angle multiplication by
##     by <A>degree</A>.
## <Example><![CDATA[
## gap> ExternalAnglesRelation(2,3);
## <equivalence relation on Rationals >
## gap> EquivalenceRelationPartition(last);
## [ [ 1/7, 2/7 ], [ 1/3, 2/3 ], [ 3/7, 4/7 ], [ 5/7, 6/7 ] ]
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareGlobalFunction("ExternalAngle");
## <ManSection>
##   <Func Name="ExternalAngle" Arg="machine"/>
##   <Returns>The external angle identifying <A>machine</A>.</Returns>
##   <Description>
##     In case <A>machine</A> is the sphere machine of a unicritical
##     polynomial, this function computes the external angle landing at the
##     critical value. More precisely, it computes the equivalence class of
##     that external angle under <Ref Func="ExternalAnglesRelation"/>.
## <Example><![CDATA[
## gap> ExternalAngle(PolynomialSphereMachine(2,[1/7])); # the rabbit
## {2/7}
## gap> Elements(last);
## [ 1/7, 2/7 ]
## ]]></Example>
##   </Description>
## </ManSection>
## <#/GAPDoc>
##
#############################################################################

#E machine.gd . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
