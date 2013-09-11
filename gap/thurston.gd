#############################################################################
##
#W thurston.gd                                              Laurent Bartholdi
##
#Y Copyright (C) 2011-2013, Laurent Bartholdi
##
#############################################################################
##
##  Thurston's algorithm
##
#############################################################################

## <#GAPDoc Label="Thurston">
DeclareGlobalFunction("NormalizedP1Map");
DeclareProperty("IsBicritical", IsObject);
## <ManSection>
##   <Func Name="NormalizedP1Map" Arg="f M param"/>
##   <Returns>[A canonical conjugate of <A>f</A>,the conjugator].</Returns>
##   <Description>
##     The last argument <A>param</A> is either <C>IsPolynomial</C>,
##     <C>IsBicritical</C> or a positive integer.
##
##     <P/> In the first case, the map <A>f</A> is assumed conjugate to
##     a polynomial. It is conjugated by a Möbius transformation that makes
##     it a <E>centered</E> polynomial, namely a polynomial of the form
##     <M>z^d+a_{d-2}z^{d-2}+\dots+a_0</M>.
##
##     <P/> In the second case, the map <A>f</A> is assumed to have only two
##     critical values; it is normalized as <M>(az^d+b)/(cz^d+e)</M>.
##
##     <P/> In the third case, the map <A>f</A> is assumed to have degree
##     <M>2</M>; it is normalized in the form <M>1+a/z+b/z^2</M>, such that
##     <M>0</M> is on a cycle of length <A>param</A>.
##   </Description>
## </ManSection>
##
DeclareOperation("ThurstonAlgorithm", [IsSphereMachine]);
## <ManSection>
##   <Oper Name="ThurstonAlgorithm" Arg="m"/>
##   <Returns>rec(map := f, machine := M, markedsphere := s).</Returns>
##   <Description>
##     This command runs Thurston's algorithm on the sphere machine
##     <A>m</A>. It either returns a record with the P1 map <C>f</C> to
##     which the algorithm converged, as well as the marked sphere with
##     <C>f</C>'s post-critical set and a simplified machine equivalent to
##     <A>m</A>; or a record returned by <Ref Oper="ThurstonObstruction"/>.
##   </Description>
## </ManSection>
##
DeclareOperation("P1MapBySphereMachine", [IsSphereMachine]);
## <ManSection>
##   <Oper Name="P1MapBySphereMachine" Arg="m"/>
##   <Returns>Either a map or an obstruction.</Returns>
##   <Description>
##     This command returs either the map computed by
##     <Ref Oper="ThurstonAlgorithm"/> or a Thurston obstruction.
##
##   <P/> It runs a modification of Hubbard and Schleicher's
##   "spider algorithm" <Cite Key="MR1315537"/> on the sphere machine <A>m</A>.
##
##   <P/> The command accepts the following options, to return a map in a given normalization: <List>
##   <Mark><C>P1MapBySphereMachine(m:param:=IsPolynomial)</C></Mark>
##         <Item>returns <M>f=z^d+A_{d-2}z^{d-2}+\cdots+A_0</M>;</Item>
##   <Mark><C>P1MapBySphereMachine(m:param:=IsBicritical)</C></Mark>
##         <Item>returns <M>f=((pz+q)/(rz+s)^d</M>, with
##               <M>1</M>postcritical;</Item>
##   <Mark><C>P1MapBySphereMachine(m:param:=n)</C></Mark>
##         <Item>returns <M>f=1+a/z+b/z^2</M> or <M>f=a/(z^2+2z)</M>
##               if <C>n=2</C>.</Item>
##   </List>
## <Example><![CDATA[
## gap> m := PolynomialSphereMachine(2,[1/3],[]);
## <FR machine with alphabet [ 1, 2 ] on Group( [ f1, f2, f3 ] )/[ f3*f2*f1 ]>
## gap> P1MapBySphereMachine(m);
## 0.866025*z^2+(-1)*z+(-0.288675)
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareOperation("ThurstonMatrix", [IsSphereMachine,IsMulticurve]);
## <ManSection>
##   <Oper Name="ThurstonMatrix" Arg="m multicurve"/>
##   <Returns>The transition matrix of the multicurve.</Returns>
##   <Description>
##     This command computes the iterated preimages of the multicurve
##     <A>multicurve</A> till it obtains a backwards-invariant multicurve
##     or some preimages intersect. In the latter case, <K>fail</K> returned,
##     while in the former case the Thurston matrix of the multicurve is
##     returned.
## <Example><![CDATA[
## gap> r := PolynomialSphereMachine(2,[],[1/6]);;
## gap> F := StateSet(r);;
## gap> twist := GroupHomomorphismByImages(F,F,GeneratorsOfGroup(F),[F.1,F.2^(F.3*F.2),F.3^F.2,F.4]);;
## gap> SupportingRays(r*twist^-1);
## rec( machine := <FR machine with alphabet [ 1, 2 ] on F/[ f4*f1*f2*f3 ]>,
##      twist := [ f1, f2, f3, f4 ] -> [ f1, f3^-1*f2*f3, f3^-1*f2^-1*f3*f2*f3, f4 ],
##      obstruction := "Dehn twist" )
## gap> ThurstonMatrix(last.machine,[ConjugacyClass(F.2*F.3)]);
## [ [ 1 ] ]
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareOperation("ThurstonObstruction", [IsSphereMachine,IsMarkedSphere]);
## <ManSection>
##   <Oper Name="ThurstonObstruction" Arg="m sphere"/>
##   <Returns>rec(multicurve := mc, matrix := m) or <K>fail</K>.</Returns>
##   <Description>
##     This command tries to find a Thurston obstruction (multicurve such
##     that its Thurston matrix has spectral radius at least <M>1</M>);
##     it either returns <K>fail</K> if the search was inconclusive, or
##     a record describing the obstruction.
##
##     <P/> The obstruction is searched for by considering small subtrees of
##     the minimal spanning tree of <A>sphere</A>, computing loops surrounding
##     these subtrees, and saturating them into a multicurve by taking their
##     iterated preimages, see <Ref Oper="ThurstonMatrix"/>.
##   </Description>
## </ManSection>
## <#/GAPDoc>

#E thurston.gd . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
