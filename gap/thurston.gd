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
## <ManSection>
##   <Func Name="NormalizedQuadraticP1Map" Arg="f M param"/>
DeclareGlobalFunction("NormalizedQuadraticP1Map");
DeclareProperty("IsBicritical", IsObject);
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
## <ManSection>
##   <Oper Name="ThurstonAlgorithm" Arg="m"/>
DeclareOperation("ThurstonAlgorithm", [IsSphereMachine]);
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
## <ManSection>
##   <Oper Name="P1MapBySphereMachine" Arg="m"/>
DeclareOperation("P1MapBySphereMachine", [IsSphereMachine]);
##   <Returns>Either a map or an obstruction.</Returns>
##   <Description>
##     This command returs either the map computed by
##     <Ref Oper="ThurstonAlgorithm"/> or a Thurston obstruction.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Oper Name="ThurstonMatrix" Arg="m multicurve"/>
DeclareOperation("ThurstonMatrix", [IsSphereMachine,IsMulticurve]);
##   <Returns>The transition matrix of the multicurve.</Returns>
##   <Description>
##     This command computes the iterated preimages of the multicurve
##     <A>multicurve</A> till it obtains a backwards-invariant multicurve
##     or some preimages intersect. In the latter case, <K>fail</K> returned,
##     while in the former case the Thurston matrix of the multicurve is
##     returned.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Oper Name="ThurstonObstruction" Arg="m sphere"/>
DeclareOperation("ThurstonObstruction", [IsSphereMachine,IsMarkedSphere]);
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
##
## <#/GAPDoc>

#E thurston.gd . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
