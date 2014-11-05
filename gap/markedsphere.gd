#############################################################################
##
#W markedsphere.gd                                          Laurent Bartholdi
##
#Y Copyright (C) 2013, Laurent Bartholdi
##
#############################################################################
##
##  Marked spheres
##
#############################################################################

## <#GAPDoc Label="MarkedSpheres">
DeclareCategory("IsMarkedSphere", IsObject);
BindGlobal("MARKEDSPHERES_FAMILY",
        NewFamily("MarkedSpheres", IsMarkedSphere));
BindGlobal("TYPE_MARKEDSPHERE",
        NewType(MARKEDSPHERES_FAMILY, IsMarkedSphere));

DeclareAttribute("MarkedSphere", IsSphereMachine);
DeclareAttribute("MarkedSphere", IsP1Map);
# undocumented for now
DeclareAttribute("VerticesOfMarkedSphere", IsMarkedSphere);
DeclareAttribute("SpanningTreeBoundary", IsMarkedSphere);
## <ManSection>
##   <Filt Name="IsMarkedSphere"/>
##   <Description>
##     The category of marked, triangulated spheres
##     (points in Teichmüller space).
##   </Description>
## </ManSection>
##
DeclareOperation("NewMarkedSphere", [IsP1PointCollection,IsSphereGroup]);
DeclareOperation("NewMarkedSphere", [IsP1PointCollection]);
## <ManSection>
##   <Oper Name="NewMarkedSphere" Arg="points [group]"/>
##   <Returns>A new marked sphere on points <A>points</A>.</Returns>
##   <Description>
##     This function creates a new marked sphere, based on the
##     Delaunay triangulation on <A>points</A>. If a sphere group <A>group</A>
##     is specified, it is used to mark the sphere; otherwise a new
##     sphere group is created.
##   </Description>
## </ManSection>
##
DeclareOperation("Draw", [IsMarkedSphere]);
## <ManSection>
##   <Oper Name="Draw" Arg="s" Label="spider"/>
##   <Description>
##     This command plots the marked sphere <A>s</A> in a separate window.
##     It displays the complex sphere, big dots at the post-critical
##     set (feet of the spider), and the arcs and dual arcs
##     of the triangulation connecting the feet.
##
##     <P/> If the option <K>julia:=&lt;gridsize&gt;</K> (if no grid size
##     is specified, it is 500 by default), then the Julia set of the
##     map associated with the spider is also displayed. Points attracted
##     to attracting cycles are coloured in pastel tones, and unattracted
##     points are coloured black.
##
##     <P/> If the option <K>noarcs</K> is specified, the printing of the
##    arcs and dual arcs is disabled.
##
##     <P/> The options <K>upper</K>, <K>lower</K> and <K>detach</K>
##     also apply.
##   </Description>
## </ManSection>
##
DeclareOperation("WiggledMarkedSphere", [IsMarkedSphere,IsObject]);
## <ManSection>
##   <Oper Name="WiggledMarkedSphere" Arg="sphere m"/>
##   <Returns>A new marked sphere.</Returns>
##   <Description>
##     This operation moves the vertices of the marked sphere <A>sphere</A>,
##     preserving its marking. The argument <A>m</A>, which specifies a
##     movement of the vertices, is either a Möbius transformation (to be
##     applied to all vertices) or a list of new positions for them.
##   </Description>
## </ManSection>
##
DeclareOperation("SphereMachineOfBranchedCovering", [IsMarkedSphere,IsMarkedSphere,IsP1Map,IsBool]);
DeclareOperation("SphereMachineOfBranchedCovering", [IsMarkedSphere,IsMarkedSphere,IsP1Map]);
DeclareOperation("SphereMachineAndSphereOfBranchedCovering", [IsMarkedSphere,IsP1Map,IsBool]);
DeclareOperation("SphereMachineAndSphereOfBranchedCovering", [IsMarkedSphere,IsP1Map]);
## <ManSection>
##   <Oper Name="SphereMachineOfBranchedCovering" Arg="down up map [poly]"/>
##   <Oper Name="SphereMachineAndSphereOfBranchedCovering" Arg="down map [poly]"/>
##   <Returns>A sphere machine or [machine,marked sphere].</Returns>
##   <Description>
##     The first function computes, out of a marked sphere <A>down</A> in
##     the range of the P1 map <A>map</A> and a marked sphere <A>up</A> in
##     its domain, the sphere machine representing the monodromy action of
##     the map. Its input stateset is the model group of <A>down</A>, while
##     its output stateset is the model group of <A>up</A>.
##
##     <P/> The second function first computes a marked sphere on the
##     full preimage by <A>map</A> of the vertices of <A>down</A>, then
##     computes the sphere machine, and finally returns a list containing
##     the machine and the sphere at the source of <A>map</A>.
##
##     <P/> The optional parameter <A>poly</A> specifies that the map <A>map</A>
##     is to be treated as a polynomial, and that the machine is to be
##     normalized so that its last generator is an adding machine in
##     standard form.
##   </Description>
## </ManSection>
##
DeclareOperation("MonodromyOfP1Map", [IsMarkedSphere,IsP1Map]);
DeclareOperation("MonodromyOfP1Map", [IsP1PointCollection,IsP1Map]);
DeclareOperation("MonodromyOfP1Map", [IsP1Map]);
## <ManSection>
##   <Oper Name="MonodromyOfP1Map" Arg="[marking] map"/>
##   <Returns>The monodromy action of <A>map</A>.</Returns>
##   <Description>
##     This function computes the monodromy of the P1 map <A>map</A>;
##     this is simply the activity of the sphere machine associated with
##     the map.
##
##     <P/> The optional first argument <A>marking</A> may be a marked sphere,
##     in which case the monodromy is returned as a homomorphism from the
##     marked sphere's marking. It may also be a list of P1 points, in which
##     case the monodromy is returned as a list of permutations, one per
##     point. If the first argument is missing, it is assumed to be the
##     list of critical values of <A>map</A>.
##   </Description>
## </ManSection>
##
DeclareAttribute("SphereMachine", IsP1Map);
## <ManSection>
##   <Oper Name="SphereMachine" Arg="f" Label="rational function"/>
##   <Returns>A sphere machine.</Returns>
##   <Description>
##   This function computes a triangulation of the sphere, on the
##   post-critical set of <A>f</A>, and lifts it through the map <A>f</A>.
##   the action of the fundamental group of the punctured sphere is
##   then read into an sphere machine <C>m</C>, which is returned.
##
##   <P/> This machine has a preset attribute <C>MarkedSphere(m)</C>.
##
##   <P/> An approximation of the Julia set of <A>f</A> can be computed,
##   and plotted on the spider, with the form <C>SphereMachine(f:julia)</C>
##   or <C>SphereMachine(f:julia:=gridsize)</C>.
## <Example><![CDATA[
## gap> SphereMachine(P1z^2-1);
## <FR machine with alphabet [ 1, 2 ] on Group( [ f1, f2, f3 ] )/[ f2*f1*f3 ]>
## gap> Display(last);
##  G  |            1        2
## ----+---------------+--------+
##  f1 |          f2,2   <id>,1
##  f2 | f3^-1*f1*f3,1   <id>,2
##  f3 |        <id>,2     f3,1
## ----+---------------+--------+
## Relator: f2*f1*f3
## ]]></Example>
##   </Description>
## </ManSection>
##
DeclareOperation("DistanceMarkedSpheres", [IsMarkedSphere, IsMarkedSphere]);
DeclareOperation("DistanceMarkedSpheres", [IsMarkedSphere, IsMarkedSphere, IsBool]);
## <ManSection>
##   <Oper Name="DistanceMarkedSpheres" Arg="sphere1 sphere2 [fast]"/>
##   <Returns>The approximate distance between the marked spheres.</Returns>
##   <Description>
##     This function approximates coarsely the Teichmüller distance between
##     marked spheres with same model group.
##     If the vertices of <A>sphere1</A> can be wiggled to
##     the vertices of <A>sphere2</A> in such a manner that the markings
##     coincide, then the distance is the sum of the movements of the
##     vertices. Otherwise, it is <M>1+</M> the sum of the lengths of the
##     images of a sphere group automorphism that carries the marking of
##     <A>sphere1</A> to that of <A>sphere2</A>.
##   </Description>
## </ManSection>
## <#/GAPDoc>

#E markedsphere.gd . . . . . . . . . . . . . . . . . . . . . . . . .ends here
