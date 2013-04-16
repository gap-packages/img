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

#############################################################################
##
#E MarkedSpheres
##
## <#GAPDoc Label="MarkedSpheres">
## <ManSection>
##   <Filt Name="IsSphereTriangulation"/>
##   <Filt Name="IsMarkedSphere"/>
##   <Attr Name="Spider" Arg="ratmap" Label="r"/>
##   <Attr Name="Spider" Arg="machine" Label="m"/>
##   <Description>
##     The category of triangulated spheres (points in Moduli space),
##     or of marked, triangulated spheres (points in Teichmüller space).
##
##     <P/> Various commands have an attribudte <C>Spider</C>, which records
##     this point in Teichmüller space.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Oper Name="RationalFunction" Arg="[z,] m"/>
##   <Returns>A rational function.</Returns>
##   <Description>
##   This command runs a modification of Hubbard and Schleicher's
##   "spider algorithm" <Cite Key="MR1315537"/> on the IMG FR machine <A>m</A>.
##   It either returns a rational function <C>f</C> whose associated machine
##   is <A>m</A>; or a record describing the Thurston obstruction to
##   realizability of <C>f</C>.
##
##   <P/> This obstruction record <C>r</C> contains a list <C>r.multicurve</C>
##   of conjugacy classes in <C>StateSet(m)</C>, which represent
##   short multicurves; a matrix <C>r.mat</C>, and a spider <C>r.spider</C>
##   on which the obstruction was discovered.
##
##   <P/> If a rational function is returned, it has preset attributes
##   <C>Spider(f)</C> and <C>IMGMachine(f)</C> which is a simplified
##   version of <A>m</A>. This rational function is also normalized so that
##   its post-critical points have barycenter=0 and has two post-critical
##   points at infinity and on the positive real axis.
##   Furthermore, if <A>m</A> is polynomial-like, then the returned map is
##   a polynomial.
##
##   <P/> The command accepts the following options, to return a map in a given normalization: <List>
##   <Mark><C>RationalFunction(m:param:=IsPolynomial)</C></Mark>
##         <Item>returns <M>f=z^d+A_{d-2}z^{d-2}+\cdots+A_0</M>;</Item>
##   <Mark><C>RationalFunction(m:param:=IsBicritical)</C></Mark>
##         <Item>returns <M>f=((pz+q)/(rz+s)^d</M>, with
##               <M>1</M>postcritical;</Item>
##   <Mark><C>RationalFunction(m:param:=n)</C></Mark>
##         <Item>returns <M>f=1+a/z+b/z^2</M> or <M>f=a/(z^2+2z)</M>
##               if <C>n=2</C>.</Item>
##   </List>
## <Example><![CDATA[
## gap> m := PolynomialIMGMachine(2,[1/3],[]);
## <FR machine with alphabet [ 1, 2 ] on Group( [ f1, f2, f3 ] )/[ f3*f2*f1 ]>
## gap> RationalFunction(m);
## 0.866025*z^2+(-1)*z+(-0.288675)
## ]]></Example>
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Oper Name="Draw" Arg="s" Label="spider"/>
##   <Description>
##     This command plots the spider <A>s</A> in a separate X window.
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

## <ManSection>
##   <Oper Name="SphereMachine" Arg="f" Label="rational function"/>
##   <Returns>A sphere machine.</Returns>
##   <Description>
##   This function computes a triangulation of the sphere, on the
##   post-critical set of <A>f</A>, and lifts it through the map <A>f</A>.
##   the action of the fundamental group of the punctured sphere is
##   then read into an IMG fr machine <C>m</C>, which is returned.
##
##   <P/> This machine has a preset attribute <C>Spider(m)</C>.
##
##   <P/> An approximation of the Julia set of <A>f</A> can be computed,
##   and plotted on the spider, with the form <C>IMGMachine(f:julia)</C>
##   or <C>IMGMachine(f:julia:=gridsize)</C>.
## <Example><![CDATA[
## gap> z := Indeterminate(COMPLEX_FIELD);;
## gap> IMGMachine(z^2-1);
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
## <ManSection>
##   <Oper Name="FindThurstonObstruction" Arg="list"/>
##   <Returns>A description of the obstruction corresponding to <A>list</A>, or <K>fail</K>.</Returns>
##   <Description>
##     This method accepts a list of IMG elements on the same underlying
##     machine, and treats these as representatives of conjugacy classes
##     defining (part of) a multicurve. It computes whether these
##     curves, when supplemented with their lifts under the recursion,
##     constitute a Thurston obstruction, by computing its transition matrix.
##
##     <P/> The method either returns <K>fail</K>, if there is no obstruction,
##     or a record with as fields <C>matrix,machine,obstruction</C> giving
##     respectively the transition matrix, a simplified machine, and the
##     curves that constitute a minimal obstruction.
## <Example><![CDATA[
## gap> r := PolynomialIMGMachine(2,[],[1/6]);;
## gap> F := StateSet(r);;
## gap> twist := GroupHomomorphismByImages(F,F,GeneratorsOfGroup(F),[F.1,F.2^(F.3*F.2),F.3^F.2,F.4]);;
## gap> SupportingRays(r*twist^-1);
## rec( machine := <FR machine with alphabet [ 1, 2 ] on F/[ f4*f1*f2*f3 ]>,
##      twist := [ f1, f2, f3, f4 ] -> [ f1, f3^-1*f2*f3, f3^-1*f2^-1*f3*f2*f3, f4 ],
##      obstruction := "Dehn twist" )
## gap> FindThurstonObstruction([FRElement(last.machine,[2,3])]);
## rec( matrix := [ [ 1 ] ], machine := <FR machine with alphabet [ 1, 2 ] on F/[ f4*f1*f2*f3 ]>, obstruction := [ f1^-1*f4^-1 ] )
## ]]></Example>
##   </Description>
## </ManSection>
## <#/GAPDoc>
##
DeclareCategory("IsMarkedSphere", IsObject);
BindGlobal("MARKEDSPHERES_FAMILY",
        NewFamily("MarkedSpheres", IsMarkedSphere));
BindGlobal("TYPE_MARKEDSPHERE",
        NewType(MARKEDSPHERES_FAMILY, IsMarkedSphere));

DeclareOperation("Draw", [IsMarkedSphere]);
DeclareAttribute("Vertices", IsMarkedSphere);
DeclareAttribute("SpanningTreeBoundary", IsMarkedSphere);
DeclareOperation("NewMarkedSphere", [IsP1PointCollection,IsSphereGroup]);
DeclareOperation("NewMarkedSphere", [IsP1PointCollection]);
DeclareOperation("WiggledMarkedSphere", [IsMarkedSphere,IsObject]);

DeclareAttribute("MarkedSphere", IsSphereMachine);
DeclareAttribute("MarkedSphere", IsP1Map);

DeclareOperation("SphereMachineOfBranchedCovering", [IsMarkedSphere,IsMarkedSphere,IsP1Map,IsBool]);
DeclareOperation("SphereMachineOfBranchedCovering", [IsMarkedSphere,IsMarkedSphere,IsP1Map]);
DeclareOperation("SphereMachineAndSphereOfBranchedCovering", [IsMarkedSphere,IsP1Map,IsBool]);
DeclareOperation("SphereMachineAndSphereOfBranchedCovering", [IsMarkedSphere,IsP1Map]);

DeclareOperation("MonodromyOfP1Map", [IsMarkedSphere,IsP1Map]);
DeclareOperation("MonodromyOfP1Map", [IsP1PointCollection,IsP1Map]);
DeclareOperation("MonodromyOfP1Map", [IsP1Map]);

DeclareAttribute("SphereMachine", IsP1Map);

DeclareOperation("DistanceMarkedSpheres", [IsMarkedSphere, IsMarkedSphere]);
DeclareOperation("DistanceMarkedSpheres", [IsMarkedSphere, IsMarkedSphere, IsBool]);
#############################################################################

#E markedsphere.gd . . . . . . . . . . . . . . . . . . . . . . . . .ends here
