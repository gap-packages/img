#############################################################################
##
#W triangulations.gd                                        Laurent Bartholdi
##
#Y Copyright (C) 2013, Laurent Bartholdi
##
#############################################################################
##
##  Declarations for triangulations
##
#############################################################################

#############################################################################
##
#E Triangulations
##
## <#GAPDoc Label="Triangulations">
## <ManSection>
##   <Oper Name="DelaunayTriangulation" Arg="points, [quality]"/>
##   <Returns>A Delaunay triangulation of the sphere.</Returns>
##   <Description>
##     If <A>points</A> is a list of points on the unit sphere, represented
##     by their 3D coordinates, this function creates a triangulation of
##     the sphere with these points as vertices. This triangulation is
##     such that the angles are as equilateral as possible.
##
##     <P/> This triangulation is a recursive collection of records, one
##     for each vertex, oriented edge or face. Each such object has a
##     <C>pos</C> component giving its coordinates; and an <C>index</C>
##     component identifying it uniquely. Additionally, vertices and
##     faces have a <C>n</C> component which lists their neighbours in CCW
##     order, and edges have <C>from,to,left,right,reverse</C> components.
##
##     <P/> If all points are aligned on a great circle, or if all points
##     are in a hemisphere, some points are added so as to make the
##     triangulation simplicial with all edges of length <M>&lt;\pi</M>.
##     These vertices additionally have a <C>fake</C> component set to
##     <K>true</K>.
##
##     <P/> A triangulation may be plotted with <C>Draw</C>; this requires
##     <Package>appletviewer</Package> to be installed. The command
##     <C>Draw(t:detach)</C> detaches the subprocess after it is started.
##     The extra arguments <C>Draw(t:lower)</C> or <C>Draw(t:upper)</C>
##     stretch the triangulation to the lower, respectively upper, hemisphere.
##
##     <P/> If the second argument <A>quality</A>, which must be a floatean,
##     is present, then all triangles in the resulting triangulation are
##     guaranteed to have circumcircle ratio / minimal edge length at most
##     <A>quality</A>. Of course, additional vertices may need to be added
##     to ensure that.
## <Example><![CDATA[
## gap> octagon := Concatenation(IdentityMat(3),-IdentityMat(3))*1.0;;
## gap> dt := DelaunayTriangulation(octagon);
## <triangulation with 6 vertices, 24 edges and 8 faces>
## gap> dt!.v;
## [ <vertex 1>, <vertex 2>, <vertex 3>, <vertex 4>, <vertex 5>, <vertex 6> ]
## gap> last[1].n;
## [ <edge 17>, <edge 1>, <edge 2>, <edge 11> ]
## gap> last[1].from;
## <vertex 1>
## ]]></Example>
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Oper Name="LocateInTriangulation" Arg="t,[seed,]point"/>
##   <Returns>The face in <A>t</A> containing <A>point</A>.</Returns>
##   <Description>
##     This command locates the face in <A>t</A> that contains <A>point</A>;
##     or, if <A>point</A> lies on an edge or a vertex, it returns that
##     edge or vertex.
##
##     <P/> The optional second argument specifies a starting vertex,
##     edge, face, or vertex index from which to start the search. Its only
##     effect is to speed up the algorithm.
## <Example><![CDATA[
## gap> cube := Tuples([-1,1],3)/Sqrt(3.0);;
## gap> dt := DelaunayTriangulation(cube);
## <triangulation with 8 vertices, 36 edges and 12 faces>
## gap> LocateInTriangulation(dt,dt!.v[1].pos);
## <vertex 1>
## gap> LocateInTriangulation(dt,[3/5,0,4/5]*1.0);
## <face 9>
## ]]></Example>
##   </Description>
## </ManSection>
## <#/GAPDoc>
##
DeclareCategory("IsSphereTriangulation", IsObject);
BindGlobal("TRIANGULATION_FAMILY",
        NewFamily("SphereTriangulations", IsSphereTriangulation));
BindGlobal("TYPE_TRIANGULATION",
        NewType(TRIANGULATION_FAMILY, IsSphereTriangulation));

DeclareRepresentation("IsTriangulationObjectRep",
        IsComponentObjectRep and IsAttributeStoringRep,[]);
DeclareCategory("IsTriangulationObject",IsTriangulationObjectRep);
DeclareCategory("IsVertex",IsTriangulationObject);
DeclareCategory("IsEdge",IsTriangulationObject);
DeclareCategory("IsFace",IsTriangulationObject);
BindGlobal("TRIANGULATIONOBJECT_FAMILY",
        NewFamily("TriangulationFamily",IsTriangulationObject,CanEasilySortElements,CanEasilySortElements));
BindGlobal("TYPE_VERTEX",
        NewType(TRIANGULATIONOBJECT_FAMILY,IsVertex));
BindGlobal("TYPE_EDGE",
        NewType(TRIANGULATIONOBJECT_FAMILY,IsEdge));
BindGlobal("TYPE_FACE",
        NewType(TRIANGULATIONOBJECT_FAMILY,IsFace));

#############################################################################

DeclareOperation("DelaunayTriangulation", [IsList]);
DeclareOperation("DelaunayTriangulation", [IsList, IsFloat]);
DeclareOperation("AddToTriangulation", [IsSphereTriangulation,IsP1Point]);
DeclareOperation("AddToTriangulation", [IsSphereTriangulation,IsP1Point,IsBool]);
DeclareOperation("AddToTriangulation", [IsSphereTriangulation,IsFace,IsP1Point]);
DeclareOperation("AddToTriangulation", [IsSphereTriangulation,IsFace,IsP1Point,IsBool]);
DeclareOperation("RemoveFromTriangulation", [IsSphereTriangulation,IsVertex]);
DeclareOperation("LocateInTriangulation", [IsSphereTriangulation,IsP1Point]);
DeclareOperation("LocateInTriangulation", [IsSphereTriangulation,IsObject,IsP1Point]);
DeclareOperation("Draw", [IsSphereTriangulation]);

DeclareAttribute("Neighbour", IsVertex);
DeclareOperation("Neighbours", [IsVertex]);
DeclareOperation("Neighbours", [IsVertex,IsEdge]);
DeclareAttribute("Pos", IsVertex);
DeclareProperty("IsFake", IsVertex);
DeclareOperation("Valency", [IsVertex]);

DeclareAttribute("Left", IsEdge);
DeclareAttribute("Right", IsEdge);
DeclareAttribute("To", IsEdge);
DeclareAttribute("From", IsEdge);
DeclareAttribute("Next", IsEdge);
DeclareAttribute("Prevopp", IsEdge);
DeclareAttribute("Opposite", IsEdge);
DeclareAttribute("Pos", IsEdge);
DeclareAttribute("FromPos", IsEdge);
DeclareAttribute("ToPos", IsEdge);
DeclareAttribute("Length", IsEdge);
DeclareAttribute("Map", IsEdge);
DeclareAttribute("GroupElement", IsEdge);

DeclareAttribute("Neighbour", IsFace);
DeclareOperation("Neighbours", [IsFace]);
DeclareOperation("Neighbours", [IsFace,IsEdge]);
DeclareAttribute("Pos", IsFace);
DeclareAttribute("Radius", IsFace);
DeclareAttribute("Centre", IsFace);
DeclareOperation("Valency", [IsFace]);

DeclareOperation("ClosestFaces", [IsTriangulationObject]);
DeclareOperation("ClosestVertices", [IsTriangulationObject]);

#E triangulations.gd . . . . . . . . . . . . . . . . . . . . . . . .ends here
