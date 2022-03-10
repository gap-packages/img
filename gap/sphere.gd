#############################################################################
##
#W sphere.gd                                                Laurent Bartholdi
##
#Y Copyright (C) 2013, Laurent Bartholdi
##
#############################################################################
##
##  Sphere groups
##
#############################################################################

DeclareProperty("IsSphereGroup", IsFpGroup);

DeclareCategory("IsElementOfSphereGroup", IsElementOfFpGroup and IsAssocWord);
DeclareCategoryCollections("IsElementOfSphereGroup");
DeclareCategoryCollections("IsElementOfSphereGroupCollection");
DeclareProperty("IsElementOfSphereGroupFamily", IsElementOfFpGroupFamily);
DeclareOperation("ElementOfSphereGroup", [IsFamily, IsAssocWordWithInverse]);

DeclareAttribute("IsomorphismSphereGroup", IsFpGroup);
DeclareAttribute("AsSphereGroup", IsFpGroup);

DeclareAttribute("EulerCharacteristic", IsGroup);

DeclareAttribute("RankOfSphereGroup", IsSphereGroup);

DeclareAttribute("OrderingOfSphereGroup", IsSphereGroup);

DeclareAttribute("ExponentsOfSphereGroup", IsSphereGroup);

DeclareAttribute("IsomorphismFreeGroup", IsSphereGroup);

DeclareGlobalFunction("SphereGroup");

DeclareProperty("IsSphereConjugacyClass", IsAssociativeElementCollection and IsMultiplicativeElementWithInverseCollection);
DeclareProperty("IsSphereConjugacyClassCollection", IsAssociativeElementCollColl and IsMultiplicativeElementWithInverseCollColl);
DeclareSynonym("IsMulticurve", IsSphereConjugacyClassCollection);
DeclareOperation("InverseMutable", [IsSphereConjugacyClass]);
DeclareOperation("InverseSameMutability", [IsSphereConjugacyClass]);
DeclareAttribute("InverseImmutable", IsSphereConjugacyClass);
DeclareOperation("^", [IsSphereConjugacyClass,IsInt]);

DeclareProperty("IsPeripheral", IsElementOfSphereGroup);
DeclareProperty("IsPeripheral", IsSphereConjugacyClass);

DeclareAttribute("PeripheralClasses", IsSphereGroup);

DeclareOperation("IntersectionNumber", [IsSphereConjugacyClass,IsSphereConjugacyClass]);
DeclareOperation("SelfIntersectionNumber", [IsSphereConjugacyClass]);

DeclareProperty("IsAutomorphismGroupOfSphereGroup", IsAutomorphismGroup);
InstallTrueMethod(IsAutomorphismGroup, IsAutomorphismGroupOfSphereGroup);
DeclareAttribute("EpimorphismToOut", IsAutomorphismGroupOfSphereGroup);

DeclareOperation("AmalgamatedFreeProduct", [IsSphereGroup,IsSphereGroup,IsElementOfSphereGroup,IsElementOfSphereGroup]);
DeclareAttribute("EmbeddingsOfAmalgamatedFreeProduct", IsSphereGroup);

## <#GAPDoc Label="SphereGroups">
## <ManSection>
##   <Filt Name="IsSphereGroup"/>
##   <Description>
##     A sphere group is a special kind of finitely presented group, in which
##     exactly one relation is a product, in some order, of all the
##     generators, and all the other relations (possibly none)
##     are powers of generators.
##
##     <P/> Sphere groups are used to represent the fundamental groups of
##     punctured spheres, or more generally orbifolds whose underlying
##     space is a sphere.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Attr Name="IsomorphismSphereGroup" Arg="g"/>
##   <Attr Name="AsSphereGroup" Arg="g"/>
##   <Description>
##     These functions compute an isomorphism from <A>g</A> to a sphere group;
##     the first form returns the isomorphism, while the second one returns
##     its image.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Attr Name="EulerCharacteristic" Arg="g"/>
##   <Returns>The Euler characteristic of <A>g</A>.</Returns>
##   <Description>
##     The Euler characteristic of a free group of rank <M>n</M> is
##     <M>1-n</M>; and it multiplies by the index on subgroups. A
##     sphere group is finite if and only if its Euler characteristic is
##     positive, and is virtually abelian if and only if its Euler
##     characteristic is <M>0</M>.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Attr Name="RankOfSphereGroup" Arg="g"/>
##   <Returns>The number of generators of <A>g</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Attr Name="OrderingOfSphereGroup" Arg="g"/>
##   <Returns>The list of the orders of the generators.</Returns>
##   <Description>
##     This attribute has the property that
##     <C>Product(GeneratorsOfGroup(g){OrderingOfSphereGroup(g)})</C> is the
##     identity.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Attr Name="ExponentsOfSphereGroup" Arg="g"/>
##   <Returns>The list of exponents of the generators.</Returns>
##   <Description>
##     This attribute has the property that
##     <C>GeneratorsOfGroup(g)[i]^ExponentsOfSphereGroup(g)[i]</C> is the
##     identity for all <C>i</C>. If an element has infinite order, the
##     value stored is <C>0</C>.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Attr Name="IsomorphismFreeGroup" Arg="g"/>
##   <Returns>An isomorphism to a free group, if it exists.</Returns>
##   <Description>
##     If <A>g</A> was created as a sphere group with all exponents
##     infinity, then <A>g</A> is isomorphic to a free group on all the
##     generators but one; this attribute stores such an isomorphism.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Func Name="SphereGroup" Arg="ordering, [exponent]"/>
##   <Returns>A new sphere group.</Returns>
##   <Description>
##     <A>ordering</A> is either a list of integers, describing the order
##     of the generators that is to be trivial; or an integer <C>m</C>,
##     in which case the ordering is <C>[m,m-1,..,1]</C>.
##
##     <P/> The optional second argument <A>exponent</A> is a list of
##     integers describing the exponents of the generators. The value
##     <C>0</C> specifies a generator of infinite order.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Filt Name="IsSphereConjugacyClass"/>
##   <Description>
##     Elements of a sphere group represent based loops on a punctured
##     sphere. Loops (without specified basepoint) are represented by
##     conjugacy classes. A <E>multicurve</E> is a collection of
##     non-intersecting loops.
##
##     <P/> Conjugacy classes may be raised to integer powers; the <M>n</M>th
##     power of a conjugacy class is the conjugacy class of the <M>n</M>th
##     power of an element.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Prop Name="IsPeripheral" Arg="c"/>
##   <Returns>Whether the conjugacy class <A>c</A> is peripheral.</Returns>
##   <Description>
##     A conjugacy class is <E>peripheral</E> if it contains a generator of
##     the sphere group.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Attr Name="PeripheralClasses" Arg="g"/>
##   <Returns>The peripheral conjugacy classes of <A>g</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Attr Name="IntersectionNumber" Arg="c, d"/>
##   <Returns>The intersection number of the conjugacy classes <A>c</A> and <A>d</A>.</Returns>
##   <Description>
##     The <E>geometric intersection number</E> of two loops is the minimal
##     number of intersections they may have. The <E>self-intersection
##     number</E> of a loop is the intersection number of the loop with
##     a small translate.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Attr Name="AutomorphismGroup" Arg="g"/>
##   <Attr Name="EpimorphismToOut" Arg="a"/>
##   <Description>
##     This function computes the <E>pure</E> automorphism group of the
##     sphere group <A>g</A>, namely the group of automorphisms that preserves
##     all the peripheral conjugacy classes (conjugacy classes of generators).
##
##     <P/> The attribute <C>EpimorphismToOut</C> stores an epimorphism
##     from the automorphism group to the group of outer automorphisms.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Oper Name="AmalgamateFreeProduct" Arg="g, h, x, y"/>
##   <Description>
##     This function computes the amalgamated free product of two
##     sphere groups <A>g</A> and <A>h</A>, along the cyclic subgroups
##     <M>\langle x\rangle</M> of <A>g</A> and <M>\langle y\rangle</M>
##     of <A>h</A>.
##
##     <P/> The attribute <C>EmbeddingsOfAmalgamatedFreeProduct</C>
##     is a list of length two, storing the embeddings of <A>g</A> and
##     <A>h</A> respectively into the amalgam.
##   </Description>
## </ManSection>
##
## <#/GAPDoc>

#E sphere.gd . . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
