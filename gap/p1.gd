#############################################################################
##
#W p1.gd                                                    Laurent Bartholdi
##
#Y Copyright (C) 2012-2013, Laurent Bartholdi
##
#############################################################################
##
##  This file declares code for P1 points
##
#############################################################################

## <#GAPDoc Label="P1Points">
## <ManSection>
##   <Filt Name="IsP1Point"/>
DeclareCategory("IsP1Point",IsObject);
DeclareCategoryCollections("IsP1Point");
DeclareCategoryCollections("IsP1PointCollection");
DeclareSynonym("IsP1PointList",IsP1PointCollection and IsList);
DeclareCategory("IsIEEE754P1Point",IsP1Point);
##   <Fam Name="P1PointsFamily"/>
BindGlobal("P1PointsFamily",NewFamily("P1PointsFamily",IsP1Point));
BindGlobal("TYPE_P1POINT",NewType(P1PointsFamily,IsP1Point and IsPositionalObjectRep));
BindGlobal("TYPE_IEEE754P1POINT",NewType(P1PointsFamily,IsIEEE754P1Point and IsDataObjectRep));
##   <Func Name="P1Point" Arg="complex"/>
##   <Func Name="P1Point" Arg="real, imag" Label="ri"/>
##   <Func Name="P1Point" Arg="string" Label="s"/>
DeclareOperation("P1Point",[IsFloat]);
DeclareOperation("P1Point",[IsRat]);
DeclareOperation("P1Point",[IsInfinity]);
DeclareOperation("P1Point",[IsFloat,IsFloat]);
##   <Description>
##     P1 points are complex numbers or infinity;
##     fast methods are implemented to compute with them, and to apply
##     rational maps to them.
##     <P/>
##     The first filter recognizes these objects. Next, the family they
##     belong to. The next methods create a new P1 point.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Coordinate" Arg="p"/>
DeclareAttribute("P1Coordinate",IsP1Point);
##   <Returns>The complex number represented by <A>p</A>.</Returns>
## </ManSection>
##

## <ManSection>
##   <Func Name="CleanedP1Point" Arg="p, [prec]"/>
DeclareOperation("CleanedP1Point",[IsP1Point,IsFloat]);
DeclareOperation("CleanedP1Point",[IsP1Point]);
##   <Returns><A>p</A>, rounded towards 0/1/infinity/reals at precision <A>prec</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Var Name="P1infinity"/>
DeclareGlobalVariable("P1infinity");
DeclareOperation("P1INFINITY@",[IsP1Point]);
##   <Var Name="P1one"/>
DeclareGlobalVariable("P1one");
##   <Var Name="P1zero"/>
DeclareGlobalVariable("P1zero");
##   <Description>The south, north and 'east' poles of the Riemann sphere.</Description>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Antipode" Arg="p"/>
DeclareAttribute("P1Antipode",IsP1Point);
##   <Returns>The antipode of <A>p</A> on the Riemann sphere.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Barycentre" Arg="points ..."/>
DeclareOperation("P1Barycentre",[IsP1PointList]);
DeclareOperation("P1Barycentre",[IsP1Point]);
DeclareOperation("P1Barycentre",[IsP1Point,IsP1Point]);
DeclareOperation("P1Barycentre",[IsP1Point,IsP1Point,IsP1Point]);
##   <Returns>The barycentre of its arguments (which can also be a list of P1 points).</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Circumcentre" Arg="p, q, r"/>
DeclareOperation("P1Circumcentre",[IsP1Point,IsP1Point,IsP1Point]);
##   <Returns>The centre of the smallest disk containing <A>p,q,r</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Distance" Arg="p, q"/>
DeclareOperation("P1Distance",[IsP1Point,IsP1Point]);
##   <Returns>The spherical distance from <A>p</A> to <A>q</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Midpoint" Arg="p, q"/>
DeclareOperation("P1Midpoint",[IsP1Point,IsP1Point]);
##   <Returns>The point between <A>p</A> to <A>q</A> (undefined if they are antipodes of each other).</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Sphere" Arg="v"/>
DeclareAttribute("P1Sphere",IsList);
##   <Returns>The P1 point corresponding to <A>v</A> in <M>\mathbb R^3</M>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="SphereP1" Arg="p"/>
DeclareAttribute("SphereP1",IsP1Point);
##   <Returns>The coordinates in <M>\mathbb R^3</M> of <A>p</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="SphereP1Y" Arg="p"/>
DeclareAttribute("SphereP1Y",IsP1Point);
##   <Returns>The Y coordinate in <M>\mathbb R^3</M> of <A>p</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1XRatio" Arg="p q r s"/>
##   <Func Name="XRatio" Arg="p q r s"/>
DeclareOperation("P1XRatio",[IsP1Point,IsP1Point,IsP1Point,IsP1Point]);
##   <Returns>The cross ratio of <A>p, q, r, s</A>.</Returns>
##   <Description>
##     The cross ratio of four points <A>p,q,r,s</A> is defined as
##     <C>(p-r)(q-s)/(p-s)(q-r)</C>. The values <C>P1zero,P1one,P1infinity</C>
##     correspond respectively to the special cases
##     <C>(p=r or q=s)</C>, <C>(p=q or r=s)</C>, <C>(p=s or q=r)</C>.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Func Name="XRatio" Arg="p q r s"/>
DeclareOperation("XRatio",[IsP1Point,IsP1Point,IsP1Point,IsP1Point]);
##   <Returns>The cross ratio of <A>p, q, r, s</A>, as a complex number.</Returns>
##   <Description>
##     The cross ratio of four points <A>p,q,r,s</A> is defined as <C>(p-r)(q-s)/(p-s)(q-r)</C>.
##     The values <C>0,1,infinity</C> correspond respectively to the special cases
##     <C>(p=r or q=s)</C>, <C>(p=q or r=s)</C>, <C>(p=s or q=r)</C>.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Oper Name="CollectedP1Points" Arg="p1points [,precision]"/>
DeclareOperation("CollectedP1Points",[IsP1PointList]);
DeclareOperation("CollectedP1Points",[IsP1PointList,IsFloat]);
##   <Returns>A list of pairs <C>[point,multiplicity]</C>.</Returns>
##   <Description>
##     Collects the points in <A>p1points</A>; points at distance at most <A>precision</A>
##     are considered equal, and the barycentre of the clustered points is returned.
##     <P/>
##     If the argument <A>precision</A> is not supplied, <C>@IMG.p1eps</C> is taken.
## </ManSection>
##
## <ManSection>
##   <Oper Name="MatchP1Points" Arg="p1pointsA, p1pointsB [,separation]"/>
DeclareOperation("MatchP1Points",[IsP1PointList,IsP1PointCollColl,IsFloat]);
DeclareOperation("MatchP1Points",[IsP1PointList,IsP1PointCollColl]);
##   <Returns>A list of giving closest point in <A>p1pointsB</A> to points in <A>p1pointsA</A>, or <K>fail</K>.</Returns>
##   <Description>
##     Finds for each point <A>p1pointsA[i]</A> the closest point <A>p1pointB[p[i]]</A>.
##     If the next-closest is at least <A>separation</A> further away for all <C>i</C>,
##     then the list <C>p</C> is returned. Otherwise, <K>fail</K> is returned.
##     <P/>
##     If the argument <A>separation</A> is not supplied, <K>2</K> is taken.
## </ManSection>
##
## <ManSection>
##   <Filt Name="IsP1Map"/>
DeclareSynonym("IsP1Map",IsUnivariateRationalFunction and IsFloatRationalFunction);
DeclareCategory("IsIEEE754P1Map",IsP1Map);
BindGlobal("TYPE_IEEE754P1MAP", NewType(RationalFunctionsFamily(PMCOMPLEX_PSEUDOFIELD), IsIEEE754P1Map and IsDataObjectRep));
##   <Description>
##     P1 maps are stored more efficiently than rational functions, but are
##     otherwise equivalent.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Func Name="MoebiusMap" Arg="[sourcelist], destlist"/>
##   <Func Name="MoebiusMap" Arg="p, q, r, s, t, u" Label="6"/>
##   <Func Name="MoebiusMap" Arg="p, q, r" Label="3"/>
##   <Func Name="MoebiusMap" Arg="p, q" Label="2"/>
##   <Func Name="MoebiusMap" Arg="p" Label="1"/>
DeclareOperation("MoebiusMap",[IsP1Point]);
DeclareOperation("MoebiusMap",[IsP1Point,IsP1Point]);
DeclareOperation("MoebiusMap",[IsP1Point,IsP1Point,IsP1Point]);
DeclareOperation("MoebiusMap",[IsP1Point,IsP1Point,IsP1Point,IsP1Point,IsP1Point,IsP1Point]);
DeclareOperation("MoebiusMap",[IsP1PointList]);
DeclareOperation("MoebiusMap",[IsP1PointList,IsP1PointList]);
##   <Returns>A new Möbius transformation.</Returns>
##   <Description>
##     In the first case, this is the Möbius transformation sending <A>p,q,r</A> to <A>P,Q,R</A>
##     respectively;
##     in the second case, the map sending <C>0,1,P1infinity</C> to  <A>p,q,r</A> respectively;
##     in the third case, the map sending <C>0,P1infinity</C> to <A>p,q</A> respectively,
##     and of the form <M>(z-p)/(z-q)</M>;
##     and in the fourth case, a rotation sending <A>P1infinity</A> to <A>p</A>.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1MapRotatingP1Points" Arg="points [oldpoints]"/>
DeclareOperation("P1ROTATION@",[IsP1Point,IsP1PointList,IsP1PointList]);
DeclareGlobalFunction("P1MapRotatingP1Points");
##   <Returns>A Möbius rotation sending the last of <A>points</A> to P1infinity.</Returns>
##   <Description>
##     A Möbius rotation is computed that sends the last of <A>points</A> to <K>P1infinity</K> and,
##     assuming the last point in <A>oldpoints</A> is also <K>P1infinity</K>, matches the points
##     in <A>points</A> and <A>oldpoints</A> as closely as possible.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1MapNormalizingP1Points" Arg="points [oldpoints]"/>
DeclareOperation("P1MapNormalizingP1Points",[IsP1PointList]);
DeclareOperation("P1MapNormalizingP1Points",[IsP1PointList,IsP1PointList]);
##   <Returns>A Möbius transformation sending the last of <A>points</A> to P1infinity.</Returns>
##   <Description>
##     A Möbius transformation is computed that sends the last of <A>points</A> to
##     <K>P1infinity</K> and makes the barycentre of the points in <M>\mathbb R^3</M>
##     as close as possible to the origin.
##     If a list <A>oldpoints</A> is also given, the Möbius transformation computed rotates
##     about infinity so as to match the points in <A>points</A> and <A>oldpoints</A> as
##     closely as possible; this then determines the transformation uniquely.
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Var Name="P1z"/>
DeclareGlobalVariable("P1z");
##   <Description>The identity Möbius transformation.</Description>
## </ManSection>
##
## <ManSection>
##   <Func Name="P1Monomial" Arg="n"/>
DeclareGlobalFunction("P1Monomial");
##   <Returns>The rational function <M>z^n</M>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="CleanedP1Map" Arg="map, [prec]"/>
DeclareOperation("CleanedP1Map",[IsP1Map,IsFloat]);
DeclareOperation("CleanedP1Map",[IsP1Map]);
##   <Returns><A>map</A>, with coefficients rounded using <A>prec</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="CompositionP1Map" Arg="map1, ..."/>
DeclareSynonym("CompositionP1Map",CompositionMapping2);
##   <Returns>The composition of the maps passed as arguments, in the functional (<A>map1</A> last) order.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="InverseP1Map" Arg="map"/>
DeclareSynonym("InverseP1Map",InverseGeneralMapping);
##   <Returns>The functional inverse of the Möbius transformation <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="ComplexConjugate" Arg="map"/>
#DeclareAttribute("ComplexConjugate",IsP1Map);
##   <Returns>The complex conjugated map.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="ConjugatedP1Map" Arg="map, mobius"/>
DeclareOperation("ConjugatedP1Map",[IsP1Map,IsP1Map]);
##   <Returns>The map <C>CompositionP1Map(InverseP1Map(mobius),map,mobius)</C>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="CoefficientsOfP1Map" Arg="map"/>
DeclareAttribute("CoefficientsOfP1Map",IsP1Map);
##   <Returns>Coefficients of numerator and denominator of <A>map</A>, lowest degree first.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1MapByCoefficients" Arg="numer, denom"/>
DeclareGlobalFunction("P1MapByCoefficients");
DeclareOperation("P1MAPBYCOEFFICIENTS2@",[IsObject,IsList,IsList]);
##   <Returns>The P1 map with numerator coefficients <A>numer</A> and denominator <A>denom</A>, lowest degree first.</Returns>
## </ManSection>
##
## <ManSection>
##   <Attr Name="NumeratorP1Map" Arg="map"/>
##   <Attr Name="DenominatorP1Map" Arg="map"/>
DeclareAttribute("NumeratorP1Map",IsP1Map);
DeclareAttribute("DenominatorP1Map",IsP1Map);
##   <Returns>The numerator/denominator of <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1MapByZerosPoles" Arg="zeros, poles, src, dest"/>
DeclareOperation("P1MapByZerosPoles",[IsP1PointList,IsP1PointList,IsP1Point,IsP1Point]);
##   <Returns>The P1 map with specified zeros and poles, and sending <A>src</A> to <A>dest</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1Path" Arg="p q"/>
DeclareOperation("P1Path",[IsP1Point,IsP1Point]);
##   <Returns>The P1 map sending <C>0</C> to <A>p</A> and <C>1</C> to <A>q</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="DegreeOfP1Map" Arg="map"/>
DeclareAttribute("DegreeOfP1Map",IsP1Map);
##   <Returns>The degree of <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1Image" Arg="map, p1point"/>
DeclareSynonym("P1Image",ImageElm);
DeclareOperation("ImageElm",[IsP1Map,IsP1Point]);
##   <Returns>The image of <A>p1point</A> under <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1PreImage" Arg="map, p1point"/>
DeclareSynonym("P1PreImage",PreImageElm);
DeclareOperation("PreImageElm",[IsP1Map,IsP1Point]);
##   <Returns>The preimage of <A>p1point</A> under <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1PreImages" Arg="map, p1point"/>
DeclareSynonym("P1PreImages",PreImagesElm);
DeclareOperation("PreImagesElm",[IsP1Map,IsP1Point]);
##   <Returns>The preimages of <A>p1point</A> under <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="Primitive" Arg="map"/>
##   <Oper Name="Derivative" Arg="map"/>
DeclareAttribute("Primitive",IsP1Map);
#DeclareAttribute("Derivative",IsP1Map);
##   <Returns>The [anti]derivative of the rational map <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Attr Name="CriticalPointsOfP1Map" Arg="map"/>
DeclareAttribute("CriticalPointsOfP1Map",IsP1Map);
##   <Returns>The critical points of <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="AsP1Map" Arg="rat"/>
DeclareAttribute("AsP1Map",IsScalar);
##   <Returns>The P1 map given by the rational function <A>rat</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="P1MapSL2" Arg="mat"/>
DeclareOperation("P1MapSL2",[IsMatrix]);
##   <Returns>The Möbius P1 map given by the 2x2 matrix <A>mat</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Oper Name="SL2P1Map" Arg="map"/>
DeclareAttribute("SL2P1Map",IsP1Map);
##   <Returns>The matrix of the Möbius P1 map <A>map</A>.</Returns>
## </ManSection>
##
## <ManSection>
##   <Func Name="SetP1Points" Arg="record [prec]"/>
DeclareGlobalFunction("SetP1Points");
##   <Description>
##     Installs a default implementation for P1 points. Fundamentally,
##     a P1 point is a complex number or infinity, with a few extra
##     methods. The argument <A>record</A> is the record describing
##     the floating-point implementation.
##     <P/>
##     Currently, one implementation (the default) is based on pairs
##     of IEEE754 floateans. It is fast, but is limited to 53 bits of
##     precision. It is loaded via <C>SetP1Points(PMCOMPLEX);</C>.
##     <P/>
##     Another implementation, in case the package
##     <Package>Float</Package> is available, is based on MPC complex
##     numbers. It offers unlimited precision, but is much slower. It is
##     loaded via <C>SetP1Points(MPC);</C> or <C>SetP1Points(MPC,prec);</C>.
##   </Description>
## </ManSection>
##
## <#/GAPDoc>

# undocumented for now --
# given three maps alpha,beta,gamma,
# find the pairs (s,t) in [0,1]^2 such that alpha(s)=beta(gamma(t))
DeclareOperation("P1INTERSECT@",[IsP1Map,IsP1Map,IsP1Map]);
#############################################################################

#E p1.gd . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
