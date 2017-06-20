#############################################################################
##
#W examples.gd                                              Laurent Bartholdi
##
#Y Copyright (C) 2012-2013, Laurent Bartholdi
##
#############################################################################
##
##  All interesting examples of IMG's I came through
##
#############################################################################

#############################################################################
##
#E PoirierExamples
##
## <#GAPDoc Label="Poirier">
DeclareGlobalFunction("PoirierExamples");
## <ManSection>
##   <Func Name="PoirierExamples" Arg="..."/>
##   <Description>
##     The examples from Poirier's paper <Cite Key="math.DS/9305207"/>.
##     See details under <Ref Oper="PolynomialSphereMachine"/>; in particular,
##     <C>PoirierExamples(1)</C> is the Douady rabbit map.
##   </Description>
## </ManSection>
## <#/GAPDoc>
##
#############################################################################

#############################################################################
##
#E SphereTorusMap
##
## <#GAPDoc Label="SphereTorusMap">
DeclareOperation("SphereTorusMap",[IsPosInt,IsPosInt]);
## <ManSection>
##   <Oper Name="SphereTorusMap" Arg="m n"/>
##   <Returns>A sphere map doubly covered by a torus endomorphism</Returns>
##   <Description>
##     This function computes the sphere machine representation of
##     the quotient of torus endomorphism <M>[[m,0],[0,n]]</M> by the 
##     Weierstrass <M>\wp</M>-function for the square lattice.
##   </Description>
## </ManSection>
## <#/GAPDoc>
##
#############################################################################

#E examples.gd. . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
