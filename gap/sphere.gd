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

DeclareCategory("IsElementOfSphereGroup", IsElementOfFpGroup and IsAssocWord);
DeclareCategoryCollections("IsElementOfSphereGroup");
DeclareCategoryCollections("IsElementOfSphereGroupCollection");
DeclareProperty("IsElementOfSphereGroupFamily", IsElementOfFpGroupFamily);
DeclareProperty("IsSphereGroup", IsFpGroup);
DeclareAttribute("IsomorphismSphereGroup", IsFpGroup);
DeclareAttribute("AsSphereGroup", IsFpGroup);
DeclareAttribute("EulerCharacteristic", IsGroup);
DeclareAttribute("RankOfSphereGroup", IsSphereGroup);
DeclareAttribute("OrderingOfSphereGroup", IsSphereGroup);
DeclareAttribute("ExponentsOfSphereGroup", IsSphereGroup);
DeclareAttribute("IsomorphismFreeGroup", IsSphereGroup);

DeclareGlobalFunction("SphereGroup");
DeclareOperation("ElementOfSphereGroup", [IsFamily, IsAssocWordWithInverse]);

DeclareProperty("IsSphereConjugacyClass", IsAssociativeElementCollection and IsMultiplicativeElementWithInverseCollection);
DeclareProperty("IsSphereConjugacyClassCollection", IsAssociativeElementCollColl and IsMultiplicativeElementWithInverseCollColl);
DeclareSynonym("IsMulticurve", IsSphereConjugacyClassCollection);
DeclareAttribute("PeripheralClasses", IsSphereGroup);
DeclareProperty("IsPeripheral", IsElementOfSphereGroup);
DeclareProperty("IsPeripheral", IsSphereConjugacyClass);
DeclareAttribute("Inverse", IsSphereConjugacyClass);

DeclareOperation("IntersectionNumber", [IsSphereConjugacyClass,IsSphereConjugacyClass]);
DeclareOperation("SelfIntersectionNumber", [IsSphereConjugacyClass]);

DeclareProperty("IsAutomorphismGroupOfSphereGroup", IsAutomorphismGroup);
InstallTrueMethod(IsAutomorphismGroup, IsAutomorphismGroupOfSphereGroup);

DeclareAttribute("EpimorphismToOut", IsAutomorphismGroupOfSphereGroup);

DeclareOperation("AmalgamatedFreeProduct", [IsSphereGroup,IsSphereGroup,IsElementOfSphereGroup,IsElementOfSphereGroup]);
DeclareAttribute("EmbeddingsOfAmalgamatedFreeProduct", IsSphereGroup);

#E sphere.gd . . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
