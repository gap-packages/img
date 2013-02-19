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

DeclareCategory("IsElementOfSphereGroup", IsElementOfFpGroup);
DeclareProperty("IsElementOfSphereGroupFamily", IsElementOfFpGroupFamily);
DeclareProperty("IsSphereGroup", IsFpGroup);

DeclareGlobalFunction("SphereGroup");
DeclareOperation("ElementOfSphereGroup", [IsFamily, IsAssocWordWithInverse]);

DeclareProperty("IsSphereConjugacyClass", IsAssociativeElementCollection and IsMultiplicativeElementWithInverseCollection);
DeclareOperation("IntersectionNumber", [IsSphereConjugacyClass,IsSphereConjugacyClass]);

#E sphere.gd . . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
