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

DeclareGlobalFunction("NormalizedQuadraticP1Map");

DeclareOperation("P1Map", [IsRingElement, IsFRMachine]);
DeclareProperty("IsBicritical", IsObject);

DeclareOperation("FindThurstonObstruction", [IsElementOfSphereGroupCollection]);

#E thurston.gd . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
