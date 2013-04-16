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

DeclareOperation("ThurstonAlgorithm", [IsSphereMachine]);
DeclareOperation("P1MapBySphereMachine", [IsSphereMachine]);

DeclareOperation("FindThurstonObstruction", [IsElementOfSphereGroupCollection]);

DeclareProperty("IsBicritical", IsObject);

#E thurston.gd . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
