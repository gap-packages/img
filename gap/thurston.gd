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

DeclareOperation("ThurstonMatrix", [IsSphereMachine,IsMulticurve]);
DeclareOperation("ThurstonObstruction", [IsSphereMachine,IsMarkedSphere]);

DeclareProperty("IsBicritical", IsObject);

#!!! maybe move somewhere else?
DeclareProperty("IsNonContractingMatrix", IsMatrix);
DeclareOperation("LiftOfConjugacyClass", [IsGroupFRMachine,IsConjugacyClassGroupRep]);

#E thurston.gd . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
