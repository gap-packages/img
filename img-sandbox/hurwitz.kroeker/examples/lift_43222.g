################ loading

LoadPackage("float");
LoadPackage("fr");

RereadPackage("fr","hurwitz/gap/utils.gd");
RereadPackage("fr","hurwitz/gap/utils.gi");

RereadPackage("fr","hurwitz/gap/padicLift.gd");
RereadPackage("fr","hurwitz/gap/padicLift.gi");
 
RereadPackage("fr","hurwitz/gap/hurwitz.gd");
RereadPackage("fr","hurwitz/gap/hurwitz.gi");


################################## [4,3,2,2,2]- example (lifting)  draft ################################## 



hmsProblem := Hurwitz@FR.HurwitzMapSearchProblem( [[4,3,2,2,2], [3,4,2,2,2], [3,2,4,2,2]], 
                                                [[infinity,infinity], [0,0], [1,0]], 
                                                true);
    
    #### init lift parameters 
    finiteField := GF(11);
    rng := PolynomialRing( finiteField  ,["x"] );
    indeterminates := IndeterminatesOfPolynomialRing(rng);
    x := indeterminates[1];
    polTuple := [];
 
    Append( polTuple, [         (x-5)^3*(x^3 +3*x^2 +2*x +3)^2] );
    Append( polTuple, [   (x)^4*(x+3)^3*(x^3        -3*x -5)^2] );
    Append( polTuple, [ (x-1)^4*(x-3)^3*(x^3         -2*x-3)^2] );   
                                               
    opts := @PadicLift.LiftOptions();   
    opts.setDecimalPrecision (60);  
                                                    
    ##### lift 
    
    lifter := Hurwitz@FR.HurwitzMapLifter(polTuple, finiteField, hmsProblem);  
    approxHurwitzMaps := lifter.computeApproxHurwitzMaps(opts);  
    
    ################ check result #########################
    for mapData in approxHurwitzMaps do    
       Assert(0, mapData.maxResidue<1.0e-15);
    od;
