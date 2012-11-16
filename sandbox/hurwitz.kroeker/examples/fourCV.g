################ loading

LoadPackage("float");
LoadPackage("fr");

RereadPackage("fr","hurwitz/gap/utils.gd");
RereadPackage("fr","hurwitz/gap/utils.gi");

RereadPackage("fr","hurwitz/gap/padicLift.gd");
RereadPackage("fr","hurwitz/gap/padicLift.gi");
 
RereadPackage("fr","hurwitz/gap/hurwitz.gd");
RereadPackage("fr","hurwitz/gap/hurwitz.gi");

    
    
################################ four CV example draft #################################################################################

    allMapCandidates := []; # collect results here
    ########### init problem parameters
    
	finiteField := GF(13);	partitions := [ [1,2], [2,1], [2,1], [2,1] ]; 	mapDegree := 3;
	# pairs of rationals approximating real and imaginary part. 
	branchValuesApprox := [ [infinity,infinity], [0,0], [1,0], [0/1, -1/2] ]; 
	
	# reduce critical values to finite field. TODO: pass minimal polynomials to c++ binary instead CV to avoid redundant computation.
	reducedCritivalValueLists := Hurwitz@FR.ReduceCriticalValuesApprox( branchValuesApprox, finiteField );
    strictNormalization := true;   
    
    liftOptions := @PadicLift.LiftOptions();
    liftOptions.setDecimalPrecision(24); 	
	
    for reducedCriticalValues in reducedCritivalValueLists do           
        
        ########## finite field search #################################################      

        mapsModPrime := Hurwitz@FR.FindHurwitzMapModPrime( finiteField  ,partitions, reducedCriticalValues, strictNormalization );
        
        if Size(mapsModPrime)>0 then 
        ########## lift #################################################

            for mapModPrime  in mapsModPrime do 
                problem       := Hurwitz@FR.HurwitzMapSearchProblem( partitions , branchValuesApprox,  strictNormalization );
                mapCandidates := Hurwitz@FR.ApproxComplexHurwitzMaps( problem, mapModPrime[2], finiteField, liftOptions );
                
                Append( allMapCandidates, mapCandidates);
           od;
        ###################################################
        fi;
   od;
  
   # look at allMapCandidates[i].maxResidue to find good maps !
###################################################################################################################################



