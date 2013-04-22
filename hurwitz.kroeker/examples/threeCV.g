################ loading

LoadPackage("float");
LoadPackage("fr");

RereadPackage("fr","hurwitz/gap/utils.gd");
RereadPackage("fr","hurwitz/gap/utils.gi");

RereadPackage("fr","hurwitz/gap/padicLift.gd");
RereadPackage("fr","hurwitz/gap/padicLift.gi");
 
RereadPackage("fr","hurwitz/gap/hurwitz.gd");
RereadPackage("fr","hurwitz/gap/hurwitz.gi");


################################ three CV example #############################################################################

    approxHurwitzMapCandidates:=[];    #result accumulator
    
   ############# init parameters
   
	finiteField := GF(11);  permutations := [(1,2,3),(2,3),(1,2)];

	complexCriticalValuesApprox := [ [infinity,infinity],     [0,0],           [ 1/1, 0 ]     ]; 
	modPrimeCriticalValues      := [      infinity,      Zero(finiteField),  One(finiteField) ];
    hurwitzMapSearchProblem := Hurwitz@FR.HurwitzMapSearchProblem( permutations , complexCriticalValuesApprox);        
	
   ############# finite field search
   
	mapsModPrime := Hurwitz@FR.FindHurwitzMapModPrime( finiteField , permutations, modPrimeCriticalValues );
	 
   ############# lift and approximate Hurwitz map candidates
	
    for  mapModPrime  in mapsModPrime do 
        mapCandidates := Hurwitz@FR.ApproxComplexHurwitzMaps( hurwitzMapSearchProblem, 
                                                      mapModPrime[2], 
                                                      finiteField, 
                                                      @PadicLift.LiftOptions() 
                                                    );
        Append( approxHurwitzMapCandidates, mapCandidates); 
    od;
    
   ############# check if result matches expectations ("Algorithmic construction of Hurwitz maps",  page 3)
   
    z := approxHurwitzMapCandidates[1].indeterminate;
    Assert(0, Degree( ( 3*z^2+(-2.0*z^3) ) / approxHurwitzMapCandidates[1].map ) =0);  
   
#################################################################################################################################



   
   
