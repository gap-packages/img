#############################################################################
##
#W init.g                                                   Laurent Bartholdi
##
#Y Copyright (C) 2012, Laurent Bartholdi
##
#############################################################################
##
##  This file reads the declarations of the packages' new objects
##
#############################################################################

POSTHOOK@img := []; # to be processed at the end

BindGlobal("@", rec()); # a record to store locals in the package

#############################################################################
##
#I Create info class to be able to debug loading
##
InfoIMG := NewInfoClass("InfoIMG");
SetInfoLevel(InfoIMG, 1);
#############################################################################

#############################################################################
##
#R Read the declaration files.
##
ReadPackage("img", "gap/helpers.gd");
ReadPackage("img", "gap/complex.gd");
ReadPackage("img", "gap/p1.gd");
ReadPackage("img", "gap/sphere.gd");
ReadPackage("img", "gap/triangulations.gd");
ReadPackage("img", "gap/machine.gd");
ReadPackage("img", "gap/markedsphere.gd");
ReadPackage("img", "gap/hurwitz.gd");
ReadPackage("img", "gap/thurston.gd");
ReadPackage("img", "gap/examples.gd");

CallFuncList(function()
    local dirs, dll, w;
    dirs := DirectoriesPackagePrograms("img");
    dll := Filename(dirs,"img_dll.so");
    if dll=fail then
        dll := Filename(dirs[1],"img_dll.so");
        for w in ["FIND_BARYCENTER","FIND_RATIONALFUNCTION",
                "C22P1POINT","P1POINT2C2","P1POINT2STRING","EQ_P1POINT",
                "P1SPHERE","SPHEREP1","SPHEREP1Y","P1BARYCENTRE",
                "P1ANTIPODE","P1MIDPOINT","P1DISTANCE","XRATIO","P1XRATIO",
                "CLEANEDP1POINT","P1CIRCUMCENTRE","LT_P1POINT",
                "P1MAPBYCOEFFICIENTS_IEEE754","P1INTERSECT_IEEE754",
                "MAT2P1MAP","P1MAP2MAT","P1MAP3","P1MAP2","P1PATH",
                "CLEANEDP1MAP","INVERSEP1MAP","COMPOSEP1MAP","P1IMAGE",
                "P1PREIMAGES","P1MAPCRITICALPOINTS","P1MAPCONJUGATE",
                "P1MAPBYZEROSPOLES","P1MAPPRIMITIVE","P1MAPDERIVATIVE",
                "P1MAPNUMER","P1ROTATION_IEEE754","P1MAPDENOM",
                "P1MAPISPOLYNOMIAL","STRINGS2P1POINT","DEGREEOFP1MAP",
                "P1MAPNUMERATOR","P1MAPDENOMINATOR","P1MAP_SUM","P1MAP_DIFF",
                "P1MAP_PROD","P1MAP_QUO","P1MAP_INV","P1MAP_AINV"] do
            CallFuncList(function(w)
                BindGlobal(w, function(arg)
                    Error("You need to compile ",dll," before using ",w,"\nYou may compile it with './configure && make' in ",PackageInfo("img")[1].InstallationPath,"\n...");
                end);
            end,[w]);
        od;
        @.dll := false;
    else
        LoadDynamicModule(dll);
        @.dll := true;
    fi;
end,[]);

#############################################################################

#E init.g . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
