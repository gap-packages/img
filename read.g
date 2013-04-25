#############################################################################
##
#W read.g                                                   Laurent Bartholdi
##
#H   @(#)$Id$
##
#Y Copyright (C) 2006, Laurent Bartholdi
##
#############################################################################
##
##  This file reads the implementations, and in principle could be reloaded
##  during a GAP session.
#############################################################################

InstallMethod( IsSingleValued,
  "map from whole fp group or free group, given on std. gens: test relators",
  [IsFromFpGroupStdGensGeneralMappingByImages], SUM_FLAGS,
function(hom)
local s,sg,o,gi;
  s:=Source(hom);
  if not IsWholeFamily(s) then
    TryNextMethod();
  fi;
  if IsFreeGroup(s) then
    return true;
  fi;
  sg:=FreeGeneratorsOfFpGroup(s){hom!.genpositions};
  o:=One(Range(hom));
  gi:=MappingGeneratorsImages(hom)[2];
  return ForAll(RelatorsOfFpGroup(s),i->MappedWord(i,sg,gi)=o);
end);

#############################################################################
##
#R Read the install files.
##
ReadPackage("img", "gap/helpers.gi");
ReadPackage("img", "gap/complex.gi");
ReadPackage("img", "gap/p1.gi");
ReadPackage("img", "gap/p1_ieee754.gi");
if @.dll then
    SetP1Points(PMCOMPLEX);
elif IsPackageMarkedForLoading("float","") then
    if IsBound(MPC) then
        SetP1Points(MPC,100);
    elif IsBound(CXSC) then
        SetP1Points(CXSC);
    else
        Info(InfoPackageLoading, 1, "You installed the Float package but compiled neither the Float nor the IMG DLL. That's probably not what you wanted. I'll disable the P1Points code.");
        @.ro := 1.0_l; # minimal defaults to shut up warnings -- won't be usable.
        @.o := NewFloat(IsPMComplex,1);
        @.reps := IEEE754FLOAT.constants.EPSILON;
    fi;
else
    Info(InfoPackageLoading, 2, "You didn't install the Float package, and didn't compile the IMG DLL. I'll disable the P1Points code.");
    @.ro := 1.0_l; # minimal defaults to shut up warnings -- won't be usable.
    @.o := NewFloat(IsPMComplex,1);
    @.reps := IEEE754FLOAT.constants.EPSILON;
fi;
ReadPackage("img", "gap/sphere.gi");
ReadPackage("img", "gap/triangulations.gi");
ReadPackage("img", "gap/spider.gi");
ReadPackage("img", "gap/machine.gi");
ReadPackage("img", "gap/markedsphere.gi");
ReadPackage("img", "gap/hurwitz.gi");
ReadPackage("img", "gap/thurston.gi");
ReadPackage("img", "gap/examples.gi");
#############################################################################

while not IsEmpty(POSTHOOK@img) do Remove(POSTHOOK@img)(); od;
Unbind(POSTHOOK@img);

if IsBound(IO_Pickle) then
    ReadPackage("img","gap/pickle.g");
else
    if not IsBound(IO_PkgThingsToRead) then
        IO_PkgThingsToRead := [];
    fi;
    Add(IO_PkgThingsToRead, ["img","gap/pickle.g"]);
fi;

MAKEP1EPS@();

#E read.g . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
