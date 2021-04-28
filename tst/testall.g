LoadPackage("img");
SetInfoLevel(InfoIMG,1);

opts := rec(exitGAP := true, testOptions := rec(compareFunction := "uptowhitespace"));
if not IsBound(MPC) then
  opts.exclude := ["p1-mpc.tst"];
fi;

TestDirectory(DirectoriesPackageLibrary("img","tst"), opts);

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
