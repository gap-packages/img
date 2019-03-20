LoadPackage("img");
SetInfoLevel(InfoIMG,1);
TestDirectory(DirectoriesPackageLibrary( "img", "tst" ),
  rec(exitGAP     := true) );

FORCE_QUIT_GAP(1); # if we ever get here, there was an error

