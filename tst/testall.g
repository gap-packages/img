LoadPackage("img");
SetInfoLevel(InfoIMG,1);
dirs := DirectoriesPackageLibrary("img","tst");
Test(Filename(dirs,"chapter-12.tst"));
Test(Filename(dirs,"chapter-9-a.tst"));
Test(Filename(dirs,"chapter-9-b.tst"));
Test(Filename(dirs,"p1.tst"));
Test(Filename(dirs,"p1-mpc.tst"));
QUIT_GAP();
