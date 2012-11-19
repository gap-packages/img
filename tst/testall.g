LoadPackage("img");
SetInfoLevel(InfoIMG,1);
dirs := DirectoriesPackageLibrary("img","tst");
ReadTest(Filename(dirs,"chapter-12.tst"));
ReadTest(Filename(dirs,"chapter-9-a.tst"));
ReadTest(Filename(dirs,"chapter-9-b.tst"));
