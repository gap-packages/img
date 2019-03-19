LoadPackage("img");
SetInfoLevel(InfoIMG,1);
dirs := DirectoriesPackageLibrary("img","tst");
Test(Filename(dirs,"chapter-12.tst"), rec(compareFunction := "uptowhitespace"));
Test(Filename(dirs,"chapter-9-a.tst"), rec(compareFunction := "uptowhitespace"));
Test(Filename(dirs,"chapter-9-b.tst"), rec(compareFunction := "uptowhitespace"));
