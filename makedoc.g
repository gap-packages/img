#if fail = LoadPackage("AutoDoc", ">= 2016.01.21") then
#    Error("AutoDoc 2016.01.21 or newer is required");
#fi;
#AutoDoc(rec(gapdoc := rec(files:=["PackageInfo.g"])));

MakeGAPDocDoc("doc","img",
        ["../gap/complex.gd","../gap/examples.gd","../gap/helpers.gd",
         "../gap/hurwitz.gd","../gap/machine.gd","../gap/markedsphere.gd",
         "../gap/p1.gd","../gap/sphere.gd","../gap/thurston.gd",
         "../gap/triangulations.gd","../PackageInfo.g"],"img","../../..");
CopyHTMLStyleFiles("doc");
GAPDocManualLab("img");

QUIT;
