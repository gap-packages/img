if fail = LoadPackage("AutoDoc", ">= 2016.01.21") then
    Error("AutoDoc 2016.01.21 or newer is required");
fi;
AutoDoc(rec(
    gapdoc := rec(
        main:="img.xml",
        files:=["PackageInfo.g"],
    )
));

QUIT;
