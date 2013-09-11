#############################################################################
##
##  PackageInfo.g for the package `IMG'                    Laurent Bartholdi
##
SetPackageInfo( rec(
PackageName := "IMG",
Subtitle := "Computations with iterated monodromy groups",
Version := "0.0.4",
Date := "11/09/2013",
## <#GAPDoc Label="Version">
## <!ENTITY Version "0.0.4">
## <!ENTITY Date "11/09/2013">
## <#/GAPDoc>
ArchiveURL := Concatenation("https://github.com/laurentbartholdi/img/archive/",~.Version),
ArchiveFormats := ".tar.gz",
Persons := [
  rec(
    LastName      := "Bartholdi",
    FirstNames    := "Laurent",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "laurent.bartholdi@gmail.com",
    WWWHome       := "http://www.uni-math.gwdg.de/laurent",
    PostalAddress := Concatenation( [
                       "Mathematisches Institut\n",
                       "Bunsenstraße 3-5\n",
                       "D-37073 Göttingen\n",
                       "Germany" ] ),
    Place         := "Göttingen",
    Institution   := "Georg-August Universität zu Göttingen"
  )
],

Status := "deposited",
CommunicatedBy := "Götz Pfeiffer (NUI Galway)",
AcceptDate := "",

README_URL := "http://laurentbartholdi.github.com/img/README.img",
PackageInfoURL := "http://laurentbartholdi.github.com/img/PackageInfo.g",
AbstractHTML := "The <span class=\"pkgname\">IMG</span> package allows \
   GAP to manipulate iterated monodromy groups",
PackageWWWHome := "http://laurentbartholdi.github.com/img/",

PackageDoc := rec(
  BookName  := "IMG",
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Iterated monodromy groups",
  ArchiveURLSubset := ["doc"],
  Autoload  := true
),

Dependencies := rec(
  GAP := ">=4.6.0",
  NeededOtherPackages := [["FR",">=2.0.0"],
                      ["GAPDoc",">=1.0"]],
  SuggestedOtherPackages := [["Float",">=0.4"]],

  # for compilation of the external module, one needs:
  # gcc, libcblas, javac, appletviewer.
  ExternalConditions := []
),

AvailabilityTest := ReturnTrue,
                    
BannerString := Concatenation("Loading IMG ", String( ~.Version ), #CallFuncList(function() if Filename(DirectoriesPackagePrograms("img"),"img_dll.so")=fail then return ""; else return "with DLL"; fi; end,[]),
  " ...\n"),

Autoload := false,
TestFile := "tst/testall.g",
Keywords := ["iterated monodromy group"]
));
