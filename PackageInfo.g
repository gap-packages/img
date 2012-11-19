#############################################################################
##
##  PackageInfo.g for the package `IMG'                    Laurent Bartholdi
##
SetPackageInfo( rec(
PackageName := "IMG",
Subtitle := "Computations with iterated monodromy groups",
Version := "0.0.0",
## <#GAPDoc Label="Version">
## <!ENTITY0.0.0
## <#/GAPDoc>
Date := "19/11/2012",
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

README_URL := "http://www.uni-math.gwdg.de/laurent/IMG/README.img",
PackageInfoURL := "http://www.uni-math.gwdg.de/laurent/IMG/PackageInfo.g",
AbstractHTML := "The <span class=\"pkgname\">IMG</span> package allows \
   GAP to manipulate iterated monodromy groups",
PackageWWWHome := "http://www.uni-math.gwdg.de/laurent/IMG/",

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
  GAP := ">=4.5.0",
  NeededOtherPackages := [["FR",">=2.0.0"],
                      ["GAPDoc",">=1.0"]],
  SuggestedOtherPackages := [["Float",">=0.4"]],

  # for compilation of the external module, one needs:
  # gcc, gfortran, libcblas, libgsl, javac, appletviewer.
  ExternalConditions := []
),

AvailabilityTest := ReturnTrue,
                    
BannerString := Concatenation("Loading IMG ", String( ~.Version ), " ...\n"),

Autoload := false,
TestFile := "tst/testall.g",
Keywords := ["iterated monodromy group"]
));
