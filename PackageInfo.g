#############################################################################
##
##  PackageInfo.g for the package `IMG'                    Laurent Bartholdi
##
SetPackageInfo( rec(
PackageName := "IMG",
Subtitle := "Computations with iterated monodromy groups",
Version := "0.3.0",
Date := "04/02/2021",
License := "GPL-2.0-or-later",
## <#GAPDoc Label="Version">
## <!ENTITY Version "0.3.0">
## <!ENTITY Date "04/02/2021">
## <#/GAPDoc>
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

AbstractHTML := "The <span class=\"pkgname\">IMG</span> package allows \
   GAP to manipulate iterated monodromy groups",

SourceRepository:= rec(Type := "git", URL := "https://github.com/gap-packages/img"),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL, "/releases/download/v", ~.Version, "/", ~.PackageName, "-", ~.Version ),
ArchiveFormats  := ".tar.gz",
PackageWWWHome  := Concatenation( "https://gap-packages.github.io/", ~.PackageName ),
README_URL      := Concatenation( ~.PackageWWWHome, "/README.md" ),
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),

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
                      ["GAPDoc",">=1.0"],
		      ["IO",">=4.0"]],
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
