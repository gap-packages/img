#############################################################################
##
#W helpers.gd                                               Laurent Bartholdi
##
#Y Copyright (C) 2012-2013, Laurent Bartholdi
##
#############################################################################
##
##  This file contains helper code for functionally recursive groups,
##  in particular related to the geometry of groups.
##
#############################################################################

## <#GAPDoc Label="Helpers">
DeclareGlobalFunction("Mandel");
## <ManSection>
##   <Func Name="Mandel" Arg="[map]"/>
##   <Returns>Calls the external program <File>mandel</File>.</Returns>
##   <Description>
##     This function starts the external program <File>mandel</File>, by Wolf Jung.
##     The program is searched for along the standard PATH; alternatively,
##     its location can be set in the string variable EXEC@fr.mandel.
##
##     <P/> When called with no arguments, this command returns starts
##     <File>mandel</File> in its default mode. With a rational map as argument, it
##     starts <File>mandel</File> pointing at that rational map.
##
##     <P/> More information on <File>mandel</File> can be found
##     at <URL>http://www.mndynamics.com</URL>.
##   </Description>
## </ManSection>
##
DeclareOperation("NonContractingSubmatrix", [IsMatrix]);
## <ManSection>
##   <Oper Name="NonContractingSubmatrix" Arg="mat"/>
##   <Returns><K>fail</K> or a list of indices <C>l</C> such that <C>mat{l}{l}</C> is irreducible and non-contracting</Returns>
##   <Description>
##     This function computes a minimal submatrix whose spectral radius is <M>\geq1</M>. If none exists,
##     it returns <K>fail</K>.
## <Example><![CDATA[
## gap> NonContractingSubmatrix([[2]]);
## [ 1 ]
## gap> NonContractingSubmatrix([[1/2]]);
## fail
## gap> NonContractingSubmatrix([[0,1],[1,0]]);
## [ 1, 2 ]
## gap> NonContractingSubmatrix([[0,1],[0,1]]);
## [ 2 ]
## ]]></Example>
##   </Description>
## </ManSection>
## <#/GAPDoc>

#E helpers.gd . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
