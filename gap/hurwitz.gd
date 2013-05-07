Info(InfoIMG,1,"Add doc to hurwitz.gd");
#############################################################################
##
#W hurwitz.gd                                               Laurent Bartholdi
##
#Y Copyright (C) 2013, Laurent Bartholdi
##
#############################################################################
##
##  Constructing Hurwitz maps
##
#############################################################################

#############################################################################
##
#E The Hurwitz problem
##
## <#GAPDoc Label="Hurwitz">
## <ManSection>
##   <Oper Name="BranchedCoveringByMonodromy" Arg="spider, monodromy [, last]"/>
##   <Returns>A record describing a Hurwitz map.</Returns>
##   <Description>
##     If <A>spider</A> is a spider, marked by a group <M>G</M>, and
##     <A>monodromy</A> is a homomorphism from <M>G</M> to a permutation
##     group, this function computes a rational map whose critical values
##     are the vertices of <A>spider</A> and whose monodromy about these
##     critical values is given by <A>monodromy</A>.
##
##     <P/> The returned data are in a record with a field <C>degree</C>, the
##     degree of the map; two fields <C>map</C> and <C>post</C>, describing the
##     desired <M>\mathbb P^1</M>-map --- <C>post</C> is a Möbius
##     transformation, and the composition of <C>map</C> and <C>post</C> is
##     the desired map; and lists <C>zeros</C>, <C>poles</C>
##     and <C>cp</C> describing the zeros, poles and critical points of the
##     map. Each entry in these lists is a record with entries <C>degree</C>,
##     <C>pos</C> and <C>to</C> giving, for each point in the source of
##     <C>map</C>, the local degree and the vertex in <A>spider</A> it maps to.
##
##     <P/> This function requires external programs in the subdirectory
##     "hurwitz" to have been compiled.
## <Example><![CDATA[
## gap> # we'll construct 2d-2 points on the equator, and permutations
## gap> # in order (1,2),...,(d-1,d),(d-1,d),...,(1,2) for these points.
## gap> # first, the spider.
## gap> d := 20;;
## gap> z := List([0..2*d-3], i->P1Point(Exp(i*PMCOMPLEX.constants.2IPI/(2*d-2))));;
## gap> g := SphereGroup(2*d-2);;
## gap> sphere := NewMarkedSphere(z,g);;
## gap> # next, the permutation representation
## gap> perms := List([1..d-1],i->(i,i+1));;
## gap> Append(perms,Reversed(perms));
## gap> perms := GroupHomomorphismByImages(g,SymmetricGroup(d),GeneratorsOfGroup(g),perms);;
## gap> # now compute the map
## gap> BranchedCoveringByMonodromy(spider,perms);
## rec( cp := [ rec( degree := 2, pos := <1.0022-0.0099955i>, to := <vertex 19[ 9, 132, 13, 125 ]> ), 
##       rec( degree := 2, pos := <1.0022-0.0099939i>, to := <vertex 20[ 136, 128, 129, 11 ]> ), 
##       rec( degree := 2, pos := <1.0039-0.0027487i>, to := <vertex 10[ 73, 74, 16, 82 ]> ), 
##       rec( degree := 2, pos := <1.0006-0.0027266i>, to := <vertex 29[ 185, 20, 179, 21 ]> ), 
##       rec( degree := 2, pos := <1.0045-7.772e-05i>, to := <vertex 9[ 24, 77, 17, 72 ]> ), 
##       rec( degree := 2, pos := <1.1739+0.33627i>, to := <vertex 2[ 31, 32, 41, 28 ]> ), 
##       rec( degree := 2, pos := <1.0546+0.12276i>, to := <vertex 3[ 37, 38, 33, 46 ]> ), 
##       rec( degree := 2, pos := <1.026+0.061128i>, to := <vertex 4[ 43, 39, 52, 45 ]> ), 
##       rec( degree := 2, pos := <1.0148+0.03305i>, to := <vertex 5[ 49, 44, 58, 51 ]> ), 
##       rec( degree := 2, pos := <1.0098+0.018122i>, to := <vertex 6[ 55, 50, 64, 57 ]> ), 
##       rec( degree := 2, pos := <1.0071+0.0093947i>, to := <vertex 7[ 61, 62, 71, 59 ]> ), 
##       rec( degree := 2, pos := <1.0055+0.0037559i>, to := <vertex 8[ 67, 68, 63, 69 ]> ), 
##       rec( degree := 2, pos := <1.0035-0.0047633i>, to := <vertex 11[ 79, 75, 88, 81 ]> ), 
##       rec( degree := 2, pos := <1.0031-0.0062329i>, to := <vertex 12[ 85, 80, 94, 87 ]> ), 
##       rec( degree := 2, pos := <1.0029-0.0073311i>, to := <vertex 13[ 91, 86, 100, 93 ]> ), 
##       rec( degree := 2, pos := <1.0027-0.008187i>, to := <vertex 14[ 97, 92, 106, 99 ]> ), 
##       rec( degree := 2, pos := <1.0026-0.008824i>, to := <vertex 15[ 103, 98, 112, 105 ]> ), 
##       rec( degree := 2, pos := <1.0025-0.0092966i>, to := <vertex 16[ 109, 104, 118, 111 ]> ), 
##       rec( degree := 2, pos := <1.0024-0.0096345i>, to := <vertex 17[ 115, 110, 124, 117 ]> ), 
##       rec( degree := 2, pos := <1.0023-0.0098698i>, to := <vertex 18[ 121, 116, 122, 123 ]> ), 
##       rec( degree := 2, pos := <1.0021-0.0098672i>, to := <vertex 21[ 133, 127, 142, 135 ]> ), 
##       rec( degree := 2, pos := <1.002-0.0096298i>, to := <vertex 22[ 139, 134, 148, 141 ]> ), 
##       rec( degree := 2, pos := <1.002-0.0092884i>, to := <vertex 23[ 145, 140, 154, 147 ]> ), 
##       rec( degree := 2, pos := <1.0019-0.0088147i>, to := <vertex 24[ 151, 146, 160, 153 ]> ), 
##       rec( degree := 2, pos := <1.0017-0.008166i>, to := <vertex 25[ 157, 152, 166, 159 ]> ), 
##       rec( degree := 2, pos := <1.0016-0.0073244i>, to := <vertex 26[ 163, 158, 172, 165 ]> ), 
##       rec( degree := 2, pos := <1.0014-0.0061985i>, to := <vertex 27[ 169, 164, 178, 171 ]> ), 
##       rec( degree := 2, pos := <1.0011-0.0047031i>, to := <vertex 28[ 175, 170, 176, 177 ]> ), 
##       rec( degree := 2, pos := <0.99908+0.0038448i>, to := <vertex 31[ 187, 183, 196, 189 ]> ), 
##       rec( degree := 2, pos := <0.99759+0.0094326i>, to := <vertex 32[ 193, 188, 202, 195 ]> ), 
##       rec( degree := 2, pos := <0.99461+0.018114i>, to := <vertex 33[ 199, 194, 208, 201 ]> ), 
##       rec( degree := 2, pos := <0.98944+0.032796i>, to := <vertex 34[ 205, 200, 214, 207 ]> ), 
##       rec( degree := 2, pos := <0.9772+0.058259i>, to := <vertex 35[ 211, 206, 220, 213 ]> ), 
##       rec( degree := 2, pos := <0.94133+0.11243i>, to := <vertex 36[ 217, 212, 226, 219 ]> ), 
##       rec( degree := 2, pos := <0.79629+0.23807i>, to := <vertex 37[ 223, 224, 225, 221 ]> ), 
##       rec( degree := 2, pos := <1+0i>, to := <vertex 30[ 181, 182, 6, 190 ]> ) ], degree := 20, 
##   map := <((-0.32271060393507572-4.3599244721894763i_z)*z^20+(3.8941736874493795+78.415744809040405\
## i_z)*z^19+(-16.808157937605603-665.79436908026275i_z)*z^18+(2.6572296014719168+3545.861245383101i_z\
## )*z^17+(316.57668022762243-13273.931613611372i_z)*z^16+(-1801.6631038749117+37090.818733740503i_z)*\
## z^15+(5888.6033008259928-80172.972599556582i_z)*z^14+(-13500.864941314803+137069.10015838256i_z)*z^\
## 13+(23251.436304923012-187900.36507913063i_z)*z^12+(-31048.192131502536+208077.63047409133i_z)*z^11\
## +(32639.349270133433-186578.17493860485i_z)*z^10+(-27155.791223040047+135145.40893002271i_z)*z^9+(1\
## 7836.343164500577-78489.005444299968i_z)*z^8+(-9153.842142530224+36053.895961137248i_z)*z^7+(3598.6\
## 408777659944-12810.65497539577i_z)*z^6+(-1047.541279063196+3397.470068169695i_z)*z^5+(212.906725643\
## 0024-633.29691376653466i_z)*z^4+(-26.989372105307872+74.040615571896637i_z)*z^3+(1.6073346640110264\
## -4.0860733899027055i_z)*z^2)/(z^18+(-18.034645372692019-0.45671993287358581i_z)*z^17+(153.540499397\
## 49956+7.7811506405054889i_z)*z^16+(-819.9344323563339-62.384270590463998i_z)*z^15+(3077.71530771320\
## 75+312.59552100187739i_z)*z^14+(-8623.1225834872057-1096.4398001099003i_z)*z^13+(18689.34396825033+\
## 2856.8568878158458i_z)*z^12+(-32038.568184053798-5725.9186424029094i_z)*z^11+(44038.148375498437+90\
## 17.0162876593004i_z)*z^10+(-48898.555649389084-11295.156285052604i_z)*z^9+(43964.579894637543+11318\
## .997395732025i_z)*z^8+(-31931.403449371515-9074.2344933443364i_z)*z^7+(18595.347261301522+5786.6036\
## 424805825i_z)*z^6+(-8565.0823844971637-2899.3353634270734i_z)*z^5+(3051.6919509143086+1117.44496422\
## 99487i_z)*z^4+(-811.56293104533825-319.93036282549667i_z)*z^3+(151.69784956523344+64.11787684283315\
## 5i_z)*z^2+(-17.785127700028404-8.0311759305108268i_z)*z+(0.98427999507354302+0.47338721325094818i_z\
## ))>, poles := [ rec( degree := 1, pos := <0.99517+0.30343i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0021+0.11512i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0028+0.05702i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0026+0.030964i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0025+0.016951i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0024+0.0085784i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0024+0.003208i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0023-0.00046905i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0023-0.0030802i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0023-0.0049913i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0023-0.0064163i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0074855i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0082954i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0089048i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0093543i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0096742i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0098869i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0099988i>, to := <vertex 1[ 26, 27, 1, 34 ]> ), 
##       rec( degree := 2, pos := <P1infinity>, to := <vertex 1[ 26, 27, 1, 34 ]> ) ], 
##   post := <((-0.91742065452766763+0.99658449300666985i_z)*z+(0.74087581626192234-1.1339948562200648\
## i_z))/((-0.75451451285920013+0.96940026593933015i_z)*z+(0.75451451285920013-0.96940026593933015i_z)\
## )>, zeros := [ rec( degree := 1, pos := <0.92957+0.28362i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <0.99173+0.11408i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <0.99985+0.056874i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0014+0.030945i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.002+0.016938i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022+0.0085785i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022+0.0032076i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.00046827i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0030802i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0049908i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.006416i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0074855i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0082953i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0089047i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0093542i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0096742i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0098869i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 1, pos := <1.0022-0.0099988i>, to := <vertex 38[ 30, 3, 228, 7 ]> ), 
##       rec( degree := 2, pos := <0+0i>, to := <vertex 38[ 30, 3, 228, 7 ]> ) ] )
## ]]></Example>
##   </Description>
## </ManSection>
##
## <ManSection>
##   <Oper Name="DessinByPermutations" Arg="s0,s1[,sinf]"/>
##   <Returns>A rational map (see <Ref Oper="HurwitzMap"/>) with monodromies <A>s,t</A>.</Returns>
##   <Description>
##     This command computes the Hurwitz map associated with the spider
##     <M>[0,1]\cup[1,\infty]</M>; the monodromy representation is by
##     the permutation <A>s0</A> at <M>0</M> and <A>s1</A> at <M>1</M>. The
##     optional third argument <A>sinf</A> is the monodromy at <M>\infty</M>,
##     and must equal <M>s_0^{-1}s_1^{-1}</M>.
##
##     <P/> The data is returned as a record, with entries <C>degree</C>,
##     <C>map</C>, <C>post</C>, and lists <C>poles</C>, <C>zeros</C>, and
##     <C>above1</C>. Each entry in the list is a record with entries
##     <C>pos</C> and <C>degree</C>.
## <Example><![CDATA[
## gap> DessinByPermutations((1,2),(2,3));
## rec( above1 := [ rec( degree := 2, pos := <1+0i> ),
##                  rec( degree := 1, pos := <-0.5-1.808e-14i> ) ],
##      degree := 3, 
##      map := <(-1.9999999999946754+2.1696575743432764e-13i_z)*z^3+(2.9999999999946754-2.1696575743432764e-13i_z)*z^2>,
##      poles := [ rec( degree := 3, pos := <P1infinity> ) ],
##      post := <z>, 
##      zeros := [ rec( degree := 1, pos := <1.5+5.4241e-14i> ),
##                 rec( degree := 2, pos := <0+0i> ) ] )
## gap> # the Cui example
## gap> DessinByPermutations((1,3,12,4)(5,9)(6,7)(10,13,11)(2,8),
##            (1,5,13,6)(7,10)(2,3)(8,11,12)(4,9),
##            (1,7,11,2)(3,8)(4,5)(9,12,13)(6,10));
## rec( 
##   above1 := [ rec( degree := 2, pos := <1.9952-0.79619i> ), 
##       rec( degree := 2, pos := <0.43236-0.17254i> ), rec( degree := 2, pos := <-0.9863-0.16498i> ),
##       rec( degree := 3, pos := <-0.12749-0.99184i> ), rec( degree := 4, pos := <1+0i> ) ], 
##   degree := 13, 
##   map := <((-6.9809917616400366e-12+0.13002709490708636i_z)*z^13+(-0.68172329304137969-0.8451761169\
## 2062078i_z)*z^12+(4.0903397584184269+0.30979932084028583i_z)*z^11+(-6.3643009040280925+7.5930410215\
## 99336i_z)*z^10+(-5.1732765988942884-16.738910009700096i_z)*z^9+(21.528087032174511+6.11354599010482\
## 5i_z)*z^8+(-15.258776392407746+13.657687016998921i_z)*z^7+(-1.6403496019814323-13.453316297094229i_\
## z)*z^6+(4.4999999996894351+3.3633290741375781i_z)*z^5+(-0.99999999990279009+2.5538451239904557e-11i\
## _z)*z^4)/(z^9+(-4.4999999999400613+3.3633290744267983i_z)*z^8+(1.6403496020557891-13.45331629745540\
## 7i_z)*z^7+(15.258776391831654+13.657687016903173i_z)*z^6+(-21.528087030670253+6.1135459892162567i_z\
## )*z^5+(5.1732765986730511-16.738910007513041i_z)*z^4+(6.3643009027133139+7.593041020468557i_z)*z^3+\
## (-4.0903397575324512+0.30979932067785648i_z)*z^2+(0.68172329288354727-0.8451761166966415i_z)*z+(5.0\
## 734454343833585e-12+0.13002709487107747i_z))>, 
##   poles := [ rec( degree := 2, pos := <1.6127-0.49018i> ), 
##       rec( degree := 2, pos := <0.5-0.04153i> ), rec( degree := 2, pos := <-0.61269-0.49018i> ), 
##       rec( degree := 3, pos := <0.5-0.43985i> ), rec( degree := 4, pos := <P1infinity> ) ], 
##   post := <-z+1._z>, 
##   zeros := [ rec( degree := 2, pos := <1.9863-0.16498i> ), 
##       rec( degree := 3, pos := <1.1275-0.99184i> ), rec( degree := 2, pos := <0.56764-0.17254i> ), 
##       rec( degree := 2, pos := <-0.99516-0.79619i> ), rec( degree := 4, pos := <0+0i> ) ] )
## ]]></Example>
## gap> # IV.5.2 in Granboulan's PhD, the automorphism group of the Mathieu group M_22
## gap> autm22 := Group((1,2,3,4,5,6,7,8,9,10,11)(12,13,14,15,16,17,18,19,20,21,22),
##                      (1,9,3,2)(4,8,17,21)(5,20,19,6)(12,22,16,13)(7,18)(10,11)(14,15),
##                      (3,8)(4,20)(6,18)(7,17)(9,11)(13,15)(16,21));;
## gap> IsomorphismGroups(DerivedSubgroup(autm22),MathieuGroup(22))<>fail;
## true
## gap> DessinByPermutations(autm22.1,autm22.2,autm22.3);
## ...
## gap> # IV.5.3 in Granboulan's PhD, the "extraterrestrial" dessin with group M_24
## gap> m24_ET := Group((1,2,3)(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
##                     (1,7,5)(2,4,23)(3,22,8)(9,21,19)(10,18,12)(13,17,15),
##                     (1,4)(2,22)(3,7)(5,6)(8,21)(9,18)(10,11)(12,17)(13,14)(15,16)(19,20)(23,24));;
## gap> IsomorphismGroups(m24_ET,MathieuGroup(24))<>fail;
## true
## gap> @IMG.hurwitzmesh := 0.4;; # need finer precision
## gap> DessinByPermutations(m24_ET.1,m24_ET.2,m24_ET.3);
## ...
##   </Description>
## </ManSection>
## <#/GAPDoc>
##
DeclareOperation("BranchedCoveringByMonodromy", [IsMarkedSphere,IsGroupHomomorphism]);
DeclareOperation("BranchedCoveringByMonodromy", [IsMarkedSphere,IsGroupHomomorphism,IsRecord]);
DeclareOperation("DessinByPermutations", [IsPerm,IsPerm]);
DeclareOperation("DessinByPermutations", [IsPerm,IsPerm,IsPerm]);
#############################################################################

#E hurwitz.gd . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
