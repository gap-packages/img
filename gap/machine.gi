#############################################################################
##
#W machine.gi                                               Laurent Bartholdi
##
#Y Copyright (C) 2013, Laurent Bartholdi
##
#############################################################################
##
##  Iterated monodromy groups
##
#############################################################################

InstallMethod(IsSphereMachine, "(IMG) for an FR machine",
        [IsGroupFRMachine],
        function(M)
    local g;
    g := StateSet(M);
    return HasIsSphereGroup(g) and IsSphereGroup(g);
end);

if false then # optimized in C code
BindGlobal("NFFUNCTION_FR", function(rel,dir,word)
    local push_letter, match_pos, posind, negind, result, resulti, i, j, match, matchlen, n, vi;
    
    push_letter := function(v)
        if resulti>0 and v=-result[resulti] then
            resulti := resulti-1;
        else
            resulti := resulti+1;
            result[resulti] := v;
        fi;
    end;
    
    match_pos := function(v)
        if v>0 then
            return posind[v];
        else
            return negind[-v];
        fi;
    end;
    
    # word is an integer lists. dir is true/false.
    # rel is a list of lists: square of positive relator+square of negative
    # relator; positions in 1st of letter i; position in 1st of letter -i
    # if dir=true, replace all (>=1/2)-cyclic occurrences of rel in word by the shorter half
    # if dir=false, replace all occurrences of the last generator in word by the corresponding bit of rel
    
    posind := rel[2]; negind := rel[3]; rel := rel[1];
    n := Length(posind);
    
    i := 0;
    resulti := 0;
    match := 0;
    matchlen := 0;
    result := [];
    
    while i<Length(word) do
        # we produced result[1..resulti] as the compressed version of word[1..i].
        # additionally, matchlen is maximal such that
        # rel[match..match+matchlen-1] = result[resulti-matchlen+1..resulti]
        i := i+1;
        vi := word[i];
        if dir then
            if resulti>0 and vi = -result[resulti] then
                # pop letter, and update match
                resulti := resulti-1;
                matchlen := matchlen-1;
                if matchlen=0 and resulti>0 then
                    match := match_pos(result[resulti]);
                    matchlen := 1;
                    while resulti>matchlen and result[resulti-matchlen]=rel[match+n-1] do
                        matchlen := matchlen+1;
                        match := match-1;
                        if match=0 then match := n; fi;
                    od;
                else
                    match := 0;
                fi;
            else
                push_letter(vi);
                if match>0 and vi = rel[match+matchlen] then
                    matchlen := matchlen+1;
                    if matchlen >= QuoInt((n+2-QuoInt(match,2*n)),2) then # more than half, or exactly half and negatives
                        resulti := resulti-matchlen;
                        for j in [match+n-1,match+n-2..match+matchlen] do
                            push_letter(-rel[j]);
                        od;
                        matchlen := n-matchlen;
                        match := 4*n+1 - (match+n-1);
                    fi;
                else
                    matchlen := 1;
                    match := match_pos(vi);
                fi;
            fi;
        else
            if vi=n then
                match := negind[n];
                for j in [match+1..match+n-1] do
                    push_letter(rel[j]);
                od;
            elif vi=-n then
                match := posind[n];
                for j in [match+1..match+n-1] do
                    push_letter(rel[j]);
                od;
            else
                push_letter(vi);
            fi;
        fi;
    od;
    return result{[1..resulti]};
end);
fi;

# takes into account the IMG relator to try harder to express a machine
# as a subfrmachine
InstallMethod(SubFRMachine, "(IMG) for a sphere machine and a map",
        [IsSphereMachine, IsGroupHomomorphism],
        function(M,f)
    local S, trans, out, gen, pi, machine, decomp, adder;

    S := StateSet(M);
    while S<>Range(f) or not IsSphereGroup(Source(f)) do
        Error("SubFRMachine: map must be homomorphism from sphere group to stateset of machine\n");
    od;
    pi := WreathRecursion(M);
    trans := [];
    out := [];
    
    for gen in GeneratorsOfGroup(Source(f)) do
        decomp := pi(gen^f);
        Add(out,decomp[2]);
        decomp := List(decomp[1],x->PreImagesRepresentative(f,x));
        if fail in decomp then return fail; fi;
        Add(trans,decomp);
    od;
    machine := FRMachineNC(FamilyObj(M),Source(f),trans,out);
    IsSphereMachine(machine);

    if HasAddingElement(M) then
        adder := PreImagesRepresentative(f,InitialState(AddingElement(M)));
        if adder<>fail then
            SetAddingElement(machine,FRElement(machine,adder));
        fi;
    fi;
    return machine;
end);

InstallMethod(ViewString, "(IMG) for a sphere machine",
        [IsSphereMachine and IsFRMachineStdRep],
        M->CONCAT@FR("<sphere machine with alphabet ", AlphabetOfFRObject(M), " on ", StateSet(M), " / ",RelatorsOfFpGroup(StateSet(M)),">"));

InstallMethod(DisplayString, "(IMG) for a sphere machine",
        [IsPolynomialSphereMachine and IsFRMachineStdRep],
        M->CONCAT@FR(DISPLAYFRMACHINE@FR(M),"Relators: ",RelatorsOfFpGroup(StateSet(M)),"\n"));
#############################################################################

#############################################################################
##
#A AsSphereFRMachine
##
BindGlobal("ADDERPOS@", function(M)
    return Position(GeneratorsOfFRMachine(M),InitialState(AddingElement(M)));
end);

BindGlobal("ISADDER@", function(M,w)
    local r, c;
    r := WreathRecursion(M)(w);
    c := Cycles(PermList(r[2]),AlphabetOfFRObject(M));
    return Length(c)=1 and # transitive element
           IsConjugate(StateSet(M),w,Product(r[1]{c[1]}));
end);

InstallGlobalFunction(NewSphereMachine,
        function(arg)
    local r, rr, relators, machine, states, data, len, adder;
    #!!! Guess the orders, the relator, and the adder if not supplied.

    len := Length(arg);
    while not '=' in arg[len] do
        len := len-1;
        while len=0 do Error("I couldn't find any definition of generator"); od;
    od;

    while '=' in arg[Length(arg)] do Error("I couldn't find any relator"); od;
    
#    adder := ValueOption("AddingElement");
#    if adder<>fail then
#        while not IsString(adder) do
#            Error("Required form is AddingElement := <string>");
#        od;
#    fi;
    
    machine := CallFuncList(NewGroupFRMachine,arg{[1..len]});
    data := rec(holdername := RANDOMNAME@FR());
    BindGlobal(data.holdername, StateSet(machine));
    relators := List(arg{[len+1..Length(arg)]},w->STRING_WORD2GAP@FR(List(GeneratorsOfFRMachine(machine),String),"GeneratorsOfGroup",data,w));
    states := AsSphereGroup(StateSet(machine)/relators);
    
#    if adder<>fail then
#        MakeReadWriteGlobal(data.holdername);
#        UnbindGlobal(data.holdername);
#        BindGlobal(data.holdername, states);
#        adder := STRING_WORD2GAP@FR(List(GeneratorsOfGroup(states),String),"GeneratorsOfGroup",data,adder);
#    fi;
    
    MakeReadWriteGlobal(data.holdername);
    UnbindGlobal(data.holdername);
    
    for r in relators do
        rr := WreathRecursion(machine)(r);
        while not ISONE@FR(rr[2]) or ForAny(rr[1],w->not IsOne(ElementOfSphereGroup(FamilyObj(One(states)),w))) do
            Error("Relation ",r," does not hold in the machine");
        od;
    od;
    
    machine := AsSphereMachine(machine,states);
    
#    if adder<>fail then
#        while not ISADDER@(machine,adder) do
#            Error("Element ",adder," is not an adder");
#        od;
#        SetAddingElement(machine,FRElement(machine,adder));
#    fi;
    len := Length(GeneratorsOfGroup(states));
    if ISADDER@(machine,states.(len)) then
        SetAddingElement(machine,FRElement(machine,states.(len)));
    fi;
    
    return machine;
end);

InstallMethod(AsSphereMachine, "(IMG) for a group FR machine and a word",
        [IsGroupFRMachine,IsAssocWord],
        function(M,w)
    return AsSphereMachine(M,AsSphereGroup(StateSet(M)/[w]));
end);

InstallMethod(AsSphereMachine, "(IMG) for a group FR machine and a sphere group",
        [IsGroupFRMachine,IsSphereGroup],
        function(M,G)
    return M^GroupHomomorphismByImages(StateSet(M),G);
end);   

InstallMethod(AsSphereMachine, "(IMG) for a group FR machine",
        [IsGroupFRMachine],
        function(M)
    local f, g, epi, w, p, trans, perm, N;
    
    # try running a spider algorithm to discover a good ordering
    w := SPIDERALGORITHM@(M);
    f := M!.free;

    perm := [fail];
    if w<>fail and w.minimal then # try that ordering first
        Add(perm,GeneratorsOfGroup(f){w.ordering},1);
    fi;
    for p in perm do
        if p=fail then # add now all permutations fixing first letter
            for p in SymmetricGroup([2..Length(GeneratorsOfGroup(f))]) do
                Add(perm,Permuted(GeneratorsOfGroup(f),p));
            od;
            continue;
        fi;
        g := AsSphereGroup(f/[Product(p,One(StateSet(M)))]);
        epi := GroupHomomorphismByImages(f,g);
        N := M^epi;
        if N<>fail then
            SetCorrespondence(N,epi);
            return N;
        fi;
    od;
    return fail;
end);
#############################################################################

#############################################################################
##
#A Kneading machines
##
BindGlobal("ISTREELIKEPERMUTATIONLIST@", function(S,A)
    local s, t;
    S := Concatenation(List(S,x->Cycles(x,A)));
    while Size(S)>1 do
        s := Remove(S);
        t := PositionProperty(S,x->Length(Intersection(x,s))=1);
        if t=fail then return false; fi;
        S[t] := Union(s,S[t]);
    od;
    return Set(S[1])=A;
end);

InstallMethod(IsKneadingMachine, "(IMG) for a Mealy machine",
        [IsMealyMachine],
        function(M)
    local S, s, t;

    S := GeneratorsOfFRMachine(M);
    for s in S do
        if not IsOne(FRElement(M,s)) then
            if Sum(List(S,t->Number(AlphabetOfFRObject(M),a->Transition(M,t,a)=s)))>1 then
                return false;
            fi;
        fi;
        for t in Cycles(PermList(Output(M,s)),AlphabetOfFRObject(M)) do
            if Number(t,x->not IsOne(FRElement(M,Transition(M,s,t))))>1 then
                return false;
            fi;
        od;
    od;
    return ISTREELIKEPERMUTATIONLIST@(List(S,x->PermList(Output(M,x))),AlphabetOfFRObject(M));
end);

BindGlobal("PLANAREMBEDDINGMEALYMACHINE@", function(M,justone)
    local S, aS, result, a, x, perm;

    if not IsKneadingMachine(M) then
        return [];
    fi;
    S := M{GeneratorsOfFRMachine(M)};
    aS := Filtered([1..Length(S)],i->not IsOne(S[i]));
    result := [];
    perm := PermutationsList(aS);
    
    # use a spider algorithm to try to guess a good ordering
    a := SPIDERALGORITHM@(AsGroupFRMachine(M));
    if a<>fail and a.minimal then
        Add(perm,a.ordering,1); # put it at front
    fi;
    for a in PermutationsList(aS) do
        x := List([1..Length(aS)],i->Product(S{a{[i+1..Length(aS)]}},S[1]^0)
                  *Product(S{a{[1..i]}}));
        if State(x[1]^Length(AlphabetOfFRObject(M)),AlphabetOfFRObject(M)[1]) in x then
            if justone then return a; fi;
            Add(result,a);
        fi;
    od;
    return result;
end);

InstallMethod(IsPlanarKneadingMachine, "(IMG) for a Mealy machine",
        [IsMealyMachine],
        M->PLANAREMBEDDINGMEALYMACHINE@(M,true)<>[]);

InstallMethod(AsPolynomialSphereMachine, "(IMG) for a Mealy machine",
        [IsMealyMachine],
        function(M)
    local a;
    a := PLANAREMBEDDINGMEALYMACHINE@(M,true);
    if a=[] then return fail; fi;
    M := AsGroupFRMachine(M);
    SetAddingElement(M,FRElement(M,Product(a,x->x^Correspondence(M))));
    return M;
end);
#############################################################################

#############################################################################
##
#M PROD
##
BindGlobal("COPYADDER@", function(M,N)
    SetAddingElement(M,FRElement(M,InitialState(AddingElement(N))));
end);

if false then
# we don't use this anymore -- it seems actually to slow things down
BindGlobal("NORMALIZEHOMOMORPHISM@", function(f)
    local g, mapi, sf, rf, gens;
    if not HasMappingGeneratorsImages(f) then
        return f;
    fi;
    mapi := MappingGeneratorsImages(f);
    sf := Source(f);
    rf := Range(f);
    gens := GeneratorsOfGroup(sf);
    if mapi[1]=gens then
        return f;
    fi;
    g := Group(mapi[1]);
    mapi := mapi[2];
    return GroupHomomorphismByImagesNC(sf,rf,gens,List(gens,x->Product(AsWordLetterRepInGenerators(x,g),i->mapi[AbsInt(i)]^SignInt(i),One(rf))));
end);
fi;

InstallMethod(\*, "(IMG) for an FR machine and a mapping",
        [IsFRMachine and IsFRMachineStdRep, IsMapping],
        function(M,f)
    local S, N, x;
    S := StateSet(M);
    if S<>Source(f) or S<>Range(f) then
        Error("\*: source, range and stateset must be the same\n");
    fi;
    N := FRMachineNC(FamilyObj(M),S,List(M!.transitions,r->List(r,x->x^f)),M!.output);
    IsSphereMachine(N);
    if HasAddingElement(M) then
        x := InitialState(AddingElement(M));
	if x^f=x then COPYADDER@(N,M); fi;
    fi;
    return N;
end);

InstallMethod(\*, "(IMG) for a mapping and an FR machine",
        [IsMapping, IsFRMachine and IsFRMachineStdRep],
        function(f,M)
    local S, trans, out, i, pi, x, N;
    S := StateSet(M);
    if S<>Source(f) or S<>Range(f) then
        Error("\*: source, range and stateset must be the same\n");
    fi;
    pi := WreathRecursion(M);
    trans := [];
    out := [];
    
    for i in [1..Length(M!.output)] do
        x := pi(GeneratorsOfFRMachine(M)[i]^f);
        Add(trans,x[1]);
        Add(out,x[2]);
    od;
    N := FRMachineNC(FamilyObj(M),S,trans,out);
    IsSphereMachine(N);
    if HasAddingElement(M) then
        x := InitialState(AddingElement(M));
	if x^f=x then COPYADDER@(N,M); fi;
    fi;
    return N;
end);

InstallMethod(\^, "(IMG) for a group FR machine and a mapping",
        [IsGroupFRMachine, IsMapping],
        function(M,f)
    local S, newS, trans, out, i, pi, x, finv, N;
    S := StateSet(M);
    if S<>Source(f) then
        Error("\^: source and stateset must be the same\n");
    fi;
    newS := Range(f);
    pi := WreathRecursion(M);
    trans := [];
    out := [];
    finv := InverseGeneralMapping(f);
    if finv=fail then return fail; fi;
    for i in GeneratorsOfGroup(newS) do
        x := pi(ImagesRepresentative(finv,i));
        Add(trans,List(x[1],x->ImagesRepresentative(f,x)));
        Add(out,x[2]);
    od;
    N := FRMachineNC(FamilyObj(M),newS,trans,out);
    IsSphereMachine(N);
    if HasAddingElement(M) then
        SetAddingElement(N,FRElement(N,InitialState(AddingElement(M))^f));
    fi;
    return N;
end);

#############################################################################

#############################################################################
##
#A Polynomial machines
##
InstallMethod(ViewString, "(IMG) for a polynomial sphere machine",
        [IsPolynomialSphereMachine and IsFRMachineStdRep],
        M->CONCAT@FR("<sphere machine with alphabet ", AlphabetOfFRObject(M), " and adder ", AddingElement(M), " on ", StateSet(M), "/", RelatorsOfFpGroup(StateSet(M)),">"));

InstallMethod(DisplayString, "(IMG) for a polynomial sphere machine",
        [IsPolynomialSphereMachine and IsFRMachineStdRep],
        M->CONCAT@FR(DISPLAYFRMACHINE@FR(M),
                "Adding element: ",AddingElement(M),"\n",
                "Relators: ",RelatorsOfFpGroup(StateSet(M)),"\n"));

InstallMethod(AddingElement, "(IMG) for a sphere machine",
        [IsSphereMachine],
        function(M)
    local i, deg, cycle;

    deg := Length(AlphabetOfFRObject(M));

    for i in Reversed(GeneratorsOfFRMachine(M)) do
        cycle := Cycles(PermList(Output(M,i)),[1..deg]);
        if Length(cycle)=1 and IsConjugate(StateSet(M),Product(Transitions(M,i){cycle[1]}),i) then
            return FRElement(M,i);
        fi;
    od;
    TryNextMethod();
end);
#############################################################################

BindGlobal("REORDERREC@", function(m,perm)
    # reorder the entries of m=[trans,perm] according to the permutation perm.
    local i;
    for i in [1..Length(m[1])] do
        m[1][i] := Permuted(m[1][i],perm);
        m[2][i] := ListPerm(PermList(m[2][i])^perm,Length(m[2][i]));
    od;
end);

BindGlobal("FLIPSPIDER@", function(m,adder)
    # reverses the marking at infinity, by conjugating by the base change
    # <adder,1,...,1>[1,deg,deg-1,...,2]
    # this is a change of basis that sends the adding machine to its inverse,
    # assuming the adding machine was normalized to be <adder,1,...,1>(i->i-1)
    local j, k, deg;
    
    for j in [1..Length(m[1])] do
        m[1][j][1] := adder^-1*m[1][j][1];
        k := Position(m[2][j],1);
        m[1][j][k] := m[1][j][k]*adder;
    od;
    deg := Length(m[1][1]);
    REORDERREC@(m,PermList(Concatenation([1],[deg,deg-1..2])));
end);

InstallMethod(ComplexConjugate, "(IMG) for a sphere machine",
        [IsSphereMachine],
        function(M)
    local S, CS, N, a, m;
    S := StateSet(M);
    CS := SphereGroup(Reversed(OrderingOfSphereGroup(S)),ExponentsOfSphereGroup(S),FreeGroupOfFpGroup(S));
    N := M^GroupHomomorphismByImages(S,CS,GeneratorsOfGroup(S),List(GeneratorsOfGroup(CS),Inverse));
    if HasAddingElement(N) then
        m := [List(N!.transitions,ShallowCopy),ShallowCopy(N!.output)];
        a := InitialState(AddingElement(N));
        FLIPSPIDER@(m,a);
        a := Inverse(a);
        N!.transitions := List(m[1],x->List(x,y->y^a));
        N!.output := m[2];
        N!.AddingElement := FRElement(N,a); # make it a positive generator
    fi;
    return N;
end);

# keep track of adder
InstallMethod(ChangeFRMachineBasis, [IsPolynomialSphereMachine, IsCollection, IsPerm],
        function(m,c,p)
    local newm;
    newm := CHANGEFRMACHINEBASIS@FR(m,c,p);
    IsSphereMachine(newm);
    SetAddingElement(newm,FRElement(newm,InitialState(AddingElement(m))));
    return newm;
end);

InstallMethod(ChangeFRMachineBasis, [IsSphereMachine, IsCollection, IsPerm],
        function(m,c,p)
    local newm;
    newm := CHANGEFRMACHINEBASIS@FR(m,c,p);
    IsSphereMachine(newm);
    return newm;
end);

InstallMethod(TensorProductOp, "(IMG) for two polynomial machines",
        [IsList, IsPolynomialSphereMachine and IsFRMachineStdRep],
        function(M, N)
    local R;
    R := ApplicableMethod(TensorProductOp, [M,N], 0, 2)(M,N);
    if InitialState(AddingElement(M))=InitialState(AddingElement(M)) then
        SetAddingElement(R,FRElement(R,InitialState(AddingElement(M))));
    fi;
    return R;
end);

InstallMethod(RotatedSpider, "(IMG) for a polynomial IMG machine",
        [IsPolynomialSphereMachine,IsInt],
        function(M,p)
    local adder;
    adder := DecompositionOfFRElement(AddingElement(M)^p);
    return ChangeFRMachineBasis(M,List(adder[1],InitialState),PermList(adder[2]));
end);

InstallMethod(RotatedSpider, "(IMG) for a polynomial IMG machine",
        [IsPolynomialSphereMachine],
        M->RotatedSpider(M,1));

#############################################################################
BindGlobal("TRUNC@", function(a)
    # fractional part of rational, contained in [0,1)
    a := a-Int(a);
    if a<0 then return a+1; else return a; fi;
end);

BindGlobal("COMPOSERECURSION@", function(trans,out,pre,post)
    # twist the recursion [trans,out] by precomposing by pre^-1, and
    # post-composing by post. These are homomorphisms with range
    # the (source, range) group of [trans,out].
    # returns the new [trans,out].
    local i, w, newout, newtrans, deg, psi, source, gens;
    deg := Length(trans[1]);
    source := Source(pre);
    gens := GeneratorsOfGroup(source);

    newout := [];
    newtrans := [];

    if deg=1 then # bug with wreath products on the trivial group
        psi := GroupHomomorphismByImagesNC(source,Source(post),gens,List([1..Length(gens)],i->trans[i][1]));
        for i in GeneratorsOfGroup(Range(pre)) do
            Add(newout,[1]);
            Add(newtrans,[ImagesRepresentative(post,PreImagesRepresentative(pre,i)^psi)]);
        od;
        return [newtrans,newout];
    fi;

    w := WreathProduct(Source(post),SymmetricGroup(deg));
    psi := GroupHomomorphismByImagesNC(source,w,gens,List([1..Length(gens)],i->Product([1..deg],j->trans[i][j]^Embedding(w,j))*PermList(out[i])^Embedding(w,deg+1)));

    for i in GeneratorsOfGroup(Range(pre)) do
        w := PreImagesRepresentative(pre,i)^psi;
        Add(newout,ListPerm(w![deg+1],deg));
        Add(newtrans,List([1..deg],j->ImagesRepresentative(post,w![j])));
    od;
    
    return [newtrans,newout];
end);

BindGlobal("FJ@",["Fatou","Julia"]);

# this is a machine with no "infinity" element; the product of the generators (in some order)
# is assumed to be an adding element
InstallMethod(SupportingRays, "(IMG) for a group FR machine",
        [IsGroupFRMachine],
        function(M)
    local e;
    e := SPIDERALGORITHM@(M);
    if e=fail then
        return fail;
    elif e.minimal=false then
        return e;
    fi;
    return [Length(AlphabetOfFRObject(M)),e.supportingangles[1],e.supportingangles[2]];
end);

InstallMethod(SupportingRays, "(IMG) for a polynomial sphere machine",
        [IsPolynomialSphereMachine],
        function(M)
    local g, gens, img, adder, f, nf, trans, output, spider, newM;
    
    # put the adding element in standard form
    M := NormalizedPolynomialSphereMachine(M);
    g := StateSet(M);
    gens := GeneratorsOfGroup(g);
    img := ShallowCopy(gens);
    adder := Position(gens,InitialState(AddingElement(M)));
    Remove(img,adder);
    f := FreeGroup(Length(img));
    nf := GroupHomomorphismByImages(f,g,img);
    trans := List(M!.transitions,x->List(x,y->PreImagesRepresentative(nf,y)));
    output := ShallowCopy(M!.output);
    Remove(trans,adder);
    Remove(output,adder);
    newM := FRMachineNC(FamilyObj(M),f,trans,output);
    
    spider := SPIDERALGORITHM@(newM);
    if spider=fail then
        return fail;
    elif spider.minimal=false then
        return spider;
    fi;
    return [Length(AlphabetOfFRObject(M)),spider.supportingangles[1],spider.supportingangles[2]];
end);

InstallOtherMethod(SupportingRays, "(IMG) for a Mealy machine",
        [IsMealyMachine],
        function(M)
    return SupportingRays(AsPolynomialSphereMachine(M));
end);

InstallMethod(AsPolynomialSphereMachine, "(IMG) for a polynomial machine",
        [IsPolynomialSphereMachine],
        M->M);
#############################################################################

##############################################################################
##
#M Normalize polynomial machine
##
BindGlobal("REDUCEINNER@", function(img0,gen)
    # repeatedly apply conjugation by elements of gen to reduce img, the list
    # of images of generators under an endomorphism.
    # modifies img, and returns the conjugating element that was used to reduce;
    # i.e. after the run, List(img,x->x^elt) = img0

    local cost, oldcost, i, img, oldimg, elt, idle;

    elt := One(img0[1]);
    oldimg := img0;
    oldcost := Sum(oldimg,Length);
    gen := Concatenation(gen,List(gen,Inverse));
    
    repeat
        idle := true;
        for i in gen do
            img := List(oldimg,x->x^i);
            cost := Sum(img,Length);
            if cost < oldcost then
                oldcost := cost;
                oldimg := img;
                elt := elt*i;
                idle := false;
                break;
            fi;
        od;
    until idle;
    for i in [1..Length(oldimg)] do img0[i] := oldimg[i]; od;
    return elt;
end);

InstallMethod(NormalizedPolynomialSphereMachine, "(IMG) for a polynomial FR machine",
        [IsSphereMachine],
        function(M)
    # conjugate the recursion so that the adding element t becomes of the form (t,...,1)s,
    # where s is the cycle i|->i-1 mod d.
    # adder is an element of <model>.
    # model is the ambient fundamental group.
    local model, trans, out, adder, N, deg, perm, x, i, j, basis;
    
    model := StateSet(M);
    trans := List(M!.transitions,ShallowCopy);
    deg := Length(trans[1]);
    out := List(M!.output,ShallowCopy);
    adder := InitialState(AddingElement(M));
    
    while not ISADDER@(M,adder) do
        Error("Element ",adder," is not an adding element");
    od;
    
    perm := PermList(Concatenation([deg],[1..deg-1]));
    perm := RepresentativeAction(SymmetricGroup(deg),PermList(Output(M,adder)),perm);
    REORDERREC@([trans,out],perm);

    basis := [];
    x := One(model);
    for i in [deg,deg-1..1] do
        basis[i] := x;
        x := x*Transition(M,adder,i);
    od;
    basis := RepresentativeAction(model,adder,x)*basis;
    for i in [1..Length(trans)] do
        for j in [1..deg] do
            trans[i][j] := basis[j]*trans[i][j]/basis[out[i][j]];
        od;
    od;

    N := FRMachineNC(FamilyObj(M),M!.free,trans,out);
    COPYADDER@(N,M);
    IsSphereMachine(N);
    return N;
end);

InstallMethod(SimplifiedSphereMachine, "(IMG) for a polynomial sphere machine",
        [IsPolynomialSphereMachine],
        function(M)
    local r, i, x;
    Info(InfoIMG,1,"Simplification not yet implemented");
    r := SPIDERALGORITHM@(M);
    if r<>fail and r.minimal then
        SetCorrespondence(r.machine,r.transformation);
        return r.machine;
    fi;
    return M;
end);

InstallMethod(SimplifiedSphereMachine, "(IMG) for a sphere machine",
        [IsSphereMachine],
        function(M)
    local N;
    Info(InfoIMG,1,"Simplification not yet implemented for general sphere machines");
    return N;
end);
#############################################################################

#############################################################################
##
#M PolynomialMealyMachine
#M PolynomialSphereMachine
##
BindGlobal("KM_INPART@", function(part,x)
    # part is a list of subsets of [0,1) given as a list of pairs [lo,hi).
    # returns the index of the subset that contains x
    local i, j;
    for i in [1..Length(part)] do
        for j in part[i] do
            if j[1]<=x and x<j[2] then
                return i;
            fi;
        od;
    od;
end);

BindGlobal("KM_SPLITPART@", function(part,p,epsilon)
    # part is a list of subsets as above.
    # p is a list of points on [0,1).
    # intersect part with the partition cutting [0,1) at all p[i]-epsilon.
    local ind, i, j, newp;

    ind := List(p,x->KM_INPART@(part,x-epsilon));
    if Size(Set(ind))<>1 then
        Error("Some parts cross");
    fi;
    ind := part[ind[1]];
    i := 1;
    j := 1;
    while j<Length(p) do # split off ind\cap [p[j]..p[j+1]]
        while ind[i][2]<=p[j] do
            i := i+1;
        od;
        newp := [];
        if ind[i][1]<p[j] then
            Add(ind,[ind[i][1],p[j]],i);
            i := i+1;
            ind[i][1] := p[j];
        fi;
        j := j+1;
        while Length(ind)>=i and ind[i][2]<=p[j] do
            Add(newp,Remove(ind,i));
        od;
        if Length(ind)>= i and ind[i][1]<p[j] then
            Add(newp,[ind[i][1],p[j]]);
            ind[i][1] := p[j];
        fi;
        Add(part,newp);
    od;
end);

BindGlobal("PCPORDERS@", function(C,d,pcp)
    local pcporders, lcm, liftorder, i, j, lift, pos, idle, upperbd;

    pcporders := List(pcp,x->1); Add(pcporders,0);

    upperbd := Product(List(C,r->Length(r[1])))+1; # maximal order of finite points
    repeat
        idle := true;
        for i in [1..Length(pcp)] do
            lcm := 1;
            for j in [0..d-1] do
                liftorder := 1;
                lift := [(pcp[i][1]+j)/d,pcp[i][2]];
                pos := Position(pcp,lift);
                if pos<>fail then
                    liftorder := liftorder*pcporders[pos];
                fi;
                j := First(C,r->lift[1] in r[1] and lift[2]=r[2]);
                if j<>fail then
                    liftorder := liftorder*Length(j[1]);
                fi;
                lcm := Lcm(lcm,liftorder);
            od;
            if lcm<>pcporders[i] and pcporders[i]<upperbd then
                pcporders[i] := lcm;
                idle := false;
            fi;
        od;
    until idle;

    # all points that got bigger order are infinite

    for i in [1..Length(pcp)] do
        if pcporders[i]>=upperbd then pcporders[i] := 0; fi;
    od;
    return pcporders;
end);

InstallMethod(PolynomialSphereMachine, "(IMG) for a degree, Fatou and Julia preangles, and options",
        [IsPosInt,IsList,IsList,IsRecord],
        function(d,F,J,machtype)
    # d is the degree
    # F is a list of Fatou critical points
    # J is a list of Julia critical points
    # machtype.mealy: mealy machine, default=f
    # machtype.formal: formal construction (non-coalesced points with same itinerary), default=f
    # machtype.adding: formal construction, adding machine nice
    # machtype.orbispace: give points their minimial degree, default=f
    # returns an FR machine, and sets correspondence to [fF,fJ], where
    # these functions return, for a given angle, the corresponding generator.
    local C, V, i, j, part, pcp, f, newf, gens, trans, t, out, o, p, q, rank, epsilon, epi,
          one, machine, orders;

    for i in RecNames(machtype) do
        while not i in ["formal","adding","orbispace","mealy"] do
            Error("PolynomialSphereMachine: option '",i,"' not recognized");
        od;
        while not IsBool(machtype.(i)) do
            Error("PolynomialSphereMachine: option '",i,"' should be a boolean");
        od;
    od;
    machtype := ShallowCopy(machtype);
    if not IsBound(machtype.mealy) then machtype.mealy := false; fi;
    if not IsBound(machtype.formal) then machtype.formal := false; fi;
    if not IsBound(machtype.adding) then machtype.adding := true; fi;
    if not IsBound(machtype.orbispace) then machtype.orbispace := false; fi;

    C := Concatenation(List(F,x->[x,FJ@[1]]),List(J,x->[x,FJ@[2]]));
    for i in C do
        if IsRat(i[1]) then
            i[1] := Set([1..d],j->TRUNC@((i[1]+j)/d));
        else
            i[1] := Set(i[1],TRUNC@);
        fi;
    od;
    epsilon := 1/2/Lcm(List(C,x->Lcm(List(x[1],DenominatorRat))));
    
    if ForAny(C,x->Size(Set(d*x[1],TRUNC@))<>1) then
        Error("F and J must be ",d,"-prearguments");
    fi;

    pcp := [];
    V := List(C,i->[TRUNC@(d*i[1][1]),i[2]]); # critical values
    for i in V do
        p := i;
        while not p in pcp do
            AddSet(pcp,p);
            p := [TRUNC@(d*p[1]),p[2]];
        od;
        while Gcd(DenominatorRat(p[1]),d)<>1 do
            p := [TRUNC@(d*p[1]),p[2]];
        od;
        q := p;
        t := false;
        repeat
            p := [TRUNC@(d*p[1]),p[2]];
            t := t or (p in V); # check if there's a critical point on cycle
        until q=p;
        while (p[2]=FJ@[1] and not t) or (p[2]=FJ@[2] and t and d=2) do
            Error("critical value ",i[1]," should not be in the ",p[2]," set");
        od;
    od;
    Sort(pcp,function(x,y)
        return x[1]<y[1] or (x[1]=y[1] and x[2]=FJ@[2] and y[2]=FJ@[1]);
    end);
    rank := Length(pcp);

    part := [[[0,1]]];
    for i in C do
        if i[2]=FJ@[1] then
            KM_SPLITPART@(part,i[1],0);
        else
            KM_SPLITPART@(part,i[1],epsilon);
        fi;
    od;
        
    while Sum(C,x->Length(x[1])-1)<>d-1 do
        Error("F and J describe a map of wrong degree");
    od;
    
    part := part{List([0..d-1]/d,x->KM_INPART@(part,x))};

    if machtype.mealy then
        gens := [1..rank+1];
        one := rank+1;
    else
        if machtype.orbispace then
            orders := PCPORDERS@(C,d,pcp);
        else
            orders := List([0..rank],i->0);
        fi;
        f := SphereGroup(rank+1,orders);
        gens := GeneratorsOfGroup(f);
        one := One(f);
    fi;
    trans := [];
    out := [];
    for i in pcp do
        t := [];
        o := [];
        for j in [0..d-1] do
            q := (i[1]+j)/d;
            p := First(C,x->i[2]=x[2] and q in x[1]);
            if p=fail then
                p := q;
            else
                p := p[1][Position(p[1],q) mod Length(p[1])+1];
            fi;
            if machtype.mealy or not machtype.adding then
                o[KM_INPART@(part,q)] := KM_INPART@(part,p);

                p := Position(pcp,[q,i[2]]);
                if p<>fail then
                    p := gens[p];
                else
                    p := one;
                fi;
                t[KM_INPART@(part,q)] := p;
            else
                if p<>q then
                    Add(o,1+(j+d*(p-q)) mod d);
                    if p>q then
                        Add(t,Product(gens{Filtered([rank,rank-1..1],j->pcp[j][1]>q and pcp[j][1]<p and pcp[j][2]=i[2])},one)^-1);
                    else
                        Add(t,Product(gens{Filtered([rank,rank-1..1],j->pcp[j][1]>=p and pcp[j][1]<=q and pcp[j][2]=i[2])},one));
                    fi;
                else
                    Add(o,j+1);
                    p := Position(pcp,[q,i[2]]);
                    if p<>fail then
                        Add(t,gens[p]);
                    else
                        Add(t,one);
                    fi;
                fi;
            fi;
        od;
        Add(trans,t);
        Add(out,o);
    od;
    if machtype.mealy then
        Add(trans,List([1..d],i->Length(gens)));
        Add(out,[1..d]);
        machine := MealyMachine(trans,out);
        SetCorrespondence(machine,pcp);
        SetAddingElement(machine,Product(Reversed(gens),i->FRElement(machine,i)));
    else
        t := [gens[Length(gens)]];
        Append(t,List([1..d-1],i->one));
        Add(trans,t);
        Add(out,Concatenation([d],[1..d-1]));
        machine := FRMachine(f,trans,out);
        SetAddingElement(machine,FRElement(machine,gens[Length(gens)]));
        SetCorrespondence(machine,pcp);
        IsSphereMachine(machine);
    fi;

    if not (machtype.formal or machtype.mealy) then
        p := [[1..rank]]; # start by trivial partition
        t := List(pcp,x->x[1]);
        for i in [1..rank] do # repeat this many times
            q := [];
            for p in p do # go through all parts
                o := List([1..d],i->[]);
                for i in p do
                    if ForAny(C,r->t[i] in r[1]) then # it's critical, separate it
                        Add(o,[i]);
                    else
                        Add(o[KM_INPART@(part,t[i])],i); # see its address
                    fi;
                od;
                UniteSet(q,o);
            od;
            t := List(t,x->TRUNC@(d*x)); # mop the post-critical points by the doubling map
            p := Difference(q,[[]]);
        od;
        if Length(p)<rank then
            q := [];
            for i in p do
                o := one;
                for j in [i[Length(i)],i[Length(i)]-1..i[1]] do
                    if not j in i then o := gens[j]^-1*o; fi;
                    o := o*gens[j];
                od;
                Add(q,o);
            od;
            pcp := List(p,r->[pcp{r}[1],pcp[r[1]][2]]);
            orders := orders{List(p,Representative)}; Add(orders,0);
            newf := SphereGroup(Length(p)+1,orders);
            epi := GroupGeneralMappingByImages(f,newf,Concatenation(q,[gens[Length(gens)]]),GeneratorsOfGroup(newf));
            trans := COMPOSERECURSION@(trans,out,epi,epi);
            out := trans[2]; trans := trans[1];
            machine := FRMachine(newf,trans,out);
            SetAddingElement(machine,FRElement(machine,GeneratorsOfGroup(newf)[Length(orders)]));
            SetCorrespondence(machine,pcp);
        fi;
    fi;

    return machine;
end);

InstallMethod(PolynomialMealyMachine, "(IMG) for a degree, Fatou and Julia preangles",
        [IsPosInt,IsList,IsList],
        function(n,F,J)
    return PolynomialSphereMachine(n,F,J,rec(mealy:=true));
end);

InstallMethod(PolynomialSphereMachine, "(IMG) for a degree, Fatou and Julia preangles",
        [IsPosInt,IsList,IsList],
        function(n,F,J)
    return PolynomialSphereMachine(n,F,J,rec());
end);

BindGlobal("FATOUANGLES@", function(n,A)
    return Filtered(A,x->Gcd(DenominatorRat(x),n)=1);
end);
BindGlobal("JULIAANGLES@", function(n,A)
    return Filtered(A,x->Gcd(DenominatorRat(x),n)<>1);
end);

InstallMethod(PolynomialMealyMachine, "(IMG) for a degree and preangles",
        [IsPosInt,IsList],
        function(n,A)
    return PolynomialMealyMachine(n,FATOUANGLES@(n,A),JULIAANGLES@(n,A));
end);

InstallMethod(PolynomialSphereMachine, "(IMG) for a degree and preangles",
        [IsPosInt,IsList],
        function(n,A)
    return PolynomialSphereMachine(n,FATOUANGLES@(n,A),JULIAANGLES@(n,A));
end);

InstallMethod(PolynomialSphereMachine, "(IMG) for a degree, preangles, and bool",
        [IsPosInt,IsList,IsRecord],
        function(n,A,nice)
    return PolynomialSphereMachine(n,FATOUANGLES@(n,A),JULIAANGLES@(n,A),nice);
end);

InstallMethod(Mating, "(IMG) for two polynomial sphere machines and a boolean",
        [IsPolynomialSphereMachine,IsPolynomialSphereMachine,IsBool],
        function(M1,M2,formal)
    local machines, adders, deg, states, amalgam, i, s, trans, out, embed, sum;

    machines := [M1,M2];
    adders := List(machines,x->InitialState(AddingElement(x)));
    
    deg := List(machines,m->Length(AlphabetOfFRObject(m)));
    while deg[1]<>deg[2] do
        Error("In a mating, the machines must have same degree, not ",deg);
    od;
    deg := deg[1];
           
    machines := List(machines,NormalizedPolynomialSphereMachine);
    machines[2] := ChangeFRMachineBasis(machines[2],PermList([deg,deg-1..1])); # make it inverse

    states := List(machines,StateSet);
    amalgam := AmalgamatedFreeProduct(states[1],states[2],adders[1],adders[2]);
    embed := EmbeddingsOfAmalgamatedFreeProduct(amalgam);

    trans := [];
    out := [];
    for i in [1..2] do
        for s in GeneratorsOfGroup(states[i]) do
            if s=adders[i] then continue; fi;
            Add(trans,List(AlphabetOfFRObject(machines[i]),a->Transition(machines[i],s,a)^embed[i]));
            Add(out,Output(machines[i],s));
            Assert(0,Position(GeneratorsOfGroup(amalgam),s^embed[i])=Length(out));
        od;
    od;
    sum := FRMachineNC(FamilyObj(machines[1]),amalgam,trans,out);
    IsSphereMachine(sum);
    
    SetEquatorElement(sum,adders[1]^embed[1]);
    SetEquatorTwist(sum,GroupHomomorphismByImages(amalgam,amalgam,
            Concatenation(GeneratorsOfGroup(amalgam){[1..RankOfSphereGroup(states[1])-1]},List(GeneratorsOfGroup(amalgam){[RankOfSphereGroup(states[1])..RankOfSphereGroup(amalgam)]},x->x^EquatorElement(sum)))));
    
    if not formal then
        Error("Non-formal matings are not yet implemented. Complain to laurent.bartholdi@gmail.com");
    fi;
    
    SetCorrespondence(sum,embed);
    return sum;
end);
InstallMethod(Mating, "(IMG) for two polynomial sphere machines",
        [IsPolynomialSphereMachine,IsPolynomialSphereMachine],
        function(M1,M2)
    return Mating(M1,M2,true);
end);
#############################################################################

#############################################################################
##
#F Automorphisms of machines and virtual endomorphisms
##
BindGlobal("PUREMCG@", function(arg)
    local r, i, j, m, ni, nj, gens, img, aut, relator, G, maker;
    
    if Length(arg)=1 and IsFpGroup(arg[1]) then
        G := arg[1];
        relator := RelatorsOfFpGroup(G)[1];
        maker := w->ElementOfFpGroup(FamilyObj(One(G)),w);
    elif Length(arg)=2 and IsFreeGroup(arg[1]) and arg[2] in arg[1] then
        G := arg[1];
        relator := arg[2];
        maker := w->w;
    else
        Error("PUREMCG@FR: requires a fp group or a free group and a relator");
    fi;
    gens := GeneratorsOfGroup(G);
        
    aut := [];
    for ni in [1..Length(relator)] do
        for nj in [ni+1..Length(relator)] do
            m := maker(Subword(relator,ni+1,nj-1));
            i := LetterRepAssocWord(relator)[ni];
            j := LetterRepAssocWord(relator)[nj];
            img := ShallowCopy(gens);
            img[i] := gens[i]^(m*gens[j]/m);
            img[j] := gens[j]^(m^-1*gens[i]*m*gens[j]);
            Add(aut,GroupHomomorphismByImages(G,G,gens,img));
        od;
    od;
    if Length(gens)=2 then Add(aut,InnerAutomorphism(G,gens[1])); fi;
    return Group(aut);
end);

InstallMethod(AutomorphismVirtualEndomorphism, "(IMG) for a virtual endo",
        [IsGroupHomomorphism],
        function(vendo)
    # given a virtual endomorphism h->g, constructs the induced
    # virtual endomorphism on aut(g)
    local g, h, freeg, freeh, isog, isoh, embedding, v, vinverse,
          mcg, mch, freemcg, cch, restricth, mcv, mcvimg, act;
    
    g := Range(vendo);
    h := Source(vendo);
    isog := IsomorphismSimplifiedFpGroup(g);
    isog := isog*GroupHomomorphismByImages(Range(isog),FreeGroupOfFpGroup(Range(isog)));
    isoh := IsomorphismFpGroup(h);
    isoh := isoh*IsomorphismSimplifiedFpGroup(Range(isoh));
    isoh := isoh*GroupHomomorphismByImages(Range(isoh),FreeGroupOfFpGroup(Range(isoh)));
    freeg := Range(isog);
    freeh := Range(isoh);
    embedding := GroupHomomorphismByImages(freeh,freeg,List(GeneratorsOfGroup(h),x->x^isoh),List(GeneratorsOfGroup(h),x->x^isog));
    v := GroupHomomorphismByImages(freeh,freeg,List(GeneratorsOfGroup(h),x->x^isoh),List(GeneratorsOfGroup(h),x->(x^vendo)^isog));
    vinverse := GroupHomomorphismByImages(freeg,freeh,GeneratorsOfGroup(freeg),List(GeneratorsOfGroup(freeg),x->PreImagesRepresentative(vendo,PreImage(isog,x))^isoh));
    
    mcg := PUREMCG@(g);
    freemcg := Group(List(GeneratorsOfGroup(mcg),x->x^isog));
#    isomcg := GroupHomomorphismByImagesNC(mcg,freemcg);
    
    # first the subgroup of mcg stabilizing h
    act := function(subgroup,aut)
        return Group(List(GeneratorsOfGroup(subgroup),x->x^aut));
    end;
    mch := Stabilizer(freemcg,Image(embedding),act);

    # then the subgroup fixing parabolic conjugacy classes in h;
    cch := List(GeneratorsOfGroup(h),x->ConjugacyClass(freeh,x^isoh));
    restricth := aut->GroupHomomorphismByImages(freeh,freeh,List(GeneratorsOfGroup(freeh),g->PreImagesRepresentative(embedding,(g^embedding)^aut)));
    act := function(list,aut)
        return List(list,x->ConjugacyClass(freeh,Representative(x)^restricth(aut)));
    end;
    mch := Stabilizer(mch,cch,act);
    
    mch := GroupByGenerators(List(GeneratorsOfGroup(mch),a->GroupHomomorphismByImages(g,g,List(GeneratorsOfGroup(g),x->PreImage(isog,(x^isog)^a)))));
    return GroupHomomorphismByFunction(mch,mcg,a->GroupHomomorphismByImages(g,g,List(GeneratorsOfGroup(g),x->PreImage(isog,((PreImage(isoh,(x^isog)^vinverse)^a)^isoh)^v))));        
end);

BindGlobal("DISTILLATE@", function(dict,M)
    # distillates the machine M to a machine with at most one non-trivial
    # transition on each cycle; which is at the beginning of the cycle
    # and is a generator.
    # dict is a dictionary storing distillated machines, and a
    # standard representative for each distillation.
    # returns [standard machine,free group automorphism]
    # such that M=standard*automorphism
    local G, perm, trans, cc, zero, g, i, j, c, cg, aut;
    
    G := StateSet(M);
    cc := List(GeneratorsOfFRMachine(M),g->ConjugacyClass(G,g));
    zero := ConjugacyClass(G,One(G));
    aut := [];
    
    perm := M!.output;
    trans := List(M!.transitions,ShallowCopy);
    for i in [1..Length(trans)] do
        for c in Cycles(PermList(perm[i]),AlphabetOfFRObject(M)) do
            g := One(G);
            for j in c do
                g := g*trans[i][j];
                Unbind(trans[i][j]);
            od;
            cg := ConjugacyClass(G,g);
            if cg<>zero then
                j := Position(cc,cg);
                trans[i][c[1]] := cg;
                aut[j] := g;
            fi;
        od;
    od;
    c := LookupDictionary(dict,trans);
    if c=fail then
        AddDictionary(dict,trans,[M,aut]);
        return [M,IdentityMapping(G)];
    else
#j := GroupHomomorphismByImages(G,G,c[2],aut); if j=fail then
#    Error("1"); fi; return [c[1],j];
        return [c[1],GroupHomomorphismByImages(G,G,c[2],aut)];
    fi;
end);

BindGlobal("MOTIONGROUP@", function(G)
    local i, j, gens, img, aut;
    
    gens := GeneratorsOfGroup(G);
        
    aut := [];
    for i in [1..Length(gens)] do
        for j in [1..Length(gens)] do
            if i=j then continue; fi;
            img := ShallowCopy(gens);
            img[i] := gens[i]^gens[j];
            Add(aut,GroupHomomorphismByImages(G,G,img));
        od;
    od;
    return Group(aut);
end);

BindGlobal("MATCHMARKINGS@", function(M1,M2)
    # return a homomorphism phi from StateSet(M) to <group> such that:
    # if g[i]^k lifts to a conjugate h of g[j] for some integer k, then
    # phi(g[i]^k) = corresponding expression obtained from <M2>.
    local src, dst, transM, transR, i, c, x, gens, g, h, group;

    group := FamilyObj(Transition(M2,1,1))!.group;
    gens := GeneratorsOfGroup(StateSet(M1));
    transM := [One(StateSet(M1))];
    transR := [One(group)];
    src := [1];
    dst := Difference(AlphabetOfFRObject(M1),src);
    while dst<>[] do
        c := Cartesian(src,dst);
        for i in [1..Length(gens)] do
            x := First(c,p->Output(M1,gens[i])[p[1]]=p[2]);
            if x<>fail then break; fi;
        od;
        transM[x[2]] := Transition(M1,gens[i],x[1])^-1*transM[x[1]];
        transR[x[2]] := Transition(M2,gens[i],x[1])^-1*transR[x[1]];
        Add(src,Remove(dst,Position(dst,x[2])));
    od;

    src := [];
    dst := [];
    for i in [1..Length(gens)] do
        Assert(0,Output(M1,i)=Output(M2,i));
        for c in Cycles(PermList(Output(M1,i)),AlphabetOfFRObject(M1)) do
            g := Product(Transitions(M1,i){c})^transM[c[1]];
            h := Product(Transitions(M2,i){c})^transR[c[1]];
            Add(src,g);
            Add(dst,h);
        od;
    od;

    x := GroupHomomorphismByImagesNC(StateSet(M1),group,src,dst);
    dst := List(gens,g->g^x);
    REDUCEINNER@(dst,GeneratorsOfMonoid(group));

    return GroupHomomorphismByImagesNC(StateSet(M1),group,gens,dst);
end);

InstallMethod(AutomorphismSphereMachine, "(IMG) for an IMG machine",
        [IsSphereMachine],
        function(M)
    local orbit, oorbits, act, inneract, pmcg, states,
          output, transition, o, t,
          iso, epi, isoepi, a, b, c, d, g, newM, distillations;
    
    states := StateSet(M);
    distillations := NewDictionary([],true);
    act := function(machine,aut) return DISTILLATE@(distillations,aut^-1*machine)[1]; end;
    inneract := function(machine,g) return DISTILLATE@(distillations,InnerAutomorphism(states,g^-1)*machine)[1]; end;
    pmcg := AutomorphismGroup(states);
    orbit := Orbit(pmcg,DISTILLATE@(distillations,M)[1],act);

    # sort orbit so that consecutive blocks are related by inner automorphisms
    oorbits := OrbitsDomain(states,orbit,inneract);
    orbit := Concatenation(oorbits);
    
    # now oorbits[i][j] is a distilled machine;
    # oorbits[i][j] and oorbits[i][k] are related by inner automorphisms

if false then
    epi := IsomorphismFpGroup(pmcg); # really to a free group (!!!)
    
    # todo:
    # 1) fix bugs in MATCHMARKINGS that delete inner auts
    # 2) construct machine for quotient, on smaller alphabet; construct embedding of original machine too
    # 3) find irredundant gens of pure sphere braid group? i.e. giving correct abelianization rank?
    
    output := [];
    transition := [];
    
    for g in GeneratorsOfGroup(Range(epi)) do
        g := PreImagesRepresentative(epi,g)^-1;
        o := [];
        t := [];
        for a in [1..Length(orbit)] do
            newM := g*orbit[a];
            b := Position(orbit,DISTILLATE@(distillations,newM)[1]);
            Add(o,b);
            Add(t,ImagesRepresentative(epi,MATCHMARKINGS@(newM,orbit[b])));
        od;
        Add(output,o);
        Add(transition,t);
    od;
else
    iso := IsomorphismFpGroup(pmcg); # really to a free group
    epi := EpimorphismToOut(pmcg); # really to a free group
    isoepi := GroupHomomorphismByImages(Range(iso),Range(epi),List(GeneratorsOfGroup(pmcg),x->x^iso),List(GeneratorsOfGroup(pmcg),x->x^epi));

    output := [];
    transition := [];
    
    for g in GeneratorsOfGroup(Range(epi)) do
        o := [];
        t := [];
        for a in [1..Length(oorbits)] do
            newM := PreImagesRepresentative(epi,g)^-1*oorbits[a][1];
            d := DISTILLATE@(distillations,newM);
            b := PositionProperty(oorbits,o->d[1] in o);
            c := RepresentativeAction(states,d[1],oorbits[b][1],inneract);
            newM := InnerAutomorphism(states,c^-1)*newM;
            Add(o,b);
            Add(t,ImagesRepresentative(epi,MATCHMARKINGS@(newM,oorbits[b][1])));
        od;
        Add(output,o);
        Add(transition,t);
    od;
fi;

    a := FRMachine(Range(epi),transition,output);
    SetCorrespondence(a,[oorbits,epi]);
    return a;
end);
#############################################################################

#############################################################################
# conversions
#
BindGlobal("SEPARATION@", function(set,s)
    local a, b, c, i, j, cut;
  
    cut := function(set,s)
        local bins, i;
        bins := List([0..Length(s)],i->[]);
        for i in set do
            Add(bins[PositionSorted(s,i)],i);
        od;  
        if Length(bins[1])=1 then
            UniteSet(bins[1],Remove(bins));
        fi;
        return bins;
    end;
  
    b := [];
    c := [];
    for i in set do
        a :=[];
        for j in s do
            UniteSet(a,Difference(Set(cut(j,i)),[[]]));
        od;
        s := a;
    od;  
    
    return s;
end);

BindGlobal("EXTERNALANGLESRELATION@",[]);

InstallGlobalFunction(ExternalAnglesRelation, function(degree,n)
    local a, b, c, i, j, set;
    
    if not IsBound(EXTERNALANGLESRELATION@[degree]) then
        EXTERNALANGLESRELATION@[degree] := [];
    fi;
    
    if not IsBound(EXTERNALANGLESRELATION@[degree][n]) then
        a := [0..degree-2]/(degree-1);
        set := [ShallowCopy(a)];
        for i in [2..n] do
            b := [1..degree^i-2]/(degree^i-1);
            SubtractSet(b,a);
            UniteSet(a,b);
            c := SEPARATION@(set,[b]);
            UniteSet(set,c);
        od;
        EXTERNALANGLESRELATION@[degree][n] := EquivalenceRelationByPartition(Rationals,set);
    fi;
    return EXTERNALANGLESRELATION@[degree][n];
end);

InstallGlobalFunction(ExternalAngle, function(arg)
    # convert supporting rays to external angle
    local deg, sr, n;
    
    if Length(arg)=1 and IsFRMachine(arg[1]) then
        n := SupportingRays(arg[1]);
    elif Length(arg)>=2 and IsPosInt(arg[1]) and IsList(arg[2]) then
        n := arg;
    elif Length(arg)=1 and IsList(arg[1]) then
        n := arg[1];
    else
        Error("ExternalAngle: call with FR machine or supporting rays data, not ",arg);
    fi;
    deg := n[1];
    if Length(n[2])=1 then
        sr := n[2][1];
    elif Length(n[3])=1 then
        sr := n[3][1];
    fi;
    sr := sr[1]*deg; sr := sr - Int(sr);

    if RemInt(DenominatorRat(sr),deg)=0 then
        return sr;
    fi;
    n := 1; while not IsInt((deg^n-1)*sr) do n := n+1; od;
    return EquivalenceClassOfElement(ExternalAnglesRelation(deg,n),sr);
end);

InstallMethod(KneadingSequence, "(IMG) for a rational angle",
        [IsRat],
        function(angle)
    local i, t, s, set, j, marked;
        
    s := [];
    set := [];
    t := angle;
    marked := angle>=1/7 and angle<=2/7 and ValueOption("marked")<>fail;
    
    if marked then j := "A1"; else j := 1; fi;
    while not t in s do
        Add(set,j);
        if marked then j := "C0"; else j := 0; fi;
        Add(s,t);
   
        t := 2*t - Int(2*t);
        if t>angle/2 and t<angle/2+1/2 then
            if marked then j := "C1"; else j := 1; fi;
        fi; 
        if marked and t>=1/7 and t<2/7 then j := "A1"; fi;
        if marked and t>=2/7 and t<4/7 then j := "B1"; fi;
        if marked and t>=9/14 and t<11/14 then j := "A0"; fi;
        if marked and t>=11/14 or t<1/14 then j := "B0"; fi;
    od;
    i := Position(s,t);

    if i=1 then 
        Remove(set);
    else
        s:=[[],set];
        for j in [1..i-1] do
            Add(s[1],set[1]);
            Remove(set,1);
        od; 
        set := s;
    fi;
  
    return set;
end);

InstallGlobalFunction(AllInternalAddresses, function(n)
    local a, b, c, i, j, set, s, external, isseparate;
    
    external := function(n)
        local a, b, c, i, j, set, s;
        a := [];
        set := [];
        s := [[]];
        for i in [2..n] do
            b := [1..2^i-2]/(2^i-1);
            SubtractSet(b,a);
            UniteSet(a,b);
            c := SEPARATION@(set,[b]);
            UniteSet(set,c);
            Add(s,c);
        od;
        return [a,s];
    end;

    isseparate := function(a,b)
        return not (a[1]<b[1] or a[2]>b[2]);
    end;

    a := external(n);
    s := a[2];
    a := a[1]; 
    for i in [2..n] do
        for j in [1..Length(s[i])]  do
            Add(s[i][j],i);
        od;   
    od;
    
    for i in [2..n] do
        for j in [1..Length(s[i])] do
            a:=i-1;
            while a>1 and Length(s[i][j])<4 do
                for b in s[a] do
                    if isseparate(s[i][j],b) then
                        for c in b do
                            Add(s[i][j],c);
                        od; 
                    fi;
                od; 
                a :=a-1;  
            od;
        od;   
    od;

    return s;
end);
#############################################################################

#E machine.gi . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
