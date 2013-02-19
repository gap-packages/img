#############################################################################
##
#W sphere.gi                                                Laurent Bartholdi
##
#Y Copyright (C) 2013, Laurent Bartholdi
##
#############################################################################
##
##  Sphere groups
##
#############################################################################

BindGlobal("MAKENFREL@", function(rel,degree)
    local i, len;
    len := Length(rel);
    if rel[1]<0 then rel := -Reversed(rel); fi;
    rel := [rel,[],[]];
    rel[2]{rel[1]} := [1..len];
    rel[3]{rel[1]} := [3*len,3*len-1..2*len+1];
    Append(rel[1],rel[1]);
    Append(rel[1],-Reversed(rel[1]));
    for i in [1..4*len] do
        if degree[AbsoluteValue(rel[1][i])]=2 then
            rel[1][i] := AbsoluteValue(rel[1][i]);
        fi;
    od;
    return rel;
end);

InstallMethod(ElementOfSphereGroup, "(IMG) for a sphere group element",
        [IsFamily,IsAssocWordWithInverse],
        function(fam,w)
    return ElementOfFpGroup(fam,FpElementNFFunction(fam)(w));
end);

BindGlobal("PREREDUCESPHERERWS@", function (kbrws,n)
    local p, i;
    if not kbrws!.reduced then
        ReduceRules(kbrws);
    fi;
    p := 1;
    while p <= Length(kbrws!.pairs2check) do
        i := kbrws!.pairs2check[p];
        if ForAll(kbrws!.tzrules[i[1]],w->Length(w)<=n) and ForAll(kbrws!.tzrules[i[2]],w->Length(w)<=n) then
            p := KBOverlaps( i[1], i[2], kbrws, p );
        fi;
        p := p+1;
    od;
    kbrws!.pairs2check := [];
end);

BindGlobal("SPHERENFFUNCTION@", function(g,ordering,power)
    local n, mon, id, rws, rules, freemon, freegroupfam, g2m, m2g, gens, weight, special;
    
    n := Length(power);
    freemon := FreeMonoid(2*n);
    gens := GeneratorsOfMonoid(freemon);
    m2g := Concatenation([1..n],[-1,-2..-n]);
    g2m := []; g2m{n+1+m2g} := [1..2*n];
    special := false;
    
    if AsSortedList(power)=[2,3,6] then # force order [3,2,6]
        id := RepresentativeAction(SymmetricGroup(3),power,[3,2,6],Permuted);
        id := id*id^(1,4)(2,5)(3,6);
        if id^3<>() then id := id*(1,4)(2,5)(3,6); fi;
        ordering := ShortLexOrdering(freemon,Permuted(gens,id));
    elif AsSortedList(power)=[2,4,4] then # force different +/- ordering
        ordering := ShortLexOrdering(freemon,Permuted(gens,(1,6)));
    elif Set(power)=[2] then # force non-consecutive indices
        if IsEvenInt(n) then
            ordering := ShortLexOrdering(freemon,gens{Concatenation(ordering{Concatenation([1,3..n-1],[2,4..n])},[n+1..2*n])});
        else
            weight := Concatenation(ListWithIdenticalEntries(n,1),ListWithIdenticalEntries(n,n));
            weight[ordering[n]] := n;
            special := ShortLexOrdering(freemon);
            ordering := WeightLexOrdering(freemon,gens{Concatenation(ordering{Concatenation([1,3..n],[2,4..n-1])},[n+1..2*n])},weight);
        fi;
    else
        ordering := ShortLexOrdering(freemon);
    fi;
    
    id := One(freemon);
    rules := List(RelatorsOfFpGroup(g),w->[AssocWordByLetterRep(FamilyObj(id),g2m{n+1+LetterRepAssocWord(w)}),id]);
    Append(rules,List([1..2*n],i->[AssocWordByLetterRep(FamilyObj(id),[i,g2m[n+1-m2g[i]]]),id]));
    mon := freemon / rules;
    rws := ReducedConfluentRewritingSystem(mon,ordering);
    if special<>false then
        special := KnuthBendixRewritingSystem(mon,special);
        PREREDUCESPHERERWS@(special,n);
        special := special!.tzrules;
    fi;
    rules := rws!.tzrules;

    freegroupfam := FamilyObj(UnderlyingElement(One(g)));
    
    if special=false then
        return w->AssocWordByLetterRep(freegroupfam,m2g{ReduceLetterRepWordsRewSys(rules,g2m{n+1+LetterRepAssocWord(w)})});
    else
        return w->AssocWordByLetterRep(freegroupfam,m2g{ReduceLetterRepWordsRewSys(special,ReduceLetterRepWordsRewSys(rules,g2m{n+1+LetterRepAssocWord(w)}))});
    fi;
end);

InstallGlobalFunction(SphereGroup, function (arg)
    local  F, G, fam, rel, n, ordering, power;
    
    while not (Length(arg) in [1,2] and (IsPosInt(arg[1]) or IsList(arg[1])) and (Length(arg)=1 or IsList(arg[2]))) do
        Error("Usage: SphereGroup(num_points or ordering [, degrees]");
    od;
    if IsPosInt(arg[1]) then
        n := arg[1];
        ordering := [n,n-1..1];
    else
        ordering := arg[1];
        n := Length(ordering);
        while Set(ordering)<>[1..n] do
            Error("<ordering> must contain once every element of ",[1..n]);
        od;
    fi;
    if Length(arg)=1 then
        power := List([1..n],x->0);
    else
        power := arg[2];
        while not (Length(power)=n and ForAll(power,x->x>=0)) do
            Error("<degrees> must be a list of non-negative integers");
        od;
    fi;    
    
    # create the underlying f.p. group
    F := FreeGroup(n);
    rel := List(Filtered([1..n],i->power[i]<>0),i->F.(i)^power[i]);
    Add(rel, Product(GeneratorsOfGroup(F){ordering},One(F)));
    G := F / rel;
    
    # adjust the family to encode better normal form calculation
    fam := FamilyObj(One(G));
    SetIsElementOfSphereGroupFamily(fam,true);
    fam!.defaultType := NewType(fam, IsElementOfSphereGroup and IsPackedElementDefaultRep);
    fam!.reduce := true;
    fam!.group := G;
    fam!.ordering := ordering;
    fam!.power := power;
    
    #rel := MAKENFREL@(ordering,power);
    #SetFpElementNFFunction(fam,w->AssocWordByLetterRep(FamilyObj(w),NFFUNCTION_FR(rel,power,true,LetterRepAssocWord(w))));
    SetFpElementNFFunction(fam,SPHERENFFUNCTION@(G,ordering,power));

    # reset some of the attributes of G that were erroneously computed
    G!.OneImmutable := ElementOfSphereGroup(fam,One(F));
    G!.GeneratorsOfMagmaWithInverses := List(GeneratorsOfGroup(F),x->ElementOfSphereGroup(fam,x));
    SetIsSphereGroup(G,true);
    
    return G;
end);

InstallMethod(EpimorphismFromFreeGroup, "(IMG) for a sphere group",
        [IsSphereGroup],
        function(g)
    local f, n;
    n := Length(GeneratorsOfGroup(g));
    f := FreeGroup(n-1);
    return GroupHomomorphismByImages(f,g,GeneratorsOfGroup(f),GeneratorsOfGroup(g){FamilyObj(One(g))!.ordering{[1..n-1]}});
end);

InstallOtherMethod(ExponentSums, "(IMG) for a sphere group element",
        [IsElementOfSphereGroup],
        function(g)
    local l, power, i;
    
    power := FamilyObj(g)!.power;
    l := ExponentSums(UnderlyingElement(g));
    if 0 in power then
        l := l-Minimum(l);
    else
        l := l-l[Length(l)];
    fi;
    for i in [1..Length(l)] do
        if power[i]<>0 then
            l[i] := l[i] mod power[i];
        fi;
    od;
    return l;
end);

################################################################
# conjugacy classes
################################################################
BindGlobal("MAKECYCLICALLYREDUCED@", function(g)
    # returns [m,x] with m minimal representative such that m^x=g
    local w, fam, length, i, m, x, minm, minx;
    
    w := UnderlyingElement(g);
    fam := FamilyObj(g);
    length := Length(w);
    minm := g;
    minx := One(w);
    for i in [1..length] do
        x := Subword(w,i+1,length);
        m := ElementOfSphereGroup(fam,x*Subword(w,1,i));
        if m<minm then minm := m; minx := x; fi;
    od;
    return [minm,ElementOfSphereGroup(fam,minx)];
end);
    
InstallOtherMethod(CyclicallyReducedWord, "(IMG) for a sphere group element",
        [IsElementOfSphereGroup],
        g->MAKECYCLICALLYREDUCED@(g)[1]);

InstallOtherMethod(ConjugacyClass, "(IMG) for a sphere group element",
        [IsElementOfSphereGroup],
        function(g)
    local c;
    c := ConjugacyClass(FamilyObj(g)!.group,CyclicallyReducedWord(g));
    SetIsSphereConjugacyClass(c,true);
    return c;
end);

InstallMethod(IsConjugate, "(IMG) for sphere group elements",
        IsCollsElmsElms,
        [IsSphereGroup,IsElementOfSphereGroup,IsElementOfSphereGroup],
        function(g,u,v)
    return CyclicallyReducedWord(u)=CyclicallyReducedWord(v);
end);

InstallOtherMethod(RepresentativeActionOp, "(IMG) for a sphere group element",
        IsCollsElmsElmsX,
        [IsSphereGroup,IsElementOfSphereGroup,IsElementOfSphereGroup,IsFunction],
        function(g,u,v,act)
    if act<>OnPoints then
        TryNextMethod();
    fi;
    
    u := MAKECYCLICALLYREDUCED@(u);
    v := MAKECYCLICALLYREDUCED@(v);
    
    if u[1]=v[1] then
        return u[2]^-1*v[2];
    else
        return fail;
    fi;
end);

InstallMethod(EQ, "(IMG) for sphere conjugacy classes",
        [IsSphereConjugacyClass,IsSphereConjugacyClass],
        function(c,d)
    return Representative(c)=Representative(d);
end);
InstallMethod(LT, "(IMG) for sphere conjugacy classes",
        [IsSphereConjugacyClass,IsSphereConjugacyClass],
        function(c,d)
    return Representative(c)<Representative(d);
end);

InstallMethod(IntersectionNumber, "(IMG) for sphere conjugacy classes",
        IsIdenticalObj,
        [IsSphereConjugacyClass,IsSphereConjugacyClass],
        function(u,v)
    # returns the geometric intersection number of u and v, following Cohen-Lustig.
    local rank, epi, order, i, j, U, V, m, n, pairs, rel;
    
    # A1. convert u,v to lists over [0,...,2*rank-1], following the ordering
    # in which x1 < x1^-1 < x2 < x2^-1 < ...
    epi := EpimorphismFromFreeGroup(FamilyObj(Representative(u))!.group);
    rank := RankOfFreeGroup(Source(epi));
    u := PreImagesRepresentative(epi,Representative(u));
    v := PreImagesRepresentative(epi,Representative(v));
    
    order := [];
    order{[rank+2..2*rank+1]} := [0,2..2*rank-2];
    order{[rank,rank-1..1]} := [1,3..2*rank-1];
    u := order{LetterRepAssocWord(u)+rank+1};
    v := order{LetterRepAssocWord(v)+rank+1};
    m := Length(u);
    n := Length(v);
    
    # A2. Construct the cyclic permutations of u,v and their inverses
    U := List([1..m],i->u{Concatenation([i..m],[1..i-1])});
    Append(U,List(U,u->Reversed(u+List(u,x->(-1)^x))));
    
    V := List([1..n],i->v{Concatenation([i..n],[1..i-1])});
    Append(V,List(V,v->Reversed(v+List(v,x->(-1)^x))));
    
    # A3. Put their union in cyclic lexicographic ordering
    order := Concatenation(Cartesian([1],[1..m]),Cartesian([1],[-1,-2..-m]),
                     Cartesian([2],[1..n]),Cartesian([2],[-1,-2..-n]));
    SortParallel(Concatenation(U,V),order,
            function(x,y)
        local l, ix, iy;
        if x=y then return false; fi;
        ix := 1; iy := 1; l := -1;
        while x[ix]=y[iy] do
            l := x[ix]; l := l+(-1)^l; # inverse of last letter
            ix := ix+1; if ix>Length(x) then ix := 1; fi;
            iy := iy+1; if iy>Length(y) then iy := 1; fi;
        od;
        return AsSortedList([l<=x[ix],y[iy]<l,x[ix]<y[iy]])=[true,true,false];
    end);

    # A4. find linking pairs
    # A5. construct the equivalence relation
    pairs := [];
    rel := [];
    for i in [1..m] do
        for j in [1..n] do
            if Position(order,[1,i]) < Position(order,[2,j]) and Position(order,[2,j]) < Position(order,[1,-i]) and Position(order,[1,-i]) < Position(order,[2,-j]) then
                AddSet(pairs,[i,j]);
            fi;
            if Position(order,[1,i]) < Position(order,[2,-j]) and Position(order,[2,-j]) < Position(order,[1,-i]) and Position(order,[1,-i]) < Position(order,[2,j]) then
                AddSet(pairs,[i,j]);
            fi;
            if u[i]=v[j] then
                Add(rel,[[i,j],[1+i mod m,1+j mod n]]);
            fi;
            if u[i]+(-1)^u[i]=v[1+j mod n] then
                Add(rel,[[i,1+j mod n],[1+i mod m,j]]);
            fi;
        od;
    od;

    rel := EquivalenceRelationByPairs(Domain(Cartesian([1..m],[1..n])),rel);
    pairs := Set(pairs,p->EquivalenceClassOfElement(rel,p));
    return Size(pairs);
end);

################################################################
# automorphism groups
################################################################

################################################################
# graphs of groups
################################################################

#E sphere.gi . . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
