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

InstallMethod(ElementOfSphereGroup, "(IMG) for a sphere element family and a word",
        [IsFamily,IsAssocWordWithInverse],
        function(fam,w)
    return ElementOfFpGroup(fam,FpElementNFFunction(fam)(w));
end);

InstallMethod(LetterRepAssocWord, "(IMG) for a sphere group element",
        [IsElementOfSphereGroup],
        function(g)
    return LetterRepAssocWord(UnderlyingElement(g));
end);

InstallMethod(AssocWordByLetterRep, "(IMG) for a sphere group family, and word",
        [IsElementOfSphereGroupFamily,IsList],
        function(F,l)
    return ElementOfSphereGroup(F,AssocWordByLetterRep(FamilyObj(One(F!.freeGroup)),l));
end);

#!!!! this needs improvement, it's too inefficient
BindGlobal("BRUTEFORCEFACTORIZATION@", function(w,g)
    local epi, l, len, newlen, v, neww;

    if not w in g then return fail; fi;

    epi := EpimorphismFromFreeGroup(g);
    l := [];
    len := Length(w);
    while len>0 do
        for v in Source(epi) do
            neww := v^epi*w;
            newlen := Length(neww);
            if newlen<len then
                w := neww;
                len := newlen;
                Append(l,LetterRepAssocWord(v^-1));
                break;
            fi;
        od;
    od;
    return l;
end);

InstallOtherMethod(AsWordLetterRepInGenerators, "(IMG) for a sphere group element and a subgroup",
        [IsElementOfSphereGroup, IsGroup and HasGeneratorsOfGroup],
        function(w,g)
    local l, G, gens, iso, i;

    G := FamilyObj(w)!.group;
    gens := GeneratorsOfGroup(g);

    # first, try to use the FGA algorithms
    if HasIsomorphismFreeGroup(G) then
        iso := IsomorphismFreeGroup(G);
        return AsWordLetterRepInGenerators(w^iso,Group(List(gens,x->x^iso)));
    fi;

    if not IsBound(g!.freelift) then
        # we should(!) add the normal closure of all relations; but this seems to work well in
        # practice.
        g!.freelift := Group(Concatenation(List(gens,UnderlyingElement),RelatorsOfFpGroup(G)));
    fi;

    # a cheap shot
    l := AsWordLetterRepInGenerators(UnderlyingElement(w),g!.freelift);
    if l<>fail then
        w := [];
        gens := [-Length(gens)..Length(gens)];
        for i in l do if i in gens then Add(w,i); fi; od;
        return w;
    fi;

    Info(InfoIMG,1,"AsWordLetterRepInGenerators: resorting to brute-force factorization; this could be slow");
    # now complicated and inefficient method
    return BRUTEFORCEFACTORIZATION@(w,g);
end);

InstallMethod(Subword, "(IMG) for a sphere group element",
        [IsElementOfSphereGroup,IsPosInt,IsPosInt],
        function(g,i,j)
    return ElementOfSphereGroup(FamilyObj(g),Subword(UnderlyingElement(g),i,j));
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

# !!!! still doesn't work for [0,2,2,2], [2,2,0,2], [0,3,2,3,3,2], [4x2, 1x0], ...
# maybe imbed in coxeter system with generators s_i and relations (s_is_{i+1})^power[i]?

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
    elif IsSubset(power,[1,2]) then # force non-consecutive indices
        if IsEvenInt(n) then
            ordering := ShortLexOrdering(freemon,gens{Concatenation(ordering{Concatenation([1,3..n-1],[2,4..n])},[n+1..2*n])});
        else
            weight := Concatenation(ListWithIdenticalEntries(n,1),ListWithIdenticalEntries(n,n));
            weight[ordering[n]] := n;
            special := ShortLexOrdering(freemon);
            ordering := WeightLexOrdering(freemon,gens{Concatenation(ordering{Concatenation([1,3..n],[2,4..n-1])},[n+1..2*n])},weight);
        fi;
    elif 0 in power and Set(power)<>[0] then
        weight := ListWithIdenticalEntries(2*n,1);
        weight[Position(power,0)] := n;
        weight[Position(power,0)+n] := n;
        ordering := WeightLexOrdering(freemon,gens,weight);
        special := ShortLexOrdering(freemon);
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

InstallMethod(IsomorphismFreeGroup, "(IMG) for a sphere group",
        [IsSphereGroup],
        function(G)
    local F, n, gens, img, pregens, preimg;
    
    gens := GeneratorsOfGroup(G);
    n := Length(gens);
    F := FreeGroup(n-1);
    pregens := GeneratorsOfGroup(F);
    img := [];
    img{OrderingOfSphereGroup(G)} := Concatenation(pregens,[Product(pregens)^-1]);
    preimg := gens{OrderingOfSphereGroup(G){[1..n-1]}};

    return GroupHomomorphismByFunction(G,F,w->MappedWord(w,gens,img),w->MappedWord(w,pregens,preimg));
end);

InstallGlobalFunction(SphereGroup, function (arg)
    local  F, G, fam, rel, n, ordering, power;
    
    while not (Length(arg) in [1..3] and (IsPosInt(arg[1]) or IsList(arg[1])) and (Length(arg)=1 or IsList(arg[2]))) do
        Error("Usage: SphereGroup(num_points or ordering [, degrees] [,group]");
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

    # forbid degenerate cases
    if n<2 then
        Error("SphereGroup needs at least 2 generators");
    fi;
    if n=2 and power[1]<>power[2] then
        Error("When SphereGroup is called with 2 generators, the orders must be the same");
    fi;

    # create the underlying f.p. group
    if Length(arg)=3 then
        while not IsFreeGroup(arg[3]) do
            Error("Spheregroup requires a free group as optional third argument");
        od;
        F := arg[3];
    else
        F := FreeGroup(n);
    fi;

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
    
    SetFpElementNFFunction(fam,SPHERENFFUNCTION@(G,ordering,power));

    # reset some of the attributes of G that were erroneously computed
    G!.OneImmutable := ElementOfSphereGroup(fam,One(F));
    fam!.OneImmutable := G!.OneImmutable;
    G!.GeneratorsOfMagmaWithInverses := List(GeneratorsOfGroup(F),x->ElementOfSphereGroup(fam,x));

    SetIsSphereGroup(G,true);

    # set isomorphism to free group, if available
    if ForAll(power,x->x=0) then
        IsomorphismFreeGroup(G);
    fi;

    return G;
end);

InstallMethod(OrderingOfSphereGroup, "(IMG) for a sphere group, get family attribute",
        [IsSphereGroup],
        g->FamilyObj(One(g))!.ordering);

InstallMethod(ExponentsOfSphereGroup, "(IMG) for a sphere group, get family attribute",
        [IsSphereGroup],
        g->FamilyObj(One(g))!.power);

InstallMethod(IsomorphismSphereGroup, "(IMG) for a f.p. group",
        [IsFpGroup],
        function(G)
    local newG;
    newG := AsSphereGroup(G);
    if newG=fail then
        return fail;
    else
        return GroupHomomorphismByImages(G,newG);
    fi;
end);
    
InstallMethod(AsSphereGroup, "(IMG) for a f.p. group",
        [IsFpGroup],
        function(G)
    local power, ordering, gens, r, i, letters;

    gens := GeneratorsOfGroup(G);
    power := List(gens,x->0);
    ordering := fail;

    for r in RelatorsOfFpGroup(G) do
        if IsOne(r) then
            continue;
        fi;
        letters := LetterRepAssocWord(r);
        if letters[1]<0 then
            letters := LetterRepAssocWord(r^-1);
        fi;
        if ForAll(letters,l->l=letters[1]) then
            power[letters[1]] := Gcd(power[letters[1]],Length(letters));
            continue;
        fi;
        if AsSortedList(letters)=[1..Length(gens)] then
            ordering := letters;
            continue;
        fi;
        Error("AsSphereGroup: illegal relator ",r);
    od;
    while ordering=fail do
        Error("AsSphereGroup: no ordering relation");
    od;
    return SphereGroup(ordering,power,FreeGroupOfFpGroup(G));
end);

InstallMethod(EpimorphismFromFreeGroup, "(IMG) for a sphere group",
        [IsSphereGroup],
        function(g)
    local f, n;
    n := Length(GeneratorsOfGroup(g));
    f := FreeGroup(n-1);
    return GroupHomomorphismByImages(f,g,GeneratorsOfGroup(f),GeneratorsOfGroup(g){OrderingOfSphereGroup(g){[1..n-1]}});
end);

InstallMethod(EulerCharacteristic, "(IMG) for a free group",
        [IsFreeGroup],
        g->1-RankOfFreeGroup(g));

InstallMethod(EulerCharacteristic, "(IMG) for a sphere group",
        [IsSphereGroup],
        function(g)
    local chi, o;
    chi := 2;
    for o in ExponentsOfSphereGroup(g) do
        if o=0 then
            chi := chi - 1;
        else
            chi := chi - (1 - 1/o);
        fi;
    od;
    if chi>0 then chi := chi/2; fi;
    return chi;
end);

InstallMethod(RankOfSphereGroup, "(IMG) for a sphere group",
        [IsSphereGroup],
        g->RankOfFreeGroup(FreeGroupOfFpGroup(g)));

InstallOtherMethod(ExponentSums, "(IMG) for a sphere group element",
        [IsElementOfSphereGroup],
        function(g)
    local l, power, i;
    
    power := ExponentsOfSphereGroup(g);
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
InstallMethod(PeripheralClasses, "(IMG) for a sphere group",
        [IsSphereGroup],
        function(g)
    return List(GeneratorsOfGroup(g),x->ConjugacyClass(g,x));
end);

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

BindGlobal("ISPERIPHERAL@",
        function(g)
    g := MAKECYCLICALLYREDUCED@(g)[1];
    return IsOne(g) or g=Subword(g,1,1)^Length(g);
end);
    
InstallMethod(IsPeripheral, "(IMG) for a sphere group element",
        [IsElementOfSphereGroup],
        ISPERIPHERAL@);


InstallMethod(IsPeripheral, "(IMG) for a sphere conjugacy class",
        [IsSphereConjugacyClass],
        c->ISPERIPHERAL@(Representative(c)));

InstallMethod(Order, "(IMG) for a sphere group element",
        [IsElementOfSphereGroup],
        function(g)
    local order, w, len;

    if EulerCharacteristic(FamilyObj(g)!.group)>0 then # finite group
        TryNextMethod();
    fi;

    g := MAKECYCLICALLYREDUCED@(g)[1];
    w := LetterRepAssocWord(g);
    len := Length(w);
    if len=0 then return 1; fi;
    if ForAll(w,x->x=w[1]) then
        order := FamilyObj(g)!.power[AbsoluteValue(w[1])];
        if order=0 then
            return infinity;
        else
            return order/Gcd(order,len);
        fi;
    else
        return infinity;
    fi;
end);

InstallOtherMethod(ConjugacyClass, "(IMG) for a sphere group element",
        [IsElementOfSphereGroup],
        function(g)
    local c;
    c := ConjugacyClass(FamilyObj(g)!.group,CyclicallyReducedWord(g));
    SetIsSphereConjugacyClass(c,true);
    return c;
end);

InstallMethod(IsSphereConjugacyClass, "(IMG) for a group conj. class",
        [IsConjugacyClassGroupRep  and IsAssociativeElementCollection and IsMultiplicativeElementWithInverseCollection],
        function(c)
    return IsElementOfSphereGroup(Representative(c));
end);

Perform([InverseImmutable,InverseSameMutability,InverseMutable],
        function(method)
    InstallMethod(method, "(IMG) for a sphere conj. class",
            [IsSphereConjugacyClass],
            function(c)
        return ConjugacyClass(Inverse(Representative(c)));
    end);
end);

InstallMethod(POW, "(IMG) for a sphere conj. class and an integer",
        [IsSphereConjugacyClass,IsInt],
        function(c,n)
    return ConjugacyClass(Representative(c)^n);
end);

InstallMethod(IsSphereConjugacyClassCollection, "(IMG) for a group conj. class",
        [IsHomogeneousList and IsAssociativeElementCollColl and IsMultiplicativeElementWithInverseCollColl],
        function(l)
    local i, g;
    if l=[] or not IsSphereConjugacyClass(l[1]) then return false; fi;
    g := FamilyObj(Representative(l[1]))!.group;
    for i in [2..Length(l)] do
        if not IsSphereConjugacyClass(l[i]) or FamilyObj(Representative(l[i]))!.group<>g then
            return false;
        fi;
    od;
    return true;
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

BindGlobal("ISCYCLICALLYORDERED@",
        function(arg)
    return ForAny(CyclicGroup(IsPermGroup,Length(arg)),g->IsSet(Permuted(arg,g)));
end);

InstallMethod(IntersectionNumber, "(IMG) for sphere conjugacy classes",
        IsIdenticalObj,
        [IsSphereConjugacyClass,IsSphereConjugacyClass],
        function(u,v)
    # returns the geometric intersection number of u and v, following Cohen-Lustig.
    local rank, epi, order, i, j, U, V, UV, m, n, pairs, rel;
    
    # A1. convert u,v to lists over [0,...,2*rank-1], following the ordering
    # in which x1 < x1^-1 < x2 < x2^-1 < ...
    epi := EpimorphismFromFreeGroup(FamilyObj(Representative(u))!.group);
    rank := RankOfFreeGroup(Source(epi));
    u := PreImagesRepresentative(epi,Representative(u));
    v := PreImagesRepresentative(epi,Representative(v));
    
    order := [];
    order{[rank+2..2*rank+1]} := [0,2..2*rank-2];
    order{[rank,rank-1..1]} := [1,3..2*rank-1];
    u := order{LetterRepAssocWord(CyclicallyReducedWord(u))+rank+1};
    v := order{LetterRepAssocWord(CyclicallyReducedWord(v))+rank+1};
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
    rel := function(x,y)
        local l, i;
        if x=y then return false; fi;
        i := 1; l := -1;
        while x[i]=y[i] do
            l := x[i]; l := l+(-1)^l; # inverse of last letter
            i := i+1;
        od;
        return AsSortedList([l<=x[i],y[i]<l,x[i]<y[i]])=[true,true,false];
    end;
    UV := List(Concatenation(U,V),w->PeriodicList([],w));
    SortParallel(UV,order,rel);
    order := List([1..2],x->List([1,-1],s->List([1..ELM_LIST([m,n],x)],i->Position(UV,UV[Position(order,[x,s*i])]))));
    
    # A4. find linking pairs
    # A5. construct the equivalence relation
    pairs := [];
    rel := [];
    for i in [1..m] do
        for j in [1..n] do
            if ISCYCLICALLYORDERED@(order[1][1][i],order[2][1][j],order[1][2][i],order[2][2][j]) or
               ISCYCLICALLYORDERED@(order[1][1][i],order[2][2][j],order[1][2][i],order[2][1][j]) then
                AddSet(pairs,[i,j]);
            fi;
            if u[i]=v[j] then
                Add(rel,[[i,j],[1+i mod m,1+j mod n]]);
            fi;
            if u[i]+(-1)^u[i]=v[j] then
                Add(rel,[[i,1+j mod n],[1+i mod m,j]]);
            fi;
        od;
    od;
    
    # A6. Count the equivalence classes
    rel := EquivalenceRelationByPairs(Domain(Cartesian([1..m],[1..n])),rel);
    pairs := Set(pairs,p->EquivalenceClassOfElement(rel,p));
    return Size(pairs);
end);

InstallMethod(SelfIntersectionNumber, "(IMG) for sphere conjugacy classes",
        [IsSphereConjugacyClass],
        function(u)
    return IntersectionNumber(u,u)/2;
end);

################################################################
# automorphism groups
################################################################
InstallMethod(AutomorphismGroup, "(IMG) for a sphere group",
        [IsSphereGroup],
        function(g)
    local r, i, j, m, ni, nj, gens, img, aut, relator, ordering, a, inner;

    ordering := OrderingOfSphereGroup(g);
    gens := GeneratorsOfGroup(g);
        
    aut := [];
    inner := [];
    for ni in [1..Length(ordering)] do
        for nj in [ni+1..Length(ordering)] do
            m := Product(gens{ordering{[ni+1..nj-1]}},One(g));
            i := ordering[ni];
            j := ordering[nj];
            img := ShallowCopy(gens);
            img[i] := gens[i]^(m*gens[j]/m);
            img[j] := gens[j]^(m^-1*gens[i]*m*gens[j]);
            Add(aut,GroupHomomorphismByImages(g,g,gens,img));
        od;
        Add(inner,InnerAutomorphism(g,gens[ordering[ni]]));
    od;
    a := Group(Concatenation(aut,inner));
    inner := Group(inner);
    SetIsAutomorphismGroupOfSphereGroup(a,true);
    SetInnerAutomorphismsAutomorphismGroup(a,inner);
    SetParent(inner,a);
    SetIsNormalInParent(inner,true);
    return a;
end);

BindGlobal("FACTORIZEAUT@", function(ggens,a,fp,gens,extgens,outergens,invertible)
    local epi, out;

    epi := EpimorphismFromFreeGroup(a);
    out := GroupHomomorphismByFunction(a,fp,function(x)
        local w, n, newx, newn, g;
        w := One(Source(epi));
        n := infinity;
        while n>0 do
            for g in Source(epi) do
                newx := x*g^epi;
                newn := Sum(ggens,s->Length(s^newx)-1);
                if newn < n then x := newx; w := g^-1*w; n := newn; break; fi;
            od;
        od;
        return MappedWord(w,GeneratorsOfGroup(Source(epi)),extgens);
    end,invertible,w->MappedWord(w,gens,outergens));

    return out;
end);

BindGlobal("SIMPLIFYBYBRAIDTWISTS@", function(g,hom)
    # hom is a homomorphism from a group h to g.
    # find a product of automorphisms of g such that hom, post-composed by these, becomes
    # simpler, for the "sum of lengths of images of generators" metric.
    # return the list of these automorphisms.
    local gens, autgens, i, len, newhom, newlen, twist, idle;

    autgens := GeneratorsOfGroup(AutomorphismGroup(g));
    gens := GeneratorsOfGroup(Source(hom));
    len := infinity;
    twist := [];
    repeat
        idle := true;
        newhom := hom * autgens;
        newlen := List(newhom,t->Sum(gens,x->Length(x^t)));
        i := Position(newlen,Minimum(newlen));
        if newlen[i] < len then
            len := newlen[i];
            hom := newhom[i];
            Add(twist,autgens[i]);
            idle := false;
        fi;
    until idle;
    
    return twist;
end);

InstallMethod(IsomorphismFpGroup, "(IMG) for a sphere automorphism group",
        [IsAutomorphismGroupOfSphereGroup],
        function(a)
    local g, ordering, i, j, fp, gens;

    g := Source(One(a));
    ordering := OrderingOfSphereGroup(g);
    gens := [];
    for i in [1..Length(ordering)] do
        for j in [i+1..Length(ordering)] do
            Add(gens,Concatenation("t",String(ordering[i]),String(ordering[j])));
        od;
    od;
    for i in [1..Length(ordering)] do
        Add(gens,Concatenation("i",String(ordering[i])));
    od;
    fp := FreeGroup(gens);
    gens := GeneratorsOfGroup(fp);
# add relations!!!

# special cases: if #ggens<=3 then fp=g.

    return FACTORIZEAUT@(GeneratorsOfGroup(g),a,fp,gens,gens,GeneratorsOfGroup(a),true);
end);

InstallMethod(EpimorphismToOut, "(IMG) for a sphere automorphism group",
        [IsAutomorphismGroupOfSphereGroup],
        function(a)
    local g, ordering, i, j, fp, gens;

    g := Source(One(a));
    ordering := OrderingOfSphereGroup(g);
    gens := [];
    for i in [1..Length(ordering)] do
        for j in [i+1..Length(ordering)] do
            Add(gens,Concatenation("t",String(ordering[i]),String(ordering[j])));
        od;
    od;
    fp := FreeGroup(gens);
    gens := GeneratorsOfGroup(fp);
# add relations!!!

# special cases: if #ggens<=3 then fp is trivial. if #ggens=4 then fp is a sphere group

    return FACTORIZEAUT@(GeneratorsOfGroup(g),a,fp,gens,Concatenation(gens,List(ordering,x->One(fp))),GeneratorsOfGroup(a){[1..Length(gens)]},false);
end);

InstallMethod(NaturalHomomorphismByNormalSubgroupOp, "(IMG) for a sphere automorphism group",
        [IsAutomorphismGroupOfSphereGroup,IsGroup],
        function(a,i)
    if not IsIdenticalObj(i,InnerAutomorphismsAutomorphismGroup(a)) then
        TryNextMethod();
    fi;

    return EpimorphismToOut(a);
end);

################################################################
# graphs of groups
################################################################
InstallMethod(AmalgamatedFreeProduct, "(IMG) for two sphere groups and two elements",
        [IsSphereGroup,IsSphereGroup,IsElementOfSphereGroup,IsElementOfSphereGroup],
        function(G1,G2,x1,x2)
    # free product of G1 and G2 modulo relation x1x2=1
    local A, G, x, gens, order, exp, embed, img, i, c, names;

    G := [G1,G2];
    x := [x1,x2];
    gens := List(G,GeneratorsOfGroup); gens := List(gens,ShallowCopy);
    order := List(G,OrderingOfSphereGroup);
    exp := List(G,ExponentsOfSphereGroup); exp := List(exp,ShallowCopy);
    while not ForAll([1..2],i->x[i] in gens[i] and Order(x[i])=infinity) do
        Error("Amalgamated free products work (for now) only along infinite-order generators");
    od;

    x := List([1..2],i->Position(gens[i],x[i]));
    img := [];
    for i in [1..2] do
        Remove(gens[i],x[i]); Remove(exp[i],x[i]);
        img[i] := [];
        img[i]{Difference([1..Length(gens[i])+1],[x[i]])} := (i-1)*Length(gens[1])+[1..Length(gens[i])];
    od;
    x := List([1..2],i->Position(order[i],x[i]));
    
    order := List([1..2],i->Concatenation(order[i]{[x[i]+1..Length(order[i])]},order[i]{[1..x[i]-1]}));

    A := SphereGroup(Concatenation(List([1..2],i->img[i]{order[i]})),Concatenation(exp));
    embed := List([1..2],i->GroupHomomorphismByImages(G[i],A,gens[i],GeneratorsOfGroup(A){(i-1)*Length(gens[1])+[1..Length(gens[i])]}));
    SetEmbeddingsOfAmalgamatedFreeProduct(A,embed);

# a hack, now: rename the generators more nicely

    names := List(gens,g->List(g,String));
    if Intersection(names[1],names[2])<>[] then
        c := names[1][1][1]; # first character
        if c<>'z' and ForAll(names,list->ForAll(list,w->w[1]=c)) then # increment letter on half
            names[2] := List(names[2],ShallowCopy);
            for i in names[2] do i[1] := CHAR_INT(INT_CHAR(c)+1); od;
        else
            names := List([1..2],i->List(names[i],s->Concatenation(s,String(i))));
        fi;
    fi;
    FamilyObj(One(A)![1])!.names := Concatenation(names);

    return A;
end);

InstallMethod(Embedding, "(IMG) for a sphere group and an index",
        [IsSphereGroup and HasEmbeddingsOfAmalgamatedFreeProduct,IsInt],
        function(G,n)
    return EmbeddingsOfAmalgamatedFreeProduct(G)[n];
end);

#E sphere.gi . . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
