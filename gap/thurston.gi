#############################################################################
##
#W thurston.gi                                              Laurent Bartholdi
##
#Y Copyright (C) 2011-2013, Laurent Bartholdi
##
#############################################################################
##
##  Thurston's algorithm
##
#############################################################################

BindGlobal("SPIDERRELATOR@", ReturnFail);
BindGlobal("IMGORDERING@", ReturnFail);
BindGlobal("NFFUNCTION@", ReturnFail);

#############################################################################
##
#M IMG Machine to Function
##
BindGlobal("MATCHPERMS@", function(M,q)
    # find a bijection of [1..n] that conjugates M!.output[i] to q[i] for all i
    local c, g, p;
    g := SymmetricGroup(Length(q[1]));
    p := List(GeneratorsOfGroup(StateSet(M)),g->PermList(Output(M,g)));
    q := List(q,PermList);
    c := RepresentativeAction(g,q,p,OnTuples);
    return c;
end);

BindGlobal("MATCHTRANS@", function(M,recur,spider,v)
    # match generators g[i] of M to elements of v.
    # returns a list <w> of elements of <v> such that:
    # if, in M, g[i]^N lifts to a conjugate of g[j] for some integer N, and
    # through <recur> g[i]^N lifts to a conjugate of generator h[k], then
    # set w[j] = v[k].
    # it is in particular assumed that recur[1] has as many lines as
    # StateSet(M) has generators; and that entries in recur[1][j] belong to a
    # free group of rank the length of v.
    local w, i, j, k, c, x, gensM, gensR;

    gensM := GeneratorsOfGroup(StateSet(M));
    gensR := List(GeneratorsOfGroup(spider!.model),x->x^spider!.marking);
    w := [];

    for i in [1..Length(gensM)] do
        x := WreathRecursion(M)(gensM[i]);
        Assert(0,x[2]=recur[2][i]);
        for c in Cycles(PermList(x[2]),AlphabetOfFRObject(M)) do
            j := CyclicallyReducedWord(Product(x[1]{c}));
            k := CyclicallyReducedWord(Product(recur[1][i]{c})^spider!.marking);
            if IsOne(j) then continue; fi;
            j := Position(gensM,j);
            k := PositionProperty(gensR,g->IsConjugate(spider!.group,k,g));
            w[j] := v[k];
        od;
    od;
    Assert(0,BoundPositions(w)=[1..Length(gensM)]);
    return w;
end);

BindGlobal("SPIDERDIST@", function(spiderA,spiderB,fast)
    local model, points, perm, dist, recur, endo, nf, g;

    model := spiderA!.model;

    # try to match feet of spiderA and spiderB
    points := Vertices(spiderA);
    perm := Vertices(spiderB);
    
    perm := MATCHPOINTS@(perm,List(perm,x->points));
    if perm=fail or Set(perm)<>[1..Length(points)] then # no match, find something coarse
        return @.ro*Sum(GeneratorsOfGroup(spiderA!.group),x->Length(PreImagesRepresentative(spiderA!.marking,x)^spiderB!.marking))/Length(points);
    fi;
    
    
    # move points of spiderB to their spiderA matches
    spiderB := COPYSPIDER@(spiderB,points{perm});
    dist := spiderB!.cut!.wiggled;
    
    if fast then # we just wiggled the points, the combinatorics didn't change
        return dist/Length(points);
    fi;
    
    recur := LIFTSPIDER@(spiderA,spiderB,P1z,false);

    if Group(Concatenation(recur[1]))<>spiderA!.group then
        Info(InfoFR,1,"The triangulation got messed up; cross your fingers");
    fi;

    endo := GroupHomomorphismByImagesNC(spiderB!.group,model,
                    GeneratorsOfGroup(spiderB!.group),
        List(recur[1],x->PreImagesRepresentative(spiderA!.marking,x[1])))*spiderB!.marking;

    endo := List(GeneratorsOfGroup(spiderB!.group),x->x^endo);
    REDUCEINNER@(endo,GeneratorsOfMonoid(spiderB!.group),x->x);
    
    for g in endo do
        dist := dist + (Length(g)-1); # if each image is a gen, then endo=1
    od;
    return dist/Length(points);
end);

BindGlobal("PUSHRECURSION@", function(map,M)
    # returns a WreathRecursion() function for Range(map), and not
    # Source(map) = StateSet(M)
    local w;
    w := WreathRecursion(M);
    return function(x)
        local l;
        l := w(PreImagesRepresentative(map,x));
        return [List(l[1],x->Image(map,x)),l[2]];
    end;
end);

BindGlobal("PULLRECURSION@", function(map,M)
    # returns a WreathRecursion() function for Source(map), and not
    # Range(map) = StateSet(M)
    local w;
    w := WreathRecursion(M);
    return function(x)
        local l;
        l := w(Image(map,x));
        return [List(l[1],x->PreImagesRepresentative(map,x)),l[2]];
    end;
end);

BindGlobal("PERRONMATRIX@", function(mat)
    local i, j, len;
    # find if there's an eigenvalue >= 1, without using numerical methods

    len := Length(mat);
    if NullspaceMat(mat-IdentityMat(len))=[] then # no 1 eigenval
        i := List([1..len],i->1);
        j := List([1..len],i->1); # first approximation to perron-frobenius vector
        repeat
            i := i*mat;
            j := j*mat*mat; # j should have all entries growing exponentially
            if ForAll([1..len],a->j[a]=0 or j[a]<i[a]) then
                return false; # perron-frobenius eigenval < 1
            fi;
        until ForAll(j-i,IsPosRat);
    fi;
    return true;
end);

BindGlobal("SURROUNDINGCURVE@", function(t,x)
    # returns a CCW sequence of edges disconnecting x from its complement in t.
    # x is a sequence of indices of vertices. t is a triangulation.
    local starte, a, c, v, e, i;

    starte := First(t!.e,j->j.from.index in x and not j.to.index in x);
    v := starte.from;
    a := [starte.left.pos];
    c := [];
    i := POSITIONID@(v.n,starte);
    repeat
        i := i+1;
        if i > Length(v.n) then i := 1; fi;
        e := v.n[i];
        if e.to.index in x then
            v := e.to;
            e := e.reverse;
            i := POSITIONID@(v.n,e);
        else
            Add(c,e);
            Add(a,e.pos);
            Add(a,e.left.pos);
        fi;
    until IsIdenticalObj(e,starte);
    return [a,c];
end);

BindGlobal("FINDOBSTRUCTION@", function(M,multicurve,spider,boundary)
    # search for an obstruction starting with the elements of M.
    # return fail or a record describing the obstruction.
    # spider and boundary may be "fail".
    local len, w, x, mat, row, i, j, c, d, group, pi, gens, peripheral;

    len := Length(multicurve);
    gens := GeneratorsOfGroup(StateSet(M));
    group := FreeGroup(Length(gens)-1);
    c := IMGRelator(M);
    pi := GroupHomomorphismByImagesNC(StateSet(M),group,List([1..Length(gens)],i->Subword(c,i,i)),Concatenation(GeneratorsOfGroup(group),[Product(List(Reversed(GeneratorsOfGroup(group)),Inverse))]));

    w := PUSHRECURSION@(pi,M);

    peripheral := List(GeneratorsOfSemigroup(StateSet(M)),x->CyclicallyReducedWord(x^pi));
    multicurve := List(multicurve,x->CyclicallyReducedWord(x^pi));
    mat := [];
    for i in multicurve do
        d := w(i);
        row := List([1..len],i->0);
        for i in Cycles(PermList(d[2]),AlphabetOfFRObject(M)) do
            c := CyclicallyReducedWord(Product(d[1]{i}));
            if ForAny(peripheral,x->IsConjugate(group,x,c)) then
                continue; # peripheral curve
            fi;
            j := First([1..len],j->IsConjugate(group,c,multicurve[j])
                       or IsConjugate(group,c^-1,multicurve[j]));
            if j=fail then # add one more curve
                for j in mat do Add(j,0); od;
                Add(row,1/Length(i));
                len := len+1;
                Add(multicurve,c);
            else
                row[j] := row[j] + 1/Length(i);
            fi;
        od;
        Add(mat,row);
    od;

    Info(InfoFR,2,"Thurston matrix is ",mat);

    x := List(EquivalenceClasses(StronglyConnectedComponents(BinaryRelationOnPoints(List([1..len],x->Filtered([1..len],y->IsPosRat(mat[x][y])))))),Elements);
    for i in x do
        if PERRONMATRIX@(mat{i}{i}) then # there's an eigenvalue >= 1
            d := rec(machine := M,
                     obstruction := [],
                     matrix := mat{i}{i});
            if spider<>fail then
                d.spider := spider;
            fi;
            for j in i do
                if spider<>fail and IsBound(boundary[j]) then
                    if not IsBound(spider!.arcs) then spider!.arcs := []; fi;
                    Add(spider!.arcs,[[0,0,255],Float(105/100),boundary[j][1]]);
                fi;
                c := [PreImagesRepresentative(pi,multicurve[j])];
                if spider<>fail then
                    REDUCEINNER@(c,GeneratorsOfMonoid(StateSet(M)),NFFUNCTION@(spider));
                fi;
                Append(d.obstruction,c);
            od;
            return d;
        fi;
    od;
    return fail;
end);

InstallOtherMethod(FindThurstonObstruction, "(IMG) for a list of IMG elements",
#        [IsIMGElementCollection], !method selection doesn't work!
        [IsFRElementCollection],
        function(elts)
    local M;
    M := UnderlyingFRMachine(elts[1]);
    while not IsIMGMachine(M) or ForAny(elts,x->not IsIdenticalObj(M,UnderlyingFRMachine(x))) do
        Error("Elements do not all have the same underlying IMG machine");
    od;
    return FINDOBSTRUCTION@(M,List(elts,InitialState),fail,fail);
end);

BindGlobal("SPIDEROBSTRUCTION@", function(spider,M)
    # check if <spider> has coalesced points; in that case, read the
    # loops around them and check if they form an obstruction
    local multicurve, boundary, i, j, c, d, x, w;

    # construct a list <x> of (lists of vertices that coalesce)
    w := Vertices(spider);
    x := Filtered(Combinations([1..Length(w)],2),p->P1Distance(w[p[1]],w[p[2]])<@.obst);
    x := EquivalenceClasses(EquivalenceRelationByPairs(Domain([1..Length(w)]),x));
    x := Filtered(List(x,Elements),c->Size(c)>1);
    if x=[] then
        return fail;
    fi;

    # replace each x by its conjugacy class
    multicurve := [];
    boundary := [];
    for i in x do
        c := One(spider!.group);
        for j in SpanningTreeBoundary(spider) do
            if (not j.from.index in i) and j.to.index in i then
                c := c*GroupElement(j);
            fi;
        od;
        Add(multicurve,c);
        Add(boundary,SURROUNDINGCURVE@(spider!.cut,i));

    od;
    Info(InfoFR,2,"Testing multicurve ",multicurve," for an obstruction");

    return FINDOBSTRUCTION@(M,List(multicurve,x->PreImagesRepresentative(spider!.marking,x)),spider,boundary);
end);

InstallGlobalFunction(NormalizedQuadraticP1Map, function(f,M,param)
    # param is IsPolynomial, IsBicritical, or a positive integer.
    # in the first case, normalize f as z^d+a_{d-2}*z^{d-2}+...+a_0
    # in the second case, normalize f as (az^d+b)/(cz^d+e)
    # in the third case, normalize f as 1+a/z+b/z^2, such that 0 is on
    # a cycle of length <param>
    # return [new map, Möbius transformation from old to new]
    local p, i, j, k, a, b, mobius, m, coeff, degree, numer, denom, z;
    z := IndeterminateOfUnivariateRationalFunction(f);
    p := POSTCRITICALPOINTS@(f);
    degree := DegreeOfP1Map(f);
    
    if param=fail then
        if IsPolynomialIMGMachine(M) then
            param := IsPolynomial;
        else
            return fail;
        fi;
    fi;
    while not IsPosInt(param) and param<>IsPolynomial and param<>IsBicritical do
        Error("NormalizedQuadraticP1Map: parameter should be 'IsPolynomial', 'IsBicritical' or a positive integer");
    od;
        
    mobius := MoebiusMap(p[2][1][1], p[2][2][1]); # send 2 first c.p. to 0,infty
    
    if param=IsBicritical then # force critical points at 0, infty; make first other point 1
        while Length(p[2])>2 do
            Error("The map is not bicritical, I don't know how to normalize it");
        od;
        j := First(p[3],z->not IsIdenticalObj(z,p[2][1][1]) and not IsIdenticalObj(z,p[2][2][1]));
        if j=fail then # no other point; then map is z^{\pm degree}
            if p[1] then
                return [z^degree,mobius];
            else
                return [z^(-degree),mobius];
            fi;
        fi;
    elif param=IsPolynomial then
        j := fail;
        for i in [1..Length(p[2])] do
            k := First(p[4],r->r[1]=-i)[2];
            if [k,k,degree] in p[4] then j := k; break; fi;
        od;
        while not p[1] or j=fail do
            Error("Map is not a polynomial");
        od;
        if Length(p[3])=2 then
            return [z^degree,mobius];
        fi;
        k := First([1..Length(p[2])],k->p[2][k][1]<>p[3][j]); # other critical point
        mobius := MoebiusMap(p[2][k][1],p[3][j]);
        i := CleanedP1Map(ConjugatedP1Map(f,mobius),@.p1eps); # polynomial
        coeff := CoefficientsOfP1Map(i);        
        mobius := CompositionP1Map(mobius,P1MapSL2([[@.o/coeff[1][degree+1]^(1/(degree-1)),-coeff[1][degree]/degree/coeff[1][degree+1]],[@.z,@.o]]));
        m := P1MapSL2([[Exp(@.2ipi/(degree-1)),@.z],[@.z,@.o]]);
    if false then #!!! disabled; this could be an infinite loop
        j := SupportingRays(M);
        
        repeat
            k := SupportingRays(IMGMachine(ConjugatedP1Map(f,mobius)));
            if k=j then break; fi; #!!! maybe should be: "equivalent"rays?
            mobius := CompositionP1Map(mobius,m);
            Info(InfoFR,1,"param:=IsPolynomial: trying to rotate around infinity");
        until false;
    fi;
    else # parameterize on slice V_n
        while degree<>2 do
            Error("'param:=<positive integer>' only makes sense for degree-2 maps");
        od;
        for i in [1..2] do
            j := [First(p[4],r->r[1]=-i)[2]];
            for k in [1..param] do
                j[k+1] := First(p[4],r->r[1]=j[k])[2];
            od;
            if j[param+1]=j[1] and not j[1] in j{[2..param]} then j := j[1]; break; fi;
        od;
        while not IsInt(j) do
            Error("I couldn't find a cycle of length ",param);
        od;
        if param=1 then # map marked point to infty, unmarked to 0, 0=>1
            if Length(p[3])=2 then
                return [z^degree,mobius];
            fi;
            j := First(p[4],r->r[1]=i-3)[2];
            mobius := MoebiusMap(p[2][3-i][1],p[3][j],p[2][i][1]);
        elif param=2 then # normalize as a/(z^2+2z), infty=>0->infty, -1=>
            if Length(p[3])=2 then # special case infty=>0=>infty or polynomial
                if First(p[4],r->r[1]=j)[2]=j then
                    return [z^degree,mobius];
                else
                    return [z^(-degree),mobius];
                fi;
            fi;
            mobius := CompositionP1Map(MoebiusMap(p[3][j],p[2][3-i][1],p[2][i][1]),-z);
        else # normalize as 1+a/z+b/z^2, 0=>infty->1
            k := First(p[4],r->r[1]=j)[2];
            mobius := MoebiusMap(p[2][i][1],p[3][k],p[3][j]);
        fi;
    fi;
    f := CleanedP1Map(ConjugatedP1Map(f,mobius),@.p1eps);
    numer := List([0..degree],i->@.z);
    denom := List([0..degree],i->@.z);
    if param=2 then # special cleanup
        coeff := CoefficientsOfP1Map(f);
        if IsZero(coeff[1][2]) and IsZero(coeff[1][3]) and IsZero(coeff[2][1]) and AbsoluteValue(coeff[2][3]-1/2)<@.ratprec then
            numer[1] := coeff[1][1]*2; # numerator is A
            denom[2] := @.o; # denominator is z+z^2/2, off by 2 now
            denom[3] := @.o/2;
            f := CleanedP1Map(P1MapByCoefficients(numer,denom),@.p1eps);
            coeff := CoefficientsOfP1Map(f);
            f := P1MapByCoefficients(coeff[1],2*coeff[2]);
        else
            Error("Cannot normalize to A/(z^2+2z)");
        fi;
    elif param=3 then # special cleanup, force 1+a+b=0
        coeff := CoefficientsOfP1Map(f);
        if AbsoluteValue(Sum(coeff[1]))<@.ratprec and AbsoluteValue(coeff[1][3]-1)<@.ratprec then
            numer[1] := coeff[1][1];
            numer[2] := -numer[1]-@.o;
            numer[3] := @.o;
            denom[3] := @.o;
            f := P1MapByCoefficients(numer,denom);
        else
            Error("Cannot normalize to (A+(-1-A)z+z^2)/z^2");
        fi;
    fi;
    return [f,mobius];
end);

BindGlobal("SIMPLIFYBYBRAIDTWISTS@", function(m,marking)
    local bg, g, i, len, newmarking, newlen, twist, idle;

    bg := BraidTwists(m);
    g := GeneratorsOfGroup(Source(marking));
    len := infinity;
    twist := [];
    repeat
        idle := true;
        newmarking := List(bg,t->t*marking);
        newlen := List(newmarking,t->Sum(g,x->Length(x^t)));
        i := Position(newlen,Minimum(newlen));
        if newlen[i] < len then
            len := newlen[i];
            marking := newmarking[i];
            Add(twist,bg[i]);
            idle := false;
        fi;
    until idle;
    
    return List(Reversed(twist),Inverse);
end);

BindGlobal("FRMACHINE2RAT@", function(z,M)
    local oldspider, spider, t, gens, n, deg, model,
          f, mobius, match, v, i, j, recf, recmobius, map,
          dist, obstruction, lifts, sublifts, fast, poly;

    MAKEP1EPS@();
    
    model := StateSet(M);
    gens := GeneratorsOfGroup(model);
    n := Length(gens);
    deg := Length(AlphabetOfFRObject(M));
    poly := IsPolynomialFRMachine(M);
    
    if n=2 then # special handling, space is not hyperbolic
        i := Sum(List(Transitions(M,1),ExponentSums));
        if i[1]-i[2]=1 then
            return z^deg;
        elif i[1]-i[2]=-1 then
            return z^(-deg);
        else
            Error(M," is not an IMG machine");
        fi;
    fi;    
    
    # create spider on +/- equidistributed points on Greenwich meridian.
    # its spanning tree will be consecutive edges from infty to 1,
    # and so its IMG ordering is predictably that of M
    v := [];
    for i in [1..n] do
        # on positive real axis, tending to infinity
        f := Exp(@.2ipi*(i-1)/(2*n-2));
        Add(v,P1Point(@.o*ImaginaryPart(f)/(@.ro+RealPart(f))));
    od;
    v := Permuted(v,PermList(IMGORDERING@(M)));
    spider := TRIVIALSPIDER@(v);
    IMGMARKING@(spider,model);

    if ValueOption("julia")<>fail then
        i := ValueOption("julia");
        if not IsInt(i) then i := 1000; fi; # number of points to trace
        spider!.points := EquidistributedP1Points(i);
    fi;
    
    lifts := fail;
    f := fail; # in the beginning, we don't know them
    fast := false;

    repeat
        oldspider := spider;
        # find a rational map that has the right critical values
        f := RATIONALMAP@(z,spider,List(gens,g->Output(M,g)),f,lifts);
        lifts := f[2]; f := f[1];
        Info(InfoFR,3,"1: found rational map ",f," on vertices ",lifts);

        if fast then # just get points closest to those in spider t
            match := MATCHPOINTS@(sublifts,List(sublifts,x->lifts));
            if match=fail then
		Info(InfoFR,3,"1.5: back to slow mode");
                fast := false; continue;
            fi;
            sublifts := lifts{match};
        else
            # create a spider on the full preimage of the points of <spider>
            t := TRIVIALSPIDER@(lifts);
            IMGMARKING@(t,FreeGroup(Length(lifts)));
            Info(InfoFR,3,"2: created liftedspider ",t);

            # lift paths in <spider> to <t>
            recf := LIFTSPIDER@(t,spider,f,poly);
            if recf=fail then return fail; fi;
            recf := IMGRECURSION@(t,spider,recf[1],recf[2],false);
            Assert(1, CHECKREC@(recf,spider!.ordering,NFFUNCTION@(t)));
            Info(InfoFR,3,"3: recursion ",recf);
            
            # find a bijection between the alphabets of <recf> and <M>
            match := MATCHPERMS@(M,recf[2]);
            if match=fail then return fail; fi;
            
            REORDERREC@(recf,match);
            Info(InfoFR,3,"4: alphabet permutation ",match);

            # extract those vertices in <v> that appear in the recursion
            sublifts := MATCHTRANS@(M,recf,t,lifts);
            Info(InfoFR,3,"5: extracted and sorted vertices ",sublifts);
        fi;
 
        # find a mobius transformation that normalizes <sublifts> wrt PSL2C
        mobius := NORMALIZINGMAP@(sublifts,Vertices(spider));
        Info(InfoFR,3,"6: normalize by mobius map ",mobius);

        # now create the new spider on the image of these points
        v := List(sublifts,p->P1Image(mobius,p));
        
	if fast then
            dist := Sum([1..Length(v)],i->P1Distance(Vertices(spider)[i],v[i]));
            if dist>@.fast*Length(v) then
                fast := false;
                Info(InfoFR,3,"7: legs moved ",dist,"; back to slow mode");
                spider := oldspider;
                continue; # restart
            fi;
            # just wiggle spider around
            spider := COPYSPIDER@(spider,v);
        else
            spider := TRIVIALSPIDER@(v);
            recmobius := LIFTSPIDER@(spider,t,InverseP1Map(mobius),poly);
            if recmobius=fail then
                return fail;
            else
                recmobius := recmobius[1];
            fi;
            Info(InfoFR,3,"7: new spider ",spider," with recursion ",recmobius);

            # compose recursion of f with that of mobius
            map := t!.marking*GroupHomomorphismByImagesNC(t!.group,spider!.group,GeneratorsOfGroup(t!.group),List(recmobius,x->x[1]));
            for i in recf[1] do
                for j in [1..Length(i)] do i[j] := i[j]^map; od;
            od;
            Info(InfoFR,3,"8: composed recursion is ",recf);
                   
            # finally set marking of new spider using M
            spider!.model := model;
            spider!.ordering := oldspider!.ordering;
            spider!.marking := MATCHMARKINGS@(M,spider!.group,recf);
            Assert(1, CHECKREC@(recf,spider!.ordering,NFFUNCTION@(t)));
            Assert(1, CHECKSPIDER@(spider));
            Info(InfoFR,3,"9: marked new spider ",spider);
        fi;
        
        dist := SPIDERDIST@(spider,oldspider,fast);
        Info(InfoFR,2,"Spider moved ",dist," steps; feet=",Vertices(spider)," marking=",spider!.marking);

        if dist<@.ratprec then
            if fast then # force one last run with the full algorithm
                fast := false;
                continue;
            fi;
            f := CleanedP1Map(CompositionP1Map(f,InverseP1Map(mobius)),@.p1eps);
            break;
        elif dist<@.fast then
            fast := true;
        else
            fast := false;
        fi;
        obstruction := SPIDEROBSTRUCTION@(spider,M);
        if obstruction<>fail then
            return obstruction;
        fi;
    until false;
    
    Info(InfoFR,2,"Spider converged");
    
    i := NormalizedQuadraticP1Map(f,M,ValueOption("param"));
    if i<>fail then
        spider := COPYSPIDER@(spider,InverseP1Map(i[2]));
        f := i[1];
    fi;
    
    # construct a new machine with simpler recursion
    for i in recf[1] do
        for j in [1..Length(i)] do
            i[j] := PreImagesRepresentative(spider!.marking,i[j]);
        od;
    od;
    IMGOPTIMIZE@(recf[1], recf[2], SPIDERRELATOR@(spider),poly);
    t := FRMachine(model, recf[1], recf[2]);
    SetIMGRelator(t, SPIDERRELATOR@(spider));
    
    # in principle, now, t is the same as M (but perhaps in a slightly
    # different basis, i.e. conjugated by an inner automorphism)
    
    # find a good twist to shorten the marking; presumably, this
    # will give a simpler machine
    i := SIMPLIFYBYBRAIDTWISTS@(M, spider!.marking);
    t := t^Product(i);
    SetCorrespondence(t,i);
    
    spider!.map := f;
    spider!.cycle := ATTRACTINGCYCLES@(POSTCRITICALPOINTS@(f));
                     
    return [f,t,spider];
end);

InstallMethod(P1Map, "(IMG) for a sphere machine",
        [IsSphereMachine],
        function(M)
    local data, f;
    data := FRMACHINE2RAT@(P1z,M);
    if not IsList(data) then return data; fi;
    f := data[1];
    SetIMGMachine(f,data[2]);
    SetSpider(f,data[3]);
    return f;
end);

#E thurston.gi . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
