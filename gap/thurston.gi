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

#############################################################################
##
#M Sphere Machine to Function
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

BindGlobal("MATCHTRANS@", function(M1,M2)
    # match generators g[i] of M to elements of v.
    # returns a list <w> of elements of <v> such that:
    # if, in M, g[i]^N lifts to a conjugate of g[j] for some integer N, and
    # through <recur> g[i]^N lifts to a conjugate of generator h[k], then
    # set w[j] = v[k].
    # it is in particular assumed that recur[1] has as many lines as
    # StateSet(M) has generators; and that entries in recur[1][j] belong to a
    # free group of rank the length of v.
    local w, i, j, k, c, g1, g2, gens1, gens2, src, dst;

    g1 := StateSet(M1);
    g2 := FamilyObj(M2!.transitions[1][1])!.group;
    gens1 := GeneratorsOfGroup(g1);
    gens2 := GeneratorsOfGroup(g2);
    src := [];
    dst := [];
    w := [];

    for i in [1..Length(gens1)] do
        Assert(0,Output(M1,i)=Output(M2,i));
        for c in Cycles(PermList(Output(M1,i)),AlphabetOfFRObject(M1)) do
            j := Product(M1!.transitions[i]{c});
            if IsOne(j) then continue; fi;
            k := Product(M2!.transitions[i]{c});
            Add(src,j);
            Add(dst,k);
            w[Position(gens1,CyclicallyReducedWord(j))] := Position(gens2,CyclicallyReducedWord(k));
        od;
    od;
    Assert(0,BoundPositions(w)=[1..Length(gens1)]);
    for i in [1..Length(gens2)] do
        if not i in w then
            Add(dst,gens2[i]);
            Add(src,One(g1));
        fi;
    od;
    return [GroupHomomorphismByImages(g2,g1,dst,src),w];
end);

InstallMethod(IsNonContractingMatrix, "(IMG) for an integer matrix",
        [IsMatrix],
        function(mat)
    local i, j, len, part, parts, submat, sublen;
    # find if there's an eigenvalue >= 1, without using numerical methods

    len := Length(mat);
    parts := List(EquivalenceClasses(StronglyConnectedComponents(BinaryRelationOnPoints(List([1..len],x->Filtered([1..len],y->IsPosRat(mat[x][y])))))),Elements);

    for part in parts do
        submat := mat{part}{part};
        sublen := Length(part);

        if NullspaceMat(submat-IdentityMat(sublen))=[] then # no 1 eigenval
            i := List([1..sublen],i->1);
            j := List([1..sublen],i->1);   # first approximation to perron-frobenius vector
            repeat
                i := i*submat;
                j := j*submat*submat; # j should have all entries growing exponentially
                if ForAll([1..sublen],a->j[a]=0 or j[a]<i[a]) then
                    return false; # perron-frobenius eigenval < 1
                fi;
            until ForAll(j-i,IsPosRat);
        fi;
    od;
    return true;
end);

#!!! still old format
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

InstallMethod(LiftOfConjugacyClass, "(IMG) for a machine and a conjugacy class",
        [IsGroupFRMachine,IsConjugacyClassGroupRep],
        function(M,c)
    local r, lift;
    r := DecompositionOfFRElement(Representative(c));
    lift := [];
    for c in Cycles(PermList(r[2]),AlphabetOfFRObject(M)) do
        Add(lift,[Length(c),ConjugacyClass(StateSet(M),Product(r[1]{c}))]);
    od;
    return lift;
end);

InstallMethod(ThurstonMatrix, "(IMG) for a machine and a list of sphere conj. classes",
        [IsSphereMachine,IsMulticurve],
        function(M,multicurve)
    # enlarge multicurve so that it becomes invariant, and return the transition matrix;
    # or return fail if the enlargement would cause curves to intersect
    local g, len, mat, lift, c, row,
           w, x, i, j, d, group, pi, gens, peripheral;

    g := StateSet(M);
    while g<>FamilyObj(Representative(multicurve[1]))!.group do
        Error("Elements do not all have the same underlying sphere group");
    od;

    len := Length(multicurve);
    mat := [];
    for c in multicurve do
        row := List([1..len],i->0);
        for lift in LiftOfConjugacyClass(M,c) do
            if IsPeripheral(lift[2]) then
                continue;
            fi;
            j := Position(multicurve,lift[2]);
            if j=fail then
                j := Position(multicurve,lift[2]^-1);
            fi;
            if j=fail then # add one more curve
                if ForAny(multicurve,c->IntersectionNumber(c,lift[2])>0) then
                    return fail;
                fi;
                for j in mat do Add(j,0); od;
                Add(row,1/lift[1]);
                len := len+1;
                Add(multicurve,c);
            else
                row[j] := row[j] + 1/lift[1];
            fi;
        od;
        Add(mat,row);
    od;

    return mat;
end);

InstallMethod(ThurstonObstruction, "(IMG) for a marked sphere and a machine",
        [IsSphereMachine,IsMarkedSphere],
        function(M,spider)
    # check if <spider> has coalesced points; in that case, find a short loop, saturate it by
    # taking its lifts, and check if they form an obstruction
    local multicurve, e, g, mat;

    #1. for each edge of the spanning tree, consider its conjugacy class c. This represents
    #   a simple closed curve, and if the edge is long then it probably is part of a long
    #   tube, and therefore its dual curve is short (we wouldn't travel twice along the
    #   long tube, in a minimal spanning tree; therefore, one cut is enough).

    #2. for each such class c, saturate it by taking preimages under M, checking all the time
    #   that the curves are all non-intersecting.

    #3. compute the Thurston matrix of the invariant multicurve, and its spectral radius.

    g := StateSet(M);
    while g<>spider!.model do
        Error("ThurstonObstruction: the marked sphere is not marked by the stateset of machine");
    od;

    for e in GeneratorsOfGroup(spider!.group) do
        multicurve := [ConjugacyClass(g,e^spider!.marking)];
        mat := ThurstonMatrix(M,multicurve);
        if mat=fail then continue; fi;
        if IsNonContractingMatrix(mat) then
#            if spider<>fail and IsBound(boundary[j]) then
#                if not IsBound(spider!.arcs) then spider!.arcs := []; fi;
#                Add(spider!.arcs,[[0,0,255],Float(105/100),boundary[j][1]]);
#            fi;
            return rec(multicurve := multicurve,
                       matrix := mat,
                       # boundary := SURROUNDINGCURVE(spider!.cut,"half of tree cut at e"),
                       machine := M);
        fi;
    od;

    return fail;
end);

InstallGlobalFunction(NormalizedQuadraticP1Map, function(f,M,param)
    # param is IsPolynomial, IsBicritical, or a positive integer.
    # in the first case, normalize f as z^d+a_{d-2}*z^{d-2}+...+a_0
    # in the second case, normalize f as (az^d+b)/(cz^d+e)
    # in the third case, normalize f as 1+a/z+b/z^2, such that 0 is on
    # a cycle of length <param>
    # return [new map, Möbius transformation from old to new]
    local p, i, j, k, a, b, mobius, m, coeff, degree, numer, denom;
    p := POSTCRITICALPOINTS@(f);
    degree := DegreeOfP1Map(f);
    
    if param=fail then
        if IsPolynomialSphereMachine(M) then
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
                return [P1Monomial(degree),mobius];
            else
                return [P1Monomial(-degree),mobius];
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
            return [P1Monomial(degree),mobius];
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
            k := SupportingRays(SphereMachine(ConjugatedP1Map(f,mobius)));
            if k=j then break; fi; #!!! maybe should be: "equivalent"rays?
            mobius := CompositionP1Map(mobius,m);
            Info(InfoIMG,1,"param:=IsPolynomial: trying to rotate around infinity");
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
                return [P1Monomial(degree),mobius];
            fi;
            j := First(p[4],r->r[1]=i-3)[2];
            mobius := MoebiusMap(p[2][3-i][1],p[3][j],p[2][i][1]);
        elif param=2 then # normalize as a/(z^2+2z), infty=>0->infty, -1=>
            if Length(p[3])=2 then # special case infty=>0=>infty or polynomial
                if First(p[4],r->r[1]=j)[2]=j then
                    return [P1Monomial(degree),mobius];
                else
                    return [P1Monomial(-degree),mobius];
                fi;
            fi;
            mobius := CompositionP1Map(MoebiusMap(p[3][j],p[2][3-i][1],p[2][i][1]),-P1z);
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

# we rewrite it, because the transition entries are not in the stateset.
# !!!! we'll have to generalize machines to have different input and output
BindGlobal("MYCHANGEFRMACHINEBASIS@", function(M,p)
    local trans, i, d, newM;
    d := Size(AlphabetOfFRObject(M));
    while LargestMovedPoint(p)>d do
	Error("Invalid permutation ",p,"\n");
    od;
    trans := [];
    for i in [1..Length(M!.transitions)] do
        Add(trans,Permuted(M!.transitions[i],p));
    od;
    newM := FRMachineNC(FamilyObj(M),StateSet(M),trans,List(M!.output,r->ListPerm(PermList(r)^p,d)));
    return newM;
end);

InstallMethod(ThurstonAlgorithm, "(IMG) for a sphere machine",
        [IsSphereMachine],
        function(M)
    local old_downsphere, downsphere, upsphere, n, deg, model,
          f, mobius, Mf, Mmobius, hom_mobius, match, v, i, j,
          dist, obstruction, lifts, sublifts, fast, poly;

    model := StateSet(M);
    n := Length(GeneratorsOfGroup(model));
    deg := Length(AlphabetOfFRObject(M));
    poly := IsPolynomialSphereMachine(M);
    
    if n=2 then # special handling, space is not hyperbolic
        i := Sum(List(Transitions(M,1),ExponentSums));
        if i[1]-i[2]=1 then
            return P1Monomial(deg);
        elif i[1]-i[2]=-1 then
            return P1Monomial(-deg);
        else
            Error(M," is not an IMG machine");
        fi;
    fi;    
    
    # create downsphere on +/- equidistributed points on Greenwich meridian.
    v := [];
    for i in [1..n] do
        # on positive real axis, tending to infinity
        f := Exp(@.2ipi*(i-1)/(2*n-2));
        Add(v,P1Point(@.o*ImaginaryPart(f)/(@.ro+RealPart(f))));
    od;
    downsphere := NewMarkedSphere(v,model);

    if ValueOption("julia")<>fail then
        i := ValueOption("julia");
        if not IsInt(i) then i := 1000; fi; # number of points to trace
        downsphere!.points := EquidistributedP1Points(i);
    fi;
    
    lifts := fail;
    f := fail; # in the beginning, we don't know them
    fast := false;

    repeat
        old_downsphere := downsphere;
        # find a rational map that has the right critical values
        f := BranchedCoveringByMonodromy(downsphere,Output(M)); #!!! use f, lifts
        lifts := Concatenation(List(f.points,x->List(x,y->y[1])));
        f := f.map;
        Info(InfoIMG,3,"1: found rational map ",f," on vertices ",lifts);

        if fast then # just get points closest to those in downsphere t
            match := MatchP1Points(sublifts,List(sublifts,x->lifts));
            if match=fail then
		Info(InfoIMG,3,"1.5: back to slow mode");
                fast := false; continue;
            fi;
        else
            # create a spider on the full preimage of the points of <downsphere>
            upsphere := NewMarkedSphere(lifts);
            Info(InfoIMG,3,"2: created liftedspider ",upsphere);

            # lift paths in <downsphere> to <upsphere>
            Mf := SphereMachineOfBranchedCovering(downsphere,upsphere,f,poly);
            Info(InfoIMG,3,"3: recursion ",Mf);
            
            # find a bijection between the alphabets of <Mf> and <M>
            match := MATCHPERMS@(M,Mf!.output);
            if match=fail then return fail; fi;
            
            Mf := MYCHANGEFRMACHINEBASIS@(Mf,match);
            Info(InfoIMG,3,"4: alphabet permutation ",match);

            # extract those vertices in <v> that appear in the recursion
            match := MATCHTRANS@(M,Mf);
            Info(InfoIMG,3,"5: extracted and sorted vertices ",match);
            upsphere!.marking := upsphere!.marking*match[1]; # so we're again marked by <model>
            upsphere!.model := model;
            match := match[2];
        fi;

        # find a mobius transformation that normalizes <sublifts> wrt PSL2C
        sublifts := lifts{match};
        mobius := P1MapNormalizingP1Points(sublifts,VerticesOfMarkedSphere(downsphere));
        Info(InfoIMG,3,"6: normalize by mobius map ",mobius);

        # now create the new spider on the image of these points
        v := List(sublifts,p->P1Image(mobius,p));
        
	if fast then
            dist := Sum([1..Length(v)],i->P1Distance(VerticesOfMarkedSphere(downsphere)[i],v[i]));
            if dist>@.fast*Length(v) then
                fast := false;
                Info(InfoIMG,3,"7: legs moved ",dist,"; back to slow mode");
                downsphere := old_downsphere;
                continue; # restart
            fi;
            # just wiggle spider around
            downsphere := WiggledMarkedSphere(downsphere,v);
        else
            downsphere := NewMarkedSphere(v,model);
            Mmobius := SphereMachineOfBranchedCovering(downsphere,upsphere,mobius,poly);
            Info(InfoIMG,3,"7: new spider ",downsphere," with recursion ",Mmobius);

            # compose recursion of f with that of mobius
            hom_mobius := GroupHomomorphismByImages(model,model,GeneratorsOfGroup(model),List(Mmobius!.transitions,x->x[1]));
            Info(InfoIMG,3,"8: composed recursion is ",hom_mobius);
                   
            # finally set marking of new spider using M
            downsphere!.marking := downsphere!.marking*hom_mobius;
            Info(InfoIMG,3,"9: marked new spider ",downsphere);
        fi;
        
        dist := DistanceMarkedSpheres(downsphere,old_downsphere,fast);
        Info(InfoIMG,2,"Spider moved ",dist," steps; feet=",VerticesOfMarkedSphere(downsphere)," marking=",downsphere!.marking);

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
        obstruction := ThurstonObstruction(M,downsphere);
        if obstruction<>fail then
            return obstruction;
        fi;
    until false;
    
    Info(InfoIMG,2,"Spider converged");
    
    i := NormalizedQuadraticP1Map(f,M,ValueOption("param"));
    if i<>fail then
        downsphere := WiggledMarkedSphere(downsphere,InverseP1Map(i[2]));
        f := i[1];
    fi;

    # find a good twist to shorten the marking; presumably, this
    # will give a simpler machine
    v := SIMPLIFYBYBRAIDTWISTS@(model, downsphere!.marking);
    for i in v do M := M^i; od;
    SetCorrespondence(M,v);
    
    downsphere!.map := f;
    downsphere!.cycle := ATTRACTINGCYCLES@(POSTCRITICALPOINTS@(f));

    return rec(map := f, machine := M, markedsphere := downsphere);
end);

InstallMethod(P1MapBySphereMachine, "(IMG) for a sphere machine",
        [IsSphereMachine],
        function(M)
    local a;
    a := ThurstonAlgorithm(M);
    if IsBound(a.map) then
        return a.map;
    else
        return a; # the obstruction
    fi;
end);

#E thurston.gi . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
