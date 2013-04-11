#############################################################################
##
#W hurwitz.gi                                               Laurent Bartholdi
##
#Y Copyright (C) 2011-2013, Laurent Bartholdi
##
#############################################################################
##
##  Solving the Hurwitz problem
##
#############################################################################

BindGlobal("LIFTBYMONODROMY@", function(spider,monodromy,d)
    # create a new triangulation lifted by the given monodromy rep'n.
    # the positions in the triangulation are the same as in the original
    # spider; i.e., geometrically one takes deg(monodromy) copies of the
    # original sphere, cuts it along the minimal spanning tree, and glues
    # the sheets to each other along the monodromy rep'n.
    #
    # some vertices of the result acquire extra fields, "degree", which is
    # the local degree of the map, and "cover", which points to the vertex
    # of "spider" that is being covered.

    local e, f, v, i, j, left, from, edgeperm, reverse,
          edges, faces, lift;
    
    # make d copies of the faces and edges, renumber their indices
    faces := List([1..d],i->List(spider!.cut!.f,f->Objectify(TYPE_FACE,
                     rec(index := f!.index+(i-1)*Length(spider!.cut!.f), sheet := i, below := f))));
    edges := List([1..d],i->List(spider!.cut!.e,e->Objectify(TYPE_EDGE,
                     rec(index := e!.index+(i-1)*Length(spider!.cut!.e), sheet := i, below := e))));

    # reattach the faces
    for f in spider!.cut!.f do
        j := Neighbour(f)!.index;
        for i in [1..d] do
            SetNeighbour(faces[i][f!.index],edges[i][j]);
        od;
    od;

    # reattach the edges according to the monodromy
    for e in spider!.cut!.e do
        j := Left(e)!.index;
        for i in [1..d] do
            SetLeft(edges[i][e!.index],faces[i][j]);
        od;
        j := Next(e)!.index;
        for i in [1..d] do
            SetNext(edges[i][e!.index],edges[i][j]);
        od;
        j := Prevopp(e)!.index;
        edgeperm := (GroupElement(Opposite(Prevopp(e)))^spider!.marking)^monodromy;
        for i in [1..d] do
            SetPrevopp(edges[i][e!.index],edges[i^edgeperm][j]);
        od;
    od;

    # create the triangulation, except for the vertices
    lift := rec(v := [],
                e := Concatenation(edges),
                f := Concatenation(faces));

    # create new vertices, attach the edges to them
    for e in lift.e do
        if not HasFrom(e) then
            from := From(e!.below);
            v := Objectify(TYPE_VERTEX, rec(index := Length(lift.v)+1, below := from));
            SetNeighbour(v,e);
            SetFrom(e,v);
            v!.degree := Valency(v) / Valency(from);
            v!.sheet := List(Neighbours(v),x->x!.sheet);
            SetPos(v,Pos(from));
            SetIsFake(v,IsFake(from));
            for e in Neighbours(v) do
                SetFrom(e,v);
            od;
            Add(lift.v,v);
        fi;
    od;
    #!TODO
    # we'd like to keep track of the cycle corresponding
    # to v, not just its length; for this, we have to
    # consider the first edge (in the boundarytree of spider)
    # that starts at v, and keep track of the sheets
    # above that edge. Since we subdivided spider, we lost
    # the relation between spider's boundarytree and what
    # we cover.

    return Objectify(TYPE_TRIANGULATION,lift);
end);

BindGlobal("REFINETRIANGULATION@", function(triangulation,maxlen)
    # refines triangulation by adding circumcenters, until
    # the length of every edge (say connecting v to w) is at most
    # maxlen^Maximum(v.degree*w.degree)
    local idle, e, f, maxdegree, len, mult, p;

    for e in triangulation!.e do
        if From(e)!.degree>1 and To(e)!.degree>1 then
            AddToTriangulation(triangulation,Left(e),Pos(e));
            triangulation!.v[Length(triangulation!.v)]!.degree := 1;
            SetIsFake(triangulation!.v[Length(triangulation!.v)],true);
        fi;
    od;
    
    maxdegree := Maximum(List(triangulation!.v,v->v!.degree));
    mult := Int(Log(7./maxlen)/Log(1.5));
    repeat
        len := List([1..maxdegree],i->(maxlen*1.5^mult)^i);
        idle := true;
        for e in triangulation!.e do
            if Length(e) > len[Maximum(From(e)!.degree,To(e)!.degree)] then
                idle := false;
                f := Left(e);
                AddToTriangulation(triangulation,f,Centre(f));
                triangulation!.v[Length(triangulation!.v)]!.degree := 1;
                SetIsFake(triangulation!.v[Length(triangulation!.v)],true);
            fi;
        od;
        if idle then mult := mult-1; fi;
    until idle and mult<0;
end);

BindGlobal("LAYOUTTRIANGULATION@", function(triangulation)
    # run C code to optimize point placement.
    # "triangulation" is topologically a triangulated sphere. Its edge
    # lengths should be adjusted, by multiplying each edge (say from v to w)
    # by u[v]*u[w] for some scaling function u defined on the vertices;
    # in such a way that the resulting metric object is a conformal sphere.
    local i, e, f, v, m, map, max, infty, sin, sout, stdin, stdout;

    sin := "";
    stdin := OutputTextString(sin,false);
    
    max := 0;
    for v in triangulation!.v do
        m := Valency(v);
        if m>max then
            infty := v!.index;
            max := m;
        fi;
    od;
    PrintTo(stdin,"VERTICES ",Length(triangulation!.v),"\n",infty,"\n");
    
    PrintTo(stdin,"FACES ",Length(triangulation!.f),"\n");
    for f in triangulation!.f do
        for e in Neighbours(f) do
            PrintTo(stdin,From(e)!.index," ");
        od;
        for e in Neighbours(f) do
            e := Next(e);
            PrintTo(stdin,Length(e)^(1/Maximum(To(e)!.degree,From(e)!.degree))," ");
        od;
        PrintTo(stdin,"\n");
    od;
    PrintTo(stdin,"END\n");
    CloseStream(stdin);
PrintTo("layout-in",sin);#!!!
    stdin := InputTextString(sin);
    sout := "";
    stdout := OutputTextString(sout,false);
    Process(DirectoryCurrent(),Filename(DirectoriesPackagePrograms("img"),"layout"),stdin,stdout,[]);
    CloseStream(stdin);
    CloseStream(stdout);

PrintTo("layout-out",sout);#!!!
    m := EvalString(sout);
    m := 1.0*m; # make sure all entries are floats
    v := List(m,P1Sphere);
    map := P1MapNormalizingP1Points(v);
    v := List(v,v->P1Image(map,v));
    
    for i in [1..Length(v)] do
        triangulation!.v[i]!.Pos := v[i];
    od;
    for e in triangulation!.e do # reset the length and position of the edge
        RESETEDGE@(e);
    od;
    for f in triangulation!.f do
        RESETFACE@(f);
    od;
end);

BindGlobal("OPTIMIZELAYOUT@", function(spider,lift)
    # "lift" is an approximate lift of "spider". refine its positions.
    
    # find Möbius transformation putting the most ramified c.v. at 0,infty,
    # and the next most ramified c.v. at 1;
    # and find another Möbius transformation so that 0,infty,1 are fixed,
    # all of maximal order.
    local cp, max, l0, l1, linf, pre, post, x, y, numzero, numcv,
          sin, sout, stdin, stdout, data, cv, i, printp1, scanp1;
    
    printp1 := function(stream,z)
        z := P1Coordinate(z);
        PrintTo(stream,"(", RealPart(z),",",ImaginaryPart(z),")");
    end;

    scanp1 := function(str)
        str := SplitString(str{[2..Length(str)-1]},",");
        return STRINGS2P1POINT(str[1],str[2]);
    end;
    
    sin := "";
    stdin := OutputTextString(sin,false);

    cp := Filtered(lift!.v,x->IsBound(x!.degree) and IsBound(x!.below));
    max := Maximum(List(cp,x->x!.degree));
    linf := First(cp,x->x!.degree=max);
    max := Maximum(List(Filtered(cp,x->x!.below<>linf!.below),x->x!.degree));
    l0 := First(cp,x->x!.below<>linf!.below and x!.degree=max);
    max := Maximum(List(Filtered(cp,x->x!.below<>linf!.below and x!.below<>l0!.below),x->x!.degree));
    l1 := First(cp,x->x!.below<>linf!.below and x!.below<>l0!.below and x!.degree=max);
    
    post := InverseP1Map(MoebiusMap(List([l0,l1,linf],x->Pos(x!.below))));
    pre := InverseP1Map(MoebiusMap(List([l0,l1,linf],Pos)));
    
    PrintTo(stdin,"DEGREE ",(Sum(cp,x->x!.degree-1)+2)/2,"\n");
    numzero := Number(cp,x->x!.below=l0!.below or x!.below=linf!.below)-2;
    PrintTo(stdin,"ZEROS/POLES ",numzero,"\n");
    for x in cp do
        if (x!.below=l0!.below or x!.below=linf!.below) and x<>l0 and x<>linf then
            if x!.below=l0!.below then
                PrintTo(stdin,x!.degree); # zero
            else
                PrintTo(stdin,-x!.degree); # pole
            fi;
            PrintTo(stdin," "); printp1(stdin,Pos(x)^pre); PrintTo(stdin,"\n");
        fi;
    od;
    PrintTo(stdin,l0!.degree,"\n");
    numcv := Number(cp,x->x!.below<>l0!.below and x!.below<>linf!.below and x!.degree>1)-1;
    cv := [];
    PrintTo(stdin,"CRITICAL ",numcv,"\n");
    for x in cp do
        if x!.below<>l0!.below and x!.below<>linf!.below and x!.degree>1 and x<>l1 then
            PrintTo(stdin,x!.degree);
            PrintTo(stdin," "); printp1(stdin,Pos(x)^pre);
            PrintTo(stdin," "); printp1(stdin,Pos(x!.below)^post);
            PrintTo(stdin,"\n");
            Add(cv,x!.below);
        fi;
    od;
    Add(cv,l1!.below);
    PrintTo(stdin,l1!.degree,"\n");
    PrintTo(stdin,"END\n");
    CloseStream(stdin);
    
    stdin := InputTextString(sin);
    sout := "";
    stdout := OutputTextString(sout,false);
    
    Process(DirectoryCurrent(),Filename(DirectoriesPackagePrograms("img"),"hsolve"),stdin,stdout,[]);
    CloseStream(stdin);
    CloseStream(stdout);
    
    sout := SplitString(sout,WHITESPACE);
    Assert(0,sout[1]="DEGREE" and Int(sout[2])=(Sum(cp,x->x!.degree-1)+2)/2);
    data := rec(degree := Int(sout[2]),
                zeros := [],
                poles := [],
                cp := [],
                post := post);
    Assert(0,sout[3]="ZEROS/POLES" and numzero=Int(sout[4]));
    for i in [0..numzero] do
        x := Int(sout[5+2*i]);
        if i=numzero then
            y := P1zero;
        else
            y := scanp1(sout[6+2*i]);
        fi;
        if x>0 then
            Add(data.zeros,rec(degree := x, pos := y, to := l0!.below));
        else
            Add(data.poles,rec(degree := -x, pos := y, to := linf!.below));
        fi;
    od;
    Assert(0,sout[2*numzero+6]="CRITICAL" and Int(sout[2*numzero+7])=numcv);
    for i in [0..numcv] do
        x := Int(sout[2*numzero+8+3*i]);
        if i=numcv then
            y := P1one;
        else
            y := scanp1(sout[2*numzero+9+3*i]);
        fi;
        Add(data.cp,rec(degree := x, pos := y, to := cv[i+1]));
    od;
    if numcv >= 0 then
        Assert(0,sout[2*numzero+3*numcv+9]="END");
    fi;
    Add(data.poles,rec(degree := 2*data.degree-1-Sum(Concatenation(data.cp,data.zeros,data.poles),x->x.degree-1), pos := P1infinity, to := linf!.below));
    
    return data;
end);

InstallMethod(HurwitzMap, "(IMG) for a spider and a homomorphism",
        [IsMarkedSphere,IsGroupHomomorphism],
        function(spider,monodromy)
    # compute the critical points, zeros and poles of a map whose
    # critical values are vertices of "spider", with monodromy given
    # by the homomorphism "monodromy".
    local t, d;

    Assert(0,Range(spider!.marking)=Source(monodromy));
    
    d := Maximum(LargestMovedPoint(Range(monodromy)),1);
    Assert(0,IsTransitive(Image(monodromy),[1..d]));

    t := LIFTBYMONODROMY@(spider,monodromy,d);
    REFINETRIANGULATION@(t,@.hurwitzmesh);
    LAYOUTTRIANGULATION@(t);
    d := OPTIMIZELAYOUT@(spider,t);
    
    d.map := P1MapByZerosPoles(Concatenation(List(d.zeros,x->ListWithIdenticalEntries(x.degree,x.pos))),
                     Concatenation(List(d.poles,x->ListWithIdenticalEntries(x.degree,x.pos))),
                     
                     P1one,P1one);

    return d;
end);

InstallMethod(DessinByPermutations, "(IMG) for three permutations",
        [IsPerm,IsPerm,IsPerm],
        function(s0,s1,sinf)
    local permrep, f, g, spider, d, i, above1, p, dist, distmin;
    
    f := SphereGroup([1,2,3]);
    g := Group(s0,s1,sinf);
    permrep := GroupHomomorphismByImages(f,g,GeneratorsOfGroup(f),GeneratorsOfGroup(g));
    
    spider := NewMarkedSphere([P1zero,P1one,P1infinity],f);    
    d := HurwitzMap(spider,permrep);
    
    for i in d.zeros do Unbind(i.to); od;
    for i in d.poles do Unbind(i.to); od;
    d.above1 := [];
    for i in d.cp do Add(d.above1,rec(degree := i.degree, pos := i.pos)); od;
    above1 := P1PreImages(d.map,P1PreImages(d.post,P1one)[1]);
    
    for p in Concatenation(List(d.cp,x->ListWithIdenticalEntries(x.degree,x.pos))) do
        dist := List(above1,q->P1Distance(p,q));
        distmin := Minimum(dist);
        i := Position(dist,distmin);
        Remove(above1,i);
    od;
    for p in above1 do
        Add(d.above1, rec(degree := 1, pos := p));
    od;
    
    Unbind(d.cp);
    
    return d;
end);

InstallMethod(DessinByPermutations, "(IMG) for two permutations",
        [IsPerm,IsPerm],
        function(s0,s1)
    return DessinByPermutations(s0,s1,(s0*s1)^-1);
end);

BindGlobal("GENERALHURWITZMAP@", function(z,spider,perm,oldf,oldlifts)
    # returns [f,full preimage of Vertices(spider)]
    local d, r, v, pre, new, old, i, dist, mindist, permrep;
    
    #!! should create a new spider with only the critical values, not the
    # whole post-critical set.
    
    permrep := GroupHomomorphismByImages(Source(spider!.marking),
                       SymmetricGroup(Length(perm[1])),
                       GeneratorsOfGroup(Source(spider!.marking)),
                       List(perm,PermList));
    d := HurwitzMap(spider,permrep);
    
    pre := [];
    for v in spider!.cut!.v do
        if IsFake(v) then continue; fi;
        new := P1PreImages(d.map,P1PreImages(d.post,v.pos)[1]);
        old := Filtered(Concatenation(d.cp,d.poles,d.zeros),r->r.to=v);
        for r in old do
            Add(pre,r.pos);
            for i in [1..r.degree] do
                dist := List(new,p->P1Distance(r.pos,p));
                mindist := Minimum(dist);
                Remove(new,Position(dist,mindist));
            od;
        od;
        Append(pre,new);
    od;
    
    return [CompositionP1Map(d.map,d.post),pre];
end);

BindGlobal("TRICRITICAL@", function(perm)
    # find a rational function with critical values 0,1,infinity
    # with monodromy actions perm[1],perm[2],perm[3]
    # return fail if it's too hard to do;
    # otherwise, return [map, critical points (on sphere),order],
    # where order is a permutation of the critical values:
    # ELM_LIST([0,1,infinity],order[i]) has permutation perm[i]
    
    # the cases covered are:
    # [[a],[b],[c]], degree=(a+b+c-1)/2
    # [[m,n],[m,n],[3]], degree=m+n
    # [[2,3],[2,3],[2,2]], degree=5
    # [[n,n],[2,...,2],[2,...,2]], degree=2n
    # [[degree],[m,degree-m+1]]
    # [[degree],[m,degree-m],[2]]
    local deg, cl, i, j, k, m, points, f, order, p;
    
    deg := Length(perm[1]);
    perm := List(perm,PermList);
    cl := List(perm,x->SortedList(CycleLengths(x,[1..deg])));
    
    points := [P1zero, P1one, P1infinity];
    
    if ForAll(cl,x->Length(DifferenceLists(x,[1]))=1) then # [[a],[b],[c]]
        cl := List(cl,x->DifferenceLists(x,[1])[1]);

        m := List([0..deg-cl[2]],row->List([0..deg],col->(-1)^(col-row)*Binomial(cl[2],col-row)));
        p := NullspaceMat(m{1+[0..deg-cl[2]]}{1+[deg-cl[3]+1..cl[1]-1]})[1];
        p := [,p*Lcm(List(p,DenominatorRat))];
        j := p[2]*m;
        j := j / Gcd(j);
        p[3] := -j{1+[0..deg-cl[3]]};
        p[1] := j{1+[cl[1]..deg]};
        for j in [1..3] do
            p[j] := P1MapByCoefficients(@.o*p[j]{[1..deg+1-cl[j]]});
        od;
        f := P1Monomial(cl[1])*p[1]/p[3];
        for j in [1..3] do
            Append(points,List(RootsFloat(p[j]),P1Point));
        od;
        return [f,points,[1,2,3]];
    fi;
    
    if Size(Set(cl))<=2 and ForAll(cl,x->DifferenceLists(x,[1])=[3] or Length(x)=2) then # [m+n,m+n,3]
        i := PositionProperty(cl,x->DifferenceLists(x,[1])=[3]);
        m := cl[1+(i mod 3)];
        f := P1z^m[2]*((m[1]-m[2])*P1z+(m[1]+m[2]))^m[1]/((m[1]+m[2])*P1z+(m[1]-m[2]))^m[1];
        Add(points,P1Point((m[1]+m[2])/(m[2]-m[1])));
        k := P1PreImages(f,P1one);
        SortParallel(List(k,x->P1Distance(x,P1one)),k);
        Append(points,k{[4..deg]});
        if i=2 then order := [1,2,3]; else order := ListPerm((i,2),3); fi;
        return [f,points,order];
    fi;
    
    if deg=5 and IsEqualSet(cl,[[2,3],[1,2,2]]) then # (1,2)(3,4,5),(1,3)(2,5,4),(1,5)(2,3)
        f := P1z^3*((4*P1z+5)/(5*P1z+4))^2;
        Add(points,P1Point(-4/5)); # to infinity
        Add(points,P1Point(-5/4)); # to 0
        Add(points,P1Point((-7+Sqrt(-15*@.o))/8)); # to 1
        Add(points,P1Point((-7-Sqrt(-15*@.o))/8)); # to 1
        order := Permuted([1,3,2],(1,2,3)^Position(cl,[1,2,2]));
        return [f,points,order];
    fi;
    
    i := First([1..3],i->cl[i]=[deg/2,deg/2]);
    if i<>fail and ForAll([1..3],j->i=j or Set(cl[j])=[2]) then
        # deg = 2n; shapes [n,n],[2,...,2],[2,...,2]
        f := 4*P1z^(deg/2)/(1+P1z^(deg/2))^2;
        order := Permuted([3,2,1],(1,2,3)^i);
        Remove(points,2); # remove 1
        Append(points,List([0..deg-1],i->P1Point(Exp(@.2ipi*i/deg))));
        return [f,points,order];
    fi;
    
    i := First([1..3],i->cl[i]=[deg]); # max. cycle
    if i=fail then return fail; fi; # now only accept polynomials
    if Product(perm)=() then
        j := i mod 3+1; k := j mod 3+1;
    else
        k := i mod 3+1; j := k mod 3+1;
    fi;
    
    m := First([j,k],i->Length(cl[i])=2);
    if m<>fail then # [d],[m,d-m], [2,1,...,1]
        order := [m,j+k-m,i];
        m := cl[m][1];
        f := (P1z*deg/m)^m*((1-z)*deg/(deg-m))^(deg-m);
        points := points{order};
        Add(points,P1Point(m/deg));
        i := P1PreImages(f,P1one);
        SortParallel(List(i,z->P1Distance(z,P1Point(m/deg))),i);
        Append(points,i{[3..deg]});
        return [f,points,order];
    fi;
    
    m := Maximum(cl[j]);
    if Set(cl[j])=[1,m] and Set(cl[k])=[1,deg-m+1] then
        # so we know the action around i is (1,...,deg), at infinity
        # the action around j is (m,m-1,...,1), at 0
        # the action around k in (deg,deg-1...,m), at 1
        f := m*Binomial(deg,m)*Primitive(P1z^(m-1)*(1-P1z)^(deg-m));
        order := [j,k,i];
        points := points{order};
        for i in [P1zero,P1one] do
            j := P1PreImages(f,i);
            k := List(j,x->P1Distance(x,i));
            SortParallel(k,j);
            if i=P1zero then
                j := j{[m+1..deg]};
            else
                j := j{[deg+2-m..deg]};
            fi;
            Append(points,j);
        od;
        return [f,points,order];
    fi;
    
    return fail;
end);

BindGlobal("QUADRICRITICAL@", function(perm,values)
    local c, w, f, m, id, aut, z;
    
    # normalize values to be 0,1,infty,w
    aut := MoebiusMap(values{[1..3]});
    w := P1Coordinate(P1Image(InverseP1Map(aut),values[4]));
    
    # which two values have same deck transformation?
    id := First(Combinations(4,2),p->perm[p[1]]=perm[p[2]]);
    
    z := Indeterminate(@.isc);
    c := RootsFloat((z-2)^3*z-w*(z+1)^3*(z-1));
    
    # find appropriate c
    f := List(c,c->P1z^2*(c*(P1z-1)+2-c)/(c*(P1z+1)-c));
    m := List(f,FRMachine);
    
    f := CompositionP1Map(aut,f[First([1..Length(c)],i->Output(m[i],id[1])=Output(m[i],id[2]))]);
    
    return [f, [P1zero, P1one, P1infinity, P1Point(c*(c-2)/(c^2-1))]];
end);

InstallMethod(LiftMarkedSphereByMonodromy, "(FR) for a marked sphere",
        [IsMarkedSphere,IsGroupHomomorphism],
        function(spider,monodromy)
    # find a rational map that has critical values at <feet(spider)>, with
    # monodromy action given by <monodromy>.
    # returns [map,lifted_sphere] where <lifted_sphere> has as vertex set the full preimage
    # of <values>

    #!!! that's really HURWITZMAP + RETURN THE UPPER SPHERE
    local cv, values, p, f, points, deg, i;
    
    values := Vertices(spider);
    cv := Filtered([1..Length(values)],i->not ISONE@FR(perm[i]));
    deg := Length(perm[1]);
    if Length(cv)=2 then # bicritical
        f := CompositionP1Map(MoebiusMap(values{cv}),P1Monomial(deg));
        points := [P1zero,P1infinity];
    elif Length(cv)=3 then # tricritical
        p := TRICRITICAL@(perm{cv});
        if p<>fail then
            f := CompositionP1Map(MoebiusMap(ELMS_LIST(values{cv},p[3])),p[1]);
            points := p[2];
        fi;
    elif deg=3 then # quadricritical, but degree 3
        p := QUADRICRITICAL@(perm{cv},values{cv});
        f := p[1];
        points := p[2];
    fi;

    if not IsBound(points) then # run hurwitz specialized code
        return GENERALHURWITZMAP@(spider,perm,oldf,oldlifts);
    fi;

    for i in [1..Length(values)] do if not i in cv then
        Append(points,P1PreImages(f,values[i]));
    fi; od;

    return [f,points];
end);

#!!! call it "LiftSpiderByPermutations" or something
    
#E hurwitz.gi . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
