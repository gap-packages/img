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

BindGlobal("TRICRITICAL@", function(deg,perm)
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
    local cl, i, j, k, m, p, data, order;
    
    cl := List(perm,x->SortedList(CycleLengths(x,[1..deg])));
    order := ();

    # first case: [[a],[b],[c]]
    if ForAll(cl,x->Length(DifferenceLists(x,[1]))=1) then
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
        data := rec(map := P1Monomial(cl[1])*p[1]/p[3],
                    points := [[[P1zero,cl[1]]],[[P1one,cl[2]]],[[P1infinity,cl[3]]]]);
        for j in [1..3] do
            Append(data.points[j],List(RootsFloat(p[j]),z->[P1Point(z),1]));
        od;

    # [[m,n],[m,n],[3]]
    elif Size(Set(cl))<=2 and ForAll(cl,x->DifferenceLists(x,[1])=[3] or Length(x)=2) then
        i := PositionProperty(cl,x->DifferenceLists(x,[1])=[3]);
        m := cl[1+(i mod 3)];
        data := rec(map := P1z^m[2]*((m[1]-m[2])*P1z+(m[1]+m[2]))^m[1]/((m[1]+m[2])*P1z+(m[1]-m[2]))^m[1],
                    points := [[[P1zero,m[2]],[P1Point((m[1]+m[2])/(m[2]-m[1])),m[1]]],
                            [[P1one,3]],
                            [[P1infinity,m[2]],[P1Point((m[2]-m[1])/(m[1]+m[2])),m[1]]]]);
        if i<>2 then order := (i,2); fi;

    # (1,2)(3,4,5),(1,3)(2,5,4),(1,5)(2,3)
    elif deg=5 and IsEqualSet(cl,[[2,3],[1,2,2]]) then
        data := rec (map := P1z^3*((4*P1z+5)/(5*P1z+4))^2,
                     points := [[[P1zero,3],[P1Point(-5/4),2]],
                             [[P1one,1],[P1Point((-7+Sqrt(-15*@.o))/8),2],[P1Point((-7-Sqrt(-15*@.o))/8),2]],
                             [[P1infinity,3],[P1Point(-4/5),2]]]);
        order := (1,2,3)^(Position(cl,[1,2,2])+1);
    else
        i := First([1..3],i->cl[i]=[deg/2,deg/2]);
        # deg = 2n; shapes [[n,n],[2,...,2],[2,...,2]]
        if i<>fail and ForAll([1..3],j->i=j or Set(cl[j])=[2]) then
            data := rec(map := 4*P1z^(deg/2)/(1+P1z^(deg/2))^2,
                        points := [[[P1zero,deg/2],[P1infinity,deg/2]],
                                List([0,2..deg-2],i->[P1Point(Exp(@.2ipi*i/deg)),2]),
                                List([1,3..deg-1],i->[P1Point(Exp(@.2ipi*i/deg)),2])]);
            order := (1,3)^((1,3,2)^i);
        else
            # find maximal cycle
            i := First([1..3],i->cl[i]=[deg]);
            if i=fail then
                return fail; # now only accept polynomials
            fi;

            if Product(perm)=() then
                j := i mod 3+1; k := j mod 3+1;
            else
                k := i mod 3+1; j := k mod 3+1;
            fi;
        
            m := First([j,k],i->Length(cl[i])=2);
            if m<>fail then # [d],[m,d-m], [2,1,...,1]
                order := PermList([m,j+k-m,i]);
                m := cl[m][1];
                data := rec(map := (P1z*deg/m)^m*((1-P1z)*deg/(deg-m))^(deg-m),
                            points := [[[P1zero,m],[P1one,deg-m]],
                                    [[P1Point(m/deg),2]],
                                    [[P1infinity,deg]]]);
            fi;
        fi;
    fi;

    if not IsBound(data) then
        return fail;
    fi;

    if order<>() then
        data.points := Permuted(data.points,order);
        data.map := CompositionP1Map(MoebiusMap(Permuted([P1zero,P1one,P1infinity],order^-1)),data.map);
    fi;

    return data;
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

InstallMethod(BranchedCoveringByMonodromy, "(IMG) for a spider and a homomorphism",
        [IsMarkedSphere,IsGroupHomomorphism],
        function(spider,monodromy)
    # compute the critical points, zeros and poles of a map whose
    # critical values are vertices of "spider", with monodromy given
    # by the homomorphism "monodromy".
    local t, d, g, gens, cv, values, data, i, j, k, dist, mindist, minpos, new, v, w;

    g := Source(monodromy);
    gens := GeneratorsOfGroup(g);
    Assert(0,g=Range(spider!.marking));
    
    d := Maximum(LargestMovedPoint(Range(monodromy)),1);
    Assert(0,IsTransitive(Image(monodromy),[1..d]));

    cv := Filtered([1..Length(gens)],i->gens[i]^monodromy<>());
    values := Vertices(spider);
    data := fail;

    if Length(cv)=2 then # bicritical
        data := rec(map := CompositionP1Map(MoebiusMap(values{cv}),P1Monomial(d)),
                    points := []);
        data.points{cv} := [[[P1zero,d]],[[P1infinity,d]]];
    elif Length(cv)=3 then # tricritical
        data := TRICRITICAL@(d,List(cv,i->gens[i]^monodromy));
        if data<>fail then
            data.map := CompositionP1Map(MoebiusMap(values{cv}),data.map);
        fi;
#    elif d=3 then # quadricritical, but degree 3
#        p := QUADRICRITICAL@(perm{cv},values{cv});
    fi;

    if data=fail then # general lifting procedure
        t := LIFTBYMONODROMY@(spider,monodromy,d);
        REFINETRIANGULATION@(t,@.hurwitzmesh);
        LAYOUTTRIANGULATION@(t);
        t := OPTIMIZELAYOUT@(spider,t);
    
        data := rec(map := CompositionP1Map(t.post,P1MapByZerosPoles(Concatenation(List(t.zeros,x->ListWithIdenticalEntries(x.degree,x.pos))),
                        Concatenation(List(t.poles,x->ListWithIdenticalEntries(x.degree,x.pos))),
                     
                        P1one,P1one)),
                    points := []);
        for v in Concatenation(t.cp,t.poles,t.zeros) do
            i := v.to!.index;
            if not IsBound(data.points[i]) then data.points[i] := []; fi;
            Add(data.points[i], [v.pos,v.degree]);
        od;
    fi;

    for v in spider!.cut!.v do
        if IsFake(v) then continue; fi;
        i := v!.index;
        if not IsBound(data.points[i]) then data.points[i] := []; fi;
        if Sum(data.points[i],x->x[2])=d then continue; fi;

        new := P1PreImages(data.map,Pos(v));
        for w in data.points[i] do
            for j in [1..w[2]] do
                # remove the point in new that is closest to w.pos
                mindist := 10*@.ro; # larger than sphere
                for k in [1..Length(new)] do
                    dist := P1Distance(new[k],w[1]);
                    if dist<mindist then mindist := dist; minpos := k; fi;
                od;
                Remove(new,minpos);
            od;
        od;
        Append(data.points[i],List(new,z->[z,1]));
    od;

    return data;
end);

InstallMethod(DessinByPermutations, "(IMG) for three permutations",
        [IsPerm,IsPerm,IsPerm],
        function(s0,s1,sinf)
    local permrep, f, g, spider, d, i, above1, p, dist, distmin;
    
    f := SphereGroup([1,2,3]);
    g := Group(s0,s1,sinf);
    permrep := GroupHomomorphismByImages(f,g,GeneratorsOfGroup(f),GeneratorsOfGroup(g));
    
    spider := NewMarkedSphere([P1zero,P1one,P1infinity],f);    
    return BranchedCoveringByMonodromy(spider,permrep);
end);

InstallMethod(DessinByPermutations, "(IMG) for two permutations",
        [IsPerm,IsPerm],
        function(s0,s1)
    return DessinByPermutations(s0,s1,(s0*s1)^-1);
end);
    
#E hurwitz.gi . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
