#############################################################################
##
#W markedsphere.gi                                          Laurent Bartholdi
##
#Y Copyright (C) 2011-2013, Laurent Bartholdi
##
#############################################################################
##
##  Points in Teichmuller space, as marked triangulated spheres, and their lifts
##
#############################################################################

##############################################################################
##
#M  MarkedSpheres
##
InstallMethod(ViewString, "(IMG) for a point in Teichmuller space",
        [IsMarkedSphere],
        s->Concatenation("<marked sphere on ",ViewString(s!.cut)," marked by ",String(s!.marking),">"));

InstallMethod(DisplayString, "(IMG) for a point in Teichmuller space",
        [IsMarkedSphere],
        s->CONCAT@FR(DisplayString(s!.cut),"Spanning tree on edges ",List(s!.treeedge,r->r!.index)," costing ",s!.treecost,"\nMarking ",s!.marking,"\n"));

InstallMethod(String, "(IMG) for a point in Teichmuller space",
        [IsMarkedSphere],
        s->Concatenation("Spider(",ViewString(s!.cut),")"));

INSTALLPRINTERS@(IsMarkedSphere);

InstallMethod(Draw, "(IMG) for a point in Teichmuller space",
        [IsMarkedSphere],
        function(spider)
    local extrapoints, extraarcs, cid;
    
    RSS.open();
    cid := RSS.newcanvas();
    if IsBound(spider!.map) and ValueOption("nojulia")=fail then
        RSS.putmap(cid,spider!.map);
    fi;
    
    if IsBound(spider!.points) then
        extrapoints := spider!.points;
    else
        extrapoints := [];
    fi;
    if IsBound(spider!.arcs) then
        extraarcs := spider!.arcs;
    else
        extraarcs := [];
    fi;
    
    DRAWPOINTS@(cid,spider!.cut,extrapoints);
    DRAWARCS@(cid,spider!.cut!.e,extraarcs);
end);

InstallOtherMethod(Draw, "(IMG) for a P1 map",
        [IsP1Map],
        function(map)
    Draw(MarkedSphere(SphereMachine(map)));
end);

BindGlobal("CHECKSPIDER@", function(s)
    return CHECKTRIANGULATION@(s!.cut);
end);

BindGlobal("CHECKREC@", function(recur,order,reduce)
    local i, j, a, result, w;
    
    result := [[],[]];
    for i in [1..Length(recur[2][1])] do
        w := One(recur[1][1][1]);
        a := i;
        for j in order do
            w := w*recur[1][j][a];
            a := recur[2][j][a];
        od;
        Add(result[1],reduce(w));
        Add(result[2],a);
    od;
    return result[2]=[1..Length(recur[2][1])] and ForAll(result[1],IsOne);
end);

InstallMethod(VerticesOfMarkedSphere, [IsMarkedSphere],
        function(spider)
    # the vertices a spider lies on
    return List(Filtered(spider!.cut!.v,v->not IsFake(v)),Pos);
end);

InstallMethod(SpanningTreeBoundary, [IsMarkedSphere],
        function(spider)
    # return a list of edges traversed when one surrounds the tree with
    # it on our right. visit vertex n first.
    local i, e, edges, n;

    n := Length(VerticesOfMarkedSphere(spider));
    e := First(spider!.cut!.e,e->spider!.intree[e!.index] and From(e)!.index=n);
    edges := [];
    repeat
        Add(edges,e);
        e := Opposite(e);
        repeat
            e := Next(Opposite(e));
        until spider!.intree[e!.index];
    until IsIdenticalObj(e,edges[1]);
    return edges;
end);

InstallMethod(NewMarkedSphere, "(IMG) for a list of points and a group",
        [IsP1PointCollection,IsSphereGroup],
        function(points,model)
    # constructs a spider with identity marking on <points>
    local n, f, r, g, edges, tree, cost, p, i, e, source, image, ordering;

    n := Length(points);
    while RankOfSphereGroup(model)<>n do
        Error("<model> must have the same rank as number of points");
    od;

    f := FreeGroup(n-1);

    r := rec(model := model,                        # marking group
             cut := DelaunayTriangulation(points,@.maxratio), # triangulation
             group := f,                            # group on spanning tree
             intree := [],                          # if an edge is in the tree
             treeedge := []);                       # for each generator, a preferred edge with that label
    
    # construct a spanning tree
    edges := List(r!.cut!.e,e->[From(e)!.index,To(e)!.index]);
    cost := List(r!.cut!.e,Length);
    tree := MINSPANTREE@(edges,cost);
    r.treecost := Remove(tree);
    tree := List(tree,p->First(Neighbours(r.cut!.v[p[1]]),e->To(e)!.index=p[2]));
    SortParallel(cost{List(tree,e->e!.index)},tree);
    
    # start by a free group on the edges of the tree
    # by convention, if the edge goes north, then the generator, with
    # positive orientation, goes from west to east.
    g := FreeGroup(Length(tree));
    p := PresentationFpGroup(g,0);
    TzOptions(p).protected := Length(tree);
    TzInitGeneratorImages(p);

    for i in r.cut!.e do SetGroupElement(i, One(g)); od;
    r.intree := ListWithIdenticalEntries(Length(edges),false);
    for i in [1..Length(tree)] do
        e := GeneratorsOfGroup(g)[i];
        SetGroupElement(tree[i], e);
        SetGroupElement(Opposite(tree[i]), e^-1);
        r.intree[tree[i]!.index] := true;
        r.intree[Opposite(tree[i])!.index] := true;
    od;
    
    # add relators saying the cycle around a fake vertex is trivial
    for i in r!.cut!.v do
        if IsFake(i) then
            AddRelator(p,Product(List(Reversed(Neighbours(i)),GroupElement)));
        fi;
    od;
    
    # eliminate useless generators, starting by the shortest
    for i in GeneratorsOfPresentation(p) do
        TzEliminate(p,i);
    od;
    for i in r.cut!.e do SetGroupElement(i, One(f)); od;
    for i in [1..Length(tree)] do
        e := MappedWord(TzImagesOldGens(p)[i],GeneratorsOfPresentation(p),GeneratorsOfGroup(f));
        SetGroupElement(tree[i], e);
        SetGroupElement(Opposite(tree[i]), e^-1);
    od;
    
    r.treeedge := List(TzPreImagesNewGens(p),w->tree[TietzeWordAbstractWord(w)[1]]);

    r := Objectify(TYPE_MARKEDSPHERE,r);

    # now set the marking (r.group -> r.model)
    ordering := [];
    source := [];

    for e in SpanningTreeBoundary(r) do
        if not IsFake(From(e)) then
            i := From(e)!.index;
            if not IsBound(source[i]) then
                source[i] := One(f);
                Add(ordering,i);
            fi;
            source[i] := LeftQuotient(GroupElement(e), source[i]);
        fi;
    od;

    f := FreeGroupOfFpGroup(model);
    g := GeneratorsOfGroup(f);
    image := GroupHomomorphismByImages(f,f,g{ordering},List([1..n],i->g[ordering[i]]^Inverse(Product(g{Difference([1..ordering[i]],ordering{[1..i]})}))));

    ordering := OrderingOfSphereGroup(model);
    image := image*GroupHomomorphismByImages(f,f,List([1..n],i->g[ordering[i]]^Inverse(Product(g{Difference([1..ordering[i]],ordering{[1..i]})}))),g{ordering});

    r!.marking := GroupHomomorphismByImages(r!.group,model,source,List(g,x->ElementOfSphereGroup(FamilyObj(One(model)),x^image)));
    
    return r;
end);
InstallMethod(NewMarkedSphere, "(IMG) for a list of points",
        [IsP1PointCollection],
        points->NewMarkedSphere(points,SphereGroup(Length(points))));

reorder := function(sigma)
    local f, g;
    f := FreeGroup(Length(sigma));
    g := GeneratorsOfGroup(f);
    return GroupHomomorphismByImages(f,f,g{sigma},List([1..Length(sigma)],i->g[sigma[i]]^(Product(Difference([1..sigma[i]],sigma{[1..i]}),j->g[j])^-1)));
end;

InstallMethod(WiggledMarkedSphere, "(IMG) for a marked sphere and a moebius xfo/new vx positions",
        [IsMarkedSphere,IsObject],
        function(spider,movement)
    # movement is either a list of new positions for the vertices
    # (in which case the vertices, edges etc. should be wiggled to
    # their new positions), or a Möbius transformation, or true.
    local r;
    r := rec(model := spider!.model,
             cut := WiggledTriangulation(spider!.cut,movement),
             group := spider!.group,
             marking := spider!.marking,
             treecost := spider!.treecost,
             intree := spider!.intree);
    r.treeedge := r.cut!.e{List(spider!.treeedge,e->e!.index)};

    return Objectify(TYPE_MARKEDSPHERE,r);
end);

InstallMethod(ShallowCopy, "(IMG) for a marked sphere",
        [IsMarkedSphere],
        spider->WiggledMarkedSphere(spider,fail));
##############################################################################

##############################################################################
##
#M  P1 map to sphere machine
##
BindGlobal("ESSDISJOINT@", function(ratmap,p0,p1,domain)
    # return true if ratmap(P1Path(p0,p1)) is essentially disjoint
    # from domain's boundary
    local e, tu, delta, d;
    
    d := P1Distance(p0,p1);
    if d<@.p1eps then return true; fi;
    
    delta := P1Path(p0,p1);
    for e in Neighbours(domain) do
MARKTIME@(1);
        tu := P1INTERSECT@(Map(e),ratmap,delta);
MARKTIME@(2);
        for tu in tu do
            if d*(@.ro/2-AbsoluteValue(@.ro/2-tu[2])) > @.p1eps then
                return false;
            fi;
        od;
    od;
    return true;
end);

BindGlobal("LIFTARC@", function(spider,ratmap,from,to,gamma,downcell)
    # <gamma> is an arc in the range, contained in face <downcell>, which we
    # want to lift through <ratmap>.
    # <from> and <to> are described in LIFTEDGE@, as is the return value
    local curtime, nexttime, curface, curpos, curelt, curbdry, curedge,
          lift, lifts, xings, candidates,
          toface, e, f, i, c,
          choosebysubdivision, getcandidates, getbdry, getxings;

    choosebysubdivision := function(p0,t0,t1,candidates,upbdry)
    # candidates is a list of records containing in particular a field pos.
    # returns the one such that gamma[t0,t1] (which stays in downcell)
    # lifts to a path from p0 to candidate.pos and staying in upcell.
    # works by subdividing the time interval.
        local c, d, l, i, p, subdiv;
    Info(InfoIMG,3,"choosebysubdivision ",[p0,t0,t1,candidates,upbdry]);
        subdiv := [rec(t := t0, pos := p0), rec(t := t1)];
        i := 2;
        while i <= Length(subdiv) do
            if not IsBound(subdiv[i].lifts) then
                subdiv[i].lifts := [];
                for p in P1PreImages(ratmap,P1Image(gamma,P1Point(subdiv[i].t))) do
                    if ForAll(upbdry,e->SphereP1Y(P1Image(InverseP1Map(Map(e)),p))>-@.p1eps) then
                        Add(subdiv[i].lifts,p);
                    fi;
                od;
            fi;
            p := subdiv[i].lifts;
            l := Length(p);
            if l=0 then
                Error("LIFTARC: I can't analytically continue, no candidate");
            elif l>1 then # many choices, find good one(s)
                p := Filtered(p,p->ESSDISJOINT@(ratmap,subdiv[i-1].pos,p,downcell));
                l := Length(p);
                if l>1 then
                    Error("LIFTARC: I can't analytically continue, too many candidates");
                elif l=0 then  # subdivide
                    Add(subdiv,rec(t := (subdiv[i-1].t+subdiv[i].t)/2),i);
                    continue;
                fi;
            fi;
            p := p[1];
            subdiv[i].pos := p;
            i := i+1;
        od;
        
        # now find closest candidate to point p.
        l := 2*@.pi; # bigger than maximal distance between P1 points
        for i in candidates do
            d := P1Distance(i.pos,p);
            if d<l then l := d; c := i; fi;
        od;
        return c;
    end;
    
    # compute a sequence of edges surrounding our current position
    getbdry := function()
        local e;
        
        if curedge<>fail then # we're parallel to an edge, i.e.
            curbdry := [];
            for e in Neighbours(curface) do
                if not IsIdenticalObj(e,Opposite(curedge)) then
                    Add(curbdry,e);
                fi;
            od;
            for e in Neighbours(Left(curedge)) do
                if not IsIdenticalObj(e,curedge) then
                    Add(curbdry,e);
                fi;
            od;
        else
            curbdry := Neighbours(curface);
        fi;
    end;
    
    # compute edge intersections on neighbours of current cell
    getxings := function()
        local e, l, r, i, tu;
        
        for e in Neighbours(curface) do
            if not IsBound(xings[e!.index]) then
                # get list of [t,u,d,p,q] such that gamma(t)=delta(u)=p,
                # e.map(u)=q, d=Im(gamma^-1*delta)'(u)
MARKTIME@(1);
                tu := P1INTERSECT@(gamma,ratmap,Map(e));
MARKTIME@(2);
                # in increasing order along the edge
                Sort(tu,function(x,y) return x[2]<y[2]; end);
                l := []; r := []; i := 1;
                for tu in tu do
                    Add(l, rec(t := tu[1], u := tu[2], d := tu[3], gammapos := tu[4], pos := tu[5], e := e));
                    Add(r, rec(t := tu[1], u := (@.ro-tu[2])/(@.ro+tu[2]), d := -tu[3], gammapos := tu[4], pos := tu[5], e := Opposite(e)));
                    l[i].reverse := r[i];
                    r[i].reverse := l[i];
                    i := i+1;
                od;
                xings[e!.index] := l;
                xings[Opposite(e)!.index] := Reversed(r);
            fi;
        od;
    end;
    
    # out of the xings, find candidates:
    # - all endpoints of paths (in list "to")
    # - among edge crossings, only those at time >= curtime
    # - if we're parallel to an edge, all on that edge
    # - for the other ("boundary") crossings, only those pointing
    #   outward, and separated by (algebraically) >= #to and <= #from
    #   on the current face(s).
    getcandidates := function()
        local c, e;
        
        if IsBound(toface[curface!.index]) then
            candidates := ShallowCopy(toface[curface!.index]);
            nexttime := @.ro;
        else
            candidates := [];
            nexttime := 2*@.ro; # out of the way
        fi;
        
        if curedge<>fail then # we're parallel to an edge, i.e.
            # inside a lozenge. take its boundary.
            for c in xings[curedge!.index] do
                if c.t > curtime then
                    if c.t < nexttime then nexttime := c.t; fi;
                    Add(candidates,c);
                fi;
            od;
        fi;
    
        for e in curbdry do
            for c in xings[e!.index] do
                if c.t > curtime and c.d >= 0 then
                    if c.t < nexttime then nexttime := c.t; fi;
                    Add(candidates,c);
                fi;
            od;
        od;
    end;
    
    # the results will go there
    lifts := [];

    # xings[edge.index] are the (left-to-right) crossings of gamma with f(edge)
    xings := [];
    
    # toface is a list indexed by faces, containing ends of lifts
    toface := [];
    for f in to do
        if not IsBound(toface[f.face!.index]) then toface[f.face!.index] := []; fi;
        Add(toface[f.face!.index],f);
    od;

    for lift in from do
        curface := lift.face;
        curedge := fail; # will keep track of edge to which we're parallel
        curpos := lift.pos;
        curelt := lift.elt;
        curtime := -@.p1eps; # just in case we're on an edge or vertex
        Info(InfoIMG,4,"Ready to lift ",lift);
        
        repeat
            getbdry();
            getxings();
            getcandidates();
            Info(InfoIMG,4,"Lifting at time ",curtime,": pos=",curpos," in ",[curface,curedge]);
            
            c := First(candidates,c->ESSDISJOINT@(ratmap,curpos,c.pos,downcell));
            if candidates=[] then # our last chance is that we're at an edge or vertex
                Info(InfoIMG,4,"Searching for endgame move to target cell");
                for f in toface do
                    for f in f do
                        if INID@(curface,ClosestFaces(f.cell)) and ESSDISJOINT@(ratmap,curpos,f.pos,downcell) then
                            curelt := curelt * Product(List(EdgePath(spider!.cut,curface,f.face),GroupElement));
                            curface := f.face;
                            c := f; # signal done
                            break;
                        fi;
                    od;
                    if c<>fail then break; fi; # "goto done"
                od;
                if c=fail then
                    Error("I'm stuck trying to advance, there are no candidates.");
                fi;
            elif c=fail then                
                c := [];
                for i in candidates do
                    if (nexttime=@.ro and not IsBound(i.t)) or i.t=nexttime then
                        Add(c,i);
                    fi;
                od;
                c := choosebysubdivision(curpos,curtime,nexttime,c,curbdry);
                if not IsBound(c.e) and not IsBound(c.face) then # got a new point in curface
                    curpos := c.pos;
                    curtime := nexttime;
                    continue;
                fi;
            fi;
            
            while not IsRecord(c) do # something went awfully wrong
                Error("I can't move any further, but haven't reached any endpoint. Repent.");
            od;
            
            if IsBound(c.face) then # "to" cell: done!
                lift := ShallowCopy(Remove(toface[curface!.index],POSITIONID@(toface[curface!.index],c)));
                lift.elt := curelt;
                break;
            fi;
            
            # if we're parallel to an edge, maybe move back to the previous cell
            if not IsIdenticalObj(Left(c.e),curface) then
                Error("This code is probably not necessary. Please contact laurent.bartholdi@gmail.com if it gets triggered");
                curelt := curelt / GroupElement(c.e);
            fi;
            
            curface := Right(c.e);
            if c.d=0 then curedge := c.e; else curedge := fail; fi;
            curelt := curelt * GroupElement(c.e);
            curtime := c.t;
            curpos := c.pos;
            # if at a vertex, allow time to go back a little, in case the
            # edges don't really match
            if c.u < @.p1eps or c.u > @.ro-@.p1eps then
                curtime := curtime - 10*@.reps;
            fi;
            c.t := -@.ro; # mark xing, and its reverse, as unusable
            c.reverse.t := -@.ro;
        until false;
        Add(lifts,lift);
    od;
    return lifts;
end);

BindGlobal("LIFTEDGE@", function(spider,ratmap,from,to,edge)
    # lifts the arc perpendicular to <edge> through <ratmap>.
    # <from> is a list of rec(pos := <p1point>, cell := <face>,
    #     elt := <gpelement>), such that the <p1point> are the
    #     preimages of edge.left.pos.
    # <to> is a lift of rec(pos := <p1point>, cell := <face>), one per
    #     preimage of e.right.pos.
    # returns list <lifts> of length Degree(ratmap), where
    #     <lifts>[i] is a rec(pos := <p1point>, cell := <face>,
    #     elt := <gpelement>); this is a reordering of <to>,
    #     such that from[i] continues to to[i], and
    #     to[i].elt = from[i].elt * (product of edges crossed along the lift)
    local mid, r, p;
    
    Info(InfoIMG,3,"Lifting edge ",edge);
    
    mid := [];
    for p in P1PreImages(ratmap,Pos(edge)) do
        r := rec(pos := p, cell := LocateInTriangulation(spider!.cut,p));
        r.face := ClosestFace(r.cell);
        Add(mid,r);
    od;
    
    mid := LIFTARC@(spider,ratmap,from,mid,P1Path(Pos(Left(edge)),Pos(edge)),Left(edge));
    return LIFTARC@(spider,ratmap,mid,to,P1Path(Pos(edge),Pos(Right(edge))),Right(edge));
end);

BindGlobal("NORMALIZEADDINGMACHINE@", function(srcmodel,liftmodel,trans,out,srcadder,liftadder)
    # conjugate the recursion so that element adder, which is checked to
    # be an adding machine, becomes of the form (t,...,1)s,
    # where s is the cycle i|->i-1 mod d.
    # adder is the position of the adding element.
    # model is the ambient fundamental group.
    local cycle, deg, perm, x, i, j, basis;
    
    deg := Length(trans[1]);
    cycle := Cycles(PermList(out[srcadder]),[1..deg]);
    while Length(cycle)<>1 or not IsConjugate(liftmodel,Product(trans[srcadder]{cycle[1]}),liftmodel.(liftadder)) do
        Error("Element ",srcmodel.(srcadder)," is not an adding element");
    od;
    
    perm := PermList(Concatenation([deg],[1..deg-1]));
    perm := RepresentativeAction(SymmetricGroup(deg),PermList(out[srcadder]),perm);
    REORDERREC@([trans,out],perm);

    basis := [];
    x := One(liftmodel);
    for i in [deg,deg-1..1] do
        basis[i] := x;
        x := x*trans[srcadder][i];
    od;
    basis := RepresentativeAction(liftmodel,liftmodel.(liftadder),x)*basis;
    for i in [1..Length(trans)] do
        for j in [1..deg] do
            trans[i][j] := basis[j]*trans[i][j]/basis[out[i][j]];
        od;
    od;
end);

InstallMethod(SphereMachineOfBranchedCovering, "(IMG) for two marked spheres, a rational map, and a boolean",
        [IsMarkedSphere,IsMarkedSphere,IsP1Map,IsBool],
        function(src,lift,ratmap,poly)
    # lifts all dual arcs in <src> through <ratmap>; rounds their endpoints
    # to faces of <lift>; and rewrites the generators of <src> as words
    # in <lift>'s group. <base> is a preferred starting face of <src>.
    # returns [face,edge] where:
    # face is a list of length Degree(ratmap), and contains lifts of faces,
    # indexed by the faces of <src>
    # face[i][j] is rec(pos, liftface, liftgpelt)
    local face, f, e, r, i, j, todo, lifts, perm, state, p, s, base, idle, machine;
    
    # first lift all face centres, and choose a face containing the lift
    face := [];
    for f in src!.cut!.f do
        lifts := [];
        for p in P1PreImages(ratmap,Pos(f)) do
            r := rec(pos := p, cell := LocateInTriangulation(lift!.cut,p));
            r.face := ClosestFace(r.cell);
            Add(lifts,r);
        od;
        Add(face,lifts);
    od;
    
    # and choose a base point
    if poly then
        base := Left(Neighbour(src!.cut!.v[Length(GeneratorsOfGroup(src!.group))])); # some face touching infinity
    else
        base := src!.cut!.f[1];
    fi;
    for f in face[base!.index] do
        f.elt := One(lift!.group);
    od;
    
    # lift edges in the dual tree. If src!.cut!.f[i] lifts to points
    # in lift!.cut!.f[j_1]...lift!.cut!.f[j_d], then face[i][k], for k=1..d,
    # is a record (cell=j_k, elt=the word obtained by lifting the geodesic
    # from the basepoint to j_i, pos=exact position of the endpoint).
    todo := NewFIFO([base]);
    for f in todo do
        for e in Neighbours(f) do
            # face[index] is a list of rec(pos := <position in P1>,
            #    cell := <cell in lift>, and maybe elt := <group element>).
            # if elt is not assigned, we haven't lifted the edge yet
            r := Right(e)!.index;
            if not src!.intree[e!.index] and not IsBound(face[r][1].elt) then
                face[r] := LIFTEDGE@(lift,ratmap,face[f!.index],face[r],e);
                Add(todo,Right(e));
            fi;
        od;
    od;

    # then lift edges cutting the tree; store group elements and permutations
    # in [perm,state]
    perm := [];
    state := [];
    for e in src!.treeedge do
        r := Right(e)!.index;
        lifts := LIFTEDGE@(lift,ratmap,face[Left(e)!.index],face[r],e);
        p := [];
        s := [];
        for i in [1..Length(lifts)] do
            j := PositionProperty(face[r],f->IsIdenticalObj(lifts[i].pos,f.pos));
            Add(p,j);
            Add(s,lifts[i].elt/face[r][j].elt);
        od;
        Add(perm,p);
        Add(state,s);
    od;

    # lift points, if present -- this should give an approximation of the measure of maximal entropy
    if IsBound(src!.points) then
        lift!.points := [];
        for i in src!.points do
            Add(lift!.points, Random(P1PreImages(ratmap,i)));
        od;
    fi;

    perm := COMPOSERECURSION@(state,perm,src!.marking,lift!.marking);
    state := perm[1];
    perm := perm[2];

    if poly then
        NORMALIZEADDINGMACHINE@(src!.model,lift!.model,state,perm,Length(state),Length(GeneratorsOfGroup(lift!.model)));
    fi;

    machine := FRMachine(src!.model,state,perm);
    IsSphereMachine(machine); # set filter

    if poly then
        SetAddingElement(machine,FRElement(machine,src!.model.(Length(state))));
    fi;

    return machine;
end);
InstallMethod(SphereMachineOfBranchedCovering, "(IMG) for two marked spheres, and a rational map",
        [IsMarkedSphere,IsMarkedSphere,IsP1Map],
        function(src,lift,ratmap)
    return SphereMachineOfBranchedCovering(src,lift,ratmap,IsPolynomial(ratmap));
end);

InstallMethod(SphereMachineAndSphereOfBranchedCovering, "(IMG) for a marked spheres, a rational map and a bool",
        [IsMarkedSphere,IsP1Map,IsBool],
        function(src,ratmap,poly)
    local lift, v, w;
    lift := ShallowCopy(CriticalPointsOfP1Map(ratmap));
    for v in VerticesOfMarkedSphere(src) do
        for w in P1PreImages(ratmap,v) do
            if ForAll(lift,x->P1Distance(w,x)>@.p1eps) then
                Add(lift,w);
            fi;
        od;
    od;
    lift := NewMarkedSphere(lift);
    return [SphereMachineOfBranchedCovering(src,lift,ratmap,poly),src];
end);
InstallMethod(SphereMachineAndSphereOfBranchedCovering, "(IMG) for a marked sphere, and a rational map",
        [IsMarkedSphere,IsP1Map],
        function(src,ratmap)
    return SphereMachineAndSphereOfBranchedCovering(src,ratmap,IsPolynomial(ratmap));
end);

InstallMethod(MonodromyOfP1Map, "(IMG) for a marked sphere and a rational map",
        [IsMarkedSphere,IsP1Map],
        function(down,ratmap)
    local machine;
    machine := SphereMachineAndSphereOfBranchedCovering(down,ratmap)[1];
    return GroupHomomorphismByImages(down!.model,SymmetricGroup(Length(machine!.trans[1])),GeneratorsOfGroup(down!.model),List(machine!.trans,PermList));
end);
InstallMethod(MonodromyOfP1Map, "(IMG) for a list of points and a rational map",
        [IsP1PointCollection,IsP1Map],
        function(downpts,ratmap)
    local down, up;
    down := NewMarkedSphere(downpts);
    up := NewMarkedSphere(List(CollectedP1Points(CriticalPointsOfP1Map(ratmap)),x->x[1]));
    return List(SphereMachineOfBranchedCovering(down,up,ratmap)!.output,PermList);
end);
InstallMethod(MonodromyOfP1Map, "(IMG) for a rational map",
        [IsP1Map],
        function(ratmap)
    local cp, cv, down, up;
    cp := List(CollectedP1Points(CriticalPointsOfP1Map(ratmap)),x->x[1]);
    cv := List(CollectedP1Points(List(cp,v->P1Image(ratmap,v))),x->x[1]);
    up := NewMarkedSphere(cp);
    down := NewMarkedSphere(cv);
    return List(SphereMachineOfBranchedCovering(down,up,ratmap)!.output,PermList);
end);

InstallMethod(SphereMachine, "(IMG) for a rational function",
        [IsP1Map],
        function(f)
    local i, poly, pcdata, pcp, spider, machine, x;
    
    pcdata := POSTCRITICALPOINTS@(f);
    poly := pcdata[1];
    pcp := pcdata[3];
    Info(InfoIMG,2,"Post-critical points at ",pcdata[3]);

    spider := NewMarkedSphere(pcp);
    spider!.map := f;
    spider!.cycle := PCDATAATTRACTINGCYCLES@(pcdata);

    machine := SphereMachineOfBranchedCovering(spider,spider,f,poly);
    SetMarkedSphere(machine,spider);

    return machine;
end);

InstallMethod(DistanceMarkedSpheres, "(IMG) for two marked spheres and a bool",
        # this is some sort of Teichmüller distance. The marked spheres
        # are assumed to be normalized, e.g. by fixing 3 points as 0,1,infty.
        [IsMarkedSphere,IsMarkedSphere,IsBool],
        function(spiderA,spiderB,fast)
    local model, points, perm, dist, recur, endo, nf, g;

    model := spiderA!.model;
    while spiderB!.model<>model do
        Error("DistanceMarkedSpheres: the spheres don't have the same model group");
    od;

    # move points of spiderB to their spiderA matches
    points := VerticesOfMarkedSphere(spiderA);    
    spiderB := WiggledMarkedSphere(spiderB,points);
    dist := spiderB!.cut!.wiggled;
    
    if dist>@.fast then # crude estimate
        return @.ro*Sum(GeneratorsOfGroup(spiderA!.group),x->Length(PreImagesRepresentative(spiderB!.marking,x^spiderA!.marking)))/Length(points);
    fi;
    
    if fast then # we just wiggled the points, the combinatorics didn't change
        return dist/Length(points);
    fi;
    
    recur := SphereMachineOfBranchedCovering(spiderA,spiderB,P1z,false);

    endo := List(recur!.transitions,x->x[1]);
    REDUCEINNER@(endo,GeneratorsOfMonoid(model));
    
    for g in endo do
        dist := dist + (Length(g)-1); # if each image is a gen, then endo=1
    od;
    return dist/Length(points);
end);
InstallMethod(DistanceMarkedSpheres, "(IMG) for two marked spheres",
        [IsMarkedSphere,IsMarkedSphere],
        function(spiderA,spiderB)
    return DistanceMarkedSpheres(spiderA,spiderB,false);
end);

##############################################################################

#E markedsphere.gi . . . . . . . . . . . . . . . . . . . . . . . . . ends here
