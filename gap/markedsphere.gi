#! probably: mess-up in IMG relations, causing the ordering to be sometimes 1,2,3, sometimes 3,2,1

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

BindGlobal("STRINGCOMPLEX@",
    z->CONCAT@FR(RealPart(z)," ",ImaginaryPart(z)));

InstallMethod(Draw, "(IMG) for a point in Teichmuller space",
        [IsMarkedSphere],
        function(spider)
    local a, i, j, k, s, f, t, points, arcs;
    s := ""; f := OUTPUTTEXTSTRING@FR(s);
    
    if ValueOption("upper")<>fail then
        PrintTo(f,"UPPER\n");
    fi;
    if ValueOption("lower")<>fail then
        PrintTo(f,"LOWER\n");
    fi;
    if IsBound(spider!.map) and ValueOption("julia")<>fail then
        t := DegreeOfP1Map(spider!.map);
        PrintTo(f,"FUNCTION");
        a := List(CoefficientsOfP1Map(spider!.map),ShallowCopy);
        for i in [1..t+1] do PrintTo(f," ",STRINGCOMPLEX@(a[1][i])); od;
        for i in [1..t+1] do PrintTo(f," ",STRINGCOMPLEX@(a[2][i])); od;
        PrintTo(f,"\nCYCLES");
        if IsBound(spider!.cycle) then
            for i in spider!.cycle do
                if i[1]=P1infinity then
                    PrintTo(f," Infinity any");
                else
                    PrintTo(f," ",STRINGCOMPLEX@(P1Coordinate(i[1])));
                fi;
                PrintTo(f," ",i[2]," ",i[3]);
            od;
        fi;
        t := ValueOption("julia");
        if IsList(t) then # size, maxiter
            i := t[1]; j := t[1];
        elif IsPosInt(t) then
            i := t; j := 100;
        else
            i := 500; j := 100;
        fi;
        PrintTo(f,"\nIMAGE ",i," ",j,"\n");
    fi;
    
    if IsBound(spider!.points) then
        points := spider!.points;
    else
        points := [];
    fi;
    if IsBound(spider!.arcs) then
        arcs := spider!.arcs;
    else
        arcs := [];
    fi;
    
    t := spider!.cut;
    PRINTPOINTS@(f, t, points);
    PRINTARCS@(f, t!.e, arcs, Float(101/100));
    
    Info(InfoIMG,3,"calling javaplot with:\n",s);
    JAVAPLOT@(InputTextString(s));
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

InstallMethod(Vertices, [IsMarkedSphere],
        function(spider)
    # the vertices a spider lies on
    return List(Filtered(spider!.cut!.v,v->not IsFake(v)),Pos);
end);

InstallMethod(SpanningTreeBoundary, [IsMarkedSphere],
        function(spider)
    # return a list of edges traversed when one surrounds the tree with
    # it on our right. visit vertex n first.
    local i, e, edges, n;

    n := Length(Vertices(spider));
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
            source[i] := source[i] / GroupElement(e);
        fi;
    od;

    f := FreeGroupOfFpGroup(model);
    g := GeneratorsOfGroup(f);
    ordering := Reversed(ordering);
    image := GroupHomomorphismByImages(f,f,g{ordering},List([1..n],i->g[ordering[i]]^Inverse(Product(g{Difference([1..ordering[i]],ordering{[1..i]})}))));
    ordering := FamilyObj(One(model))!.ordering;
    image := image*GroupHomomorphismByImages(f,f,List([1..n],i->g[ordering[i]]^Inverse(Product(g{Difference([1..ordering[i]],ordering{[1..i]})}))),g{ordering});

    r!.marking := NORMALIZEHOMOMORPHISM@(GroupHomomorphismByImages(r!.group,model,source,List(g,x->ElementOfSphereGroup(FamilyObj(One(model)),x^image))));
    
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
#M  Post-critical machine
##
BindGlobal("POSTCRITICALPOINTS@", function(f)
    # return [poly,[critical points],[post-critical points],[transitions]]
    # where poly=true/false says if there is a fixed point of maximal degree;
    # it is then the last element of <post-critical points>
    # critical points is a list of [point in P1,degree]
    # post-critical points are points in P1
    # post-critical graph is a list of [i,j,n] meaning pcp[i] maps to pcp[j]
    # with local degree n>=1; or, if i<0, then cp[-i] maps to pcp[j].

    local c, i, j, cp, pcp, n, deg, newdeg, poly, polypos,
          transitions, src, dst;

    deg := DegreeOfP1Map(f);
    cp := CollectedP1Points(CriticalPointsOfP1Map(f));
    cp := List(cp,x->[x[1],x[2]+1]);
    poly := First([1..Length(cp)],i->cp[i][2]=deg and P1Distance(P1Image(f,cp[i][1]),cp[i][1])<@.p1eps);
    
    pcp := [];
    transitions := [];
    n := 0;
    for i in [1..Length(cp)] do
        c := cp[i][1];
        src := -i;
        deg := cp[i][2];
        repeat
            c := P1Image(f,c);
            j := PositionProperty(cp,x->P1Distance(c,x[1])<@.p1eps);
            if j<>fail then
                c := cp[j][1];
                newdeg := cp[j][2];
            else
                newdeg := 1;
            fi;
            dst := PositionProperty(pcp,d->P1Distance(c,d)<@.p1eps);
            if dst=fail then
                if j=fail then
                    Add(pcp,c);
                else
                    Add(pcp,cp[j][1]);
                fi;
                if RemInt(Length(pcp),100)=0 then
                    Info(InfoIMG,2,"Post-critical set contains at least ",Length(pcp)," points");
                fi;
                dst := Length(pcp);
                Add(transitions,[src,dst,deg]);
                n := n+1;
                if IsInt(poly) and IsIdenticalObj(pcp[n],cp[poly][1]) then
                    polypos := n;
                    poly := true;
                fi;
            else
                Add(transitions,[src,dst,deg]);
                break;
            fi;
            deg := newdeg;
            src := dst;
        until false;
    od;

    if poly=fail then
        poly := false;
    else
        Add(pcp,Remove(pcp,polypos)); # force infinity to be at end
        for c in transitions do
            for i in [1..2] do
                if c[i]=polypos then
                    c[i] := n;
                elif c[i]>polypos then
                    c[i] := c[i]-1;
                fi;
            od;
        od;
    fi;

    return [poly,cp,pcp,transitions];
end);

InstallGlobalFunction(PostCriticalMachine, function(f)
    local trans, i, pcp, machine;
    trans := [];
    pcp := POSTCRITICALPOINTS@(AsP1Map(f));
    for i in pcp[4] do
        if i[1]>0 then trans[i[1]] := [i[2]]; fi;
    od;
    machine := MealyMachineNC(FRMFamily([1]),trans,List(trans,x->[1]));
    SetCorrespondence(machine,pcp[3]);
    return machine;
end);

BindGlobal("ATTRACTINGCYCLES@", function(pcdata)
    local cycle, period, len, next, i, j, jj, periodic, critical;
    
    cycle := [];
    next := [];
    period := [];
    for i in [1..Length(pcdata[3])] do
        critical := false; periodic := false;
        j := i; jj := i;
        repeat
            jj := First(pcdata[4],x->x[1]=jj)[2];
            jj := First(pcdata[4],x->x[1]=jj)[2];
            j := First(pcdata[4],x->x[1]=j)[2];
        until j=jj;
        len := 0;
        repeat
            len := len+1;
            periodic := periodic or i=j;
            j := First(pcdata[4],x->x[1]=j);
            critical := critical or j[3]>1;
            j := j[2];
        until j=jj;
        if critical and periodic then
            Add(cycle,pcdata[3][i]);
            Add(next,i);
            Add(period,len);
        fi;
    od;
    next := List(next,i->Position(next,First(pcdata[4],x->x[1]=i)[2])-1);
    return TransposedMat([cycle,next,period]);
end);
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
    
        subdiv := [rec(t := t0, pos := p0), rec(t := t1)];
        i := 2;
        while i <= Length(subdiv) do
            if not IsBound(subdiv[i].lifts) then
                subdiv[i].lifts := [];
                for p in P1PreImages(ratmap,P1Image(gamma,P1Point(subdiv[i].t))) do
                    if ForAll(upbdry,e->SphereP1Y(P1Image(InverseP1Map(Map(e)),p))>-@.reps) then
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
    
    # toface is a list indexed by faces, containing starts and ends of lifts
    toface := [];
    for f in to do
        if not IsBound(toface[f.cell!.index]) then toface[f.cell!.index] := []; fi;
        Add(toface[f.cell!.index],f);
    od;

    for lift in from do
        curface := lift.cell;
        curedge := fail; # will keep track of edge to which we're parallel
        curpos := lift.pos;
        curelt := lift.elt;
        curtime := -@.p1eps; # just in case we're on an edge or vertex

        repeat
            getbdry();
            getxings();
            getcandidates();
            
            c := First(candidates,c->ESSDISJOINT@(ratmap,curpos,c.pos,downcell));
            if candidates=[] then # our last chance is a "to" on the other side
                lift := fail;
                for e in Neighbours(curface) do if IsBound(toface[Right(e)!.index]) then
                    for c in toface[Right(e)!.index] do
                        if ESSDISJOINT@(ratmap,curpos,c.pos,downcell) then
                            curelt := curelt * GroupElement(e);
                            curface := Right(e);
                            lift := c;
                            break;
                        fi;
                    od;
                    if lift<>fail then break; fi;
                fi; od;
            elif c=fail then                
                c := [];
                for i in candidates do
                    if (nexttime=@.ro and not IsBound(i.t)) or i.t=nexttime then
                        Add(c,i);
                    fi;
                od;
                c := choosebysubdivision(curpos,curtime,nexttime,c,curbdry);
                if not IsBound(c.e) and not IsBound(c.cell) then # got a new point in curface
                    curpos := c.pos;
                    curtime := nexttime;
                    continue;
                fi;
            fi;
            
            if IsBound(c.cell) then # "to" cell: done!
                lift := ShallowCopy(Remove(toface[curface!.index],POSITIONID@(toface[curface!.index],c)));
                lift.elt := curelt;
                break;
            fi;
            
            # if we're parallel to an edge, maybe move back to the previous cell
            if not IsIdenticalObj(Left(c.e),curface) then
                Error("This code is probably not necessary @@");
                curelt := curelt / GroupElement(c.e);
            fi;
            
            curface := Right(c.e);
            if c.d=0 then curedge := c.e; else curedge := fail; fi;
            curelt := curelt * GroupElement(c.e);
            curtime := c.t;
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
    local mid;
    
    Info(InfoIMG,3,"Lifting edge ",edge);

    mid := List(P1PreImages(ratmap,Pos(edge)),y->rec(pos := y, cell := LocateInTriangulation(spider!.cut,y)));
    
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
        Error("Element #",liftadder," is not an adding element");
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
    face := List(src!.cut!.f,x->List(P1PreImages(ratmap,Pos(x)),y->rec(pos := y, cell := LocateInTriangulation(lift!.cut,y))));
    
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

    perm := COMPOSERECURSION@(state,perm,InverseGeneralMapping(src!.marking),InverseGeneralMapping(lift!.marking));
    state := perm[1];
    perm := perm[2];

    if poly then
        NORMALIZEADDINGMACHINE@(src!.model,lift!.model,state,perm,Length(state),Length(GeneratorsOfGroup(lift!.model)));
    fi;

    machine := FRMachine(src!.model,state,perm);

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
    for v in Vertices(src) do
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
    spider!.cycle := ATTRACTINGCYCLES@(pcdata);

    machine := SphereMachineOfBranchedCovering(spider,spider,f,poly);
    SetMarkedSphere(machine,spider);
    SetP1Map(machine,f);

    return machine;
end);

InstallMethod(DistanceMarkedSpheres, "(IMG) for two marked spheres and a bool",
        [IsMarkedSphere,IsMarkedSphere,IsBool],
        function(spiderA,spiderB,fast)
    local model, points, perm, dist, recur, endo, nf, g;

    model := spiderA!.model;

    # try to match feet of spiderA and spiderB
    points := Vertices(spiderA);
    perm := Vertices(spiderB);
    
    perm := MatchP1Points(perm,List(perm,x->points));
    if perm=fail or Set(perm)<>[1..Length(points)] then # no match, find something coarse
        return @.ro*Sum(GeneratorsOfGroup(spiderA!.group),x->Length(PreImagesRepresentative(spiderA!.marking,x)^spiderB!.marking))/Length(points);
    fi;
    
    
    # move points of spiderB to their spiderA matches
    spiderB := WiggledMarkedSphere(spiderB,points{perm});
    dist := spiderB!.cut!.wiggled;
    
    if fast then # we just wiggled the points, the combinatorics didn't change
        return dist/Length(points);
    fi;
    
    recur := SphereMachineOfBranchedCovering(spiderA,spiderB,P1z,false);

    endo := List(recur!.transitions,x->x[1]);
    REDUCEINNER@(endo,GeneratorsOfMonoid(spiderB!.group),x->x);
    
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
