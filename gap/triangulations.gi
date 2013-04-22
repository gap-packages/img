#############################################################################
##
#W triangulations.gi                                        Laurent Bartholdi
##
#Y Copyright (C) 2011-2013, Laurent Bartholdi
##
#############################################################################
##
##  Triangulations of spheres
##
#############################################################################

################################################################

InstallMethod(\., [IsTriangulationObject, IsInt], function(x,n)
    n := NameRNam(n);
    Info(InfoIMG,1,"Access to triangulation's ",n);
    Error("Interrupt -- giving you a chance to clean up the code");
    if n="n" then return Neighbours(x);
    elif n="pos" then return Pos(x);
    elif n="index" then return x!.index;
    elif n="from" then return From(x);
    elif n="to" then return To(x);        
    elif n="to" then return To(x);        
    elif n="left" then return Left(x);        
    elif n="right" then return Right(x);        
    elif n="reverse" then return Opposite(x);        
    elif n="map" then return Map(x);
    elif n="len" then return Length(x);
    elif n="radius" then return Radius(x);
    else Error("Don't know how to lookup ",n);
    fi;
end);

InstallMethod(ViewString, [IsTriangulationVertex],
        v->CONCAT@FR("<vertex ",v!.index,List(Neighbours(v),e->e!.index),">"));
InstallMethod(ViewString, [IsTriangulationEdge],
        e->CONCAT@FR("<edge ",e!.index,List([From(e),To(e)],v->v!.index),">"));
InstallMethod(ViewString, [IsTriangulationFace],
        f->CONCAT@FR("<face ",f!.index,List(Neighbours(f),e->e!.index),">"));
INSTALLPRINTERS@(IsTriangulationObject);

Perform([IsTriangulationVertex,IsTriangulationEdge,IsTriangulationFace], function(filter)
    InstallMethod(\=, [filter,filter], function(a,b) return a!.index=b!.index; end);
    InstallMethod(\<, [filter,filter], function(a,b) return a!.index<b!.index; end);
end);
InstallMethod(\=, [IsTriangulationVertex,IsTriangulationEdge], ReturnFalse);
InstallMethod(\<, [IsTriangulationVertex,IsTriangulationEdge], ReturnTrue);
InstallMethod(\=, [IsTriangulationVertex,IsTriangulationFace], ReturnFalse);
InstallMethod(\<, [IsTriangulationVertex,IsTriangulationFace], ReturnTrue);
InstallMethod(\=, [IsTriangulationEdge,IsTriangulationFace], ReturnFalse);
InstallMethod(\<, [IsTriangulationEdge,IsTriangulationFace], ReturnTrue);


InstallMethod(Neighbours, [IsTriangulationVertex], v->Neighbours(v,Neighbour(v)));
InstallMethod(Neighbours, [IsTriangulationVertex, IsTriangulationEdge], function(v,e)
    local n;
    n := [];
    repeat
        Add(n,e);
        e := Prevopp(e);
    until e=n[1];
    return n;
end);
InstallMethod(Valency, [IsTriangulationVertex], function(v)
    local n, e0, e;
    n := 0;
    e0 := Neighbour(v);
    e := e0;
    repeat
        n := n+1;
        e := Prevopp(e);
    until e=e0;
    return n;
end);

InstallMethod(IsFake, [IsTriangulationVertex], ReturnFalse);

InstallMethod(To, [IsTriangulationEdge], e->From(Opposite(e)));
InstallMethod(Right, [IsTriangulationEdge], e->Left(Opposite(e)));
InstallMethod(Map, [IsTriangulationEdge], e->P1Path(FromPos(e),ToPos(e)));
InstallMethod(Length, [IsTriangulationEdge], e->P1Distance(FromPos(e),ToPos(e)));
InstallMethod(FromPos, [IsTriangulationEdge], e->Pos(From(e)));
InstallMethod(ToPos, [IsTriangulationEdge], e->Pos(To(e)));
InstallMethod(Pos, [IsTriangulationEdge], e->P1Barycentre(FromPos(e),ToPos(e)));
InstallMethod(GroupElement, [IsTriangulationEdge], function(e)
    e := Opposite(e);
    if HasGroupElement(e) then
        return GroupElement(e)^-1;
    else
        TryNextMethod();
    fi;
end);
InstallMethod(Opposite, [IsTriangulationEdge and HasNext], e->Prevopp(Next(e)));

InstallMethod(Radius, [IsTriangulationFace], function(f)
    local p;
    p := CallFuncList(P1Circumcentre,List(Neighbours(f),FromPos));
    SetCentre(f, p[1]); # we computed it for free
    return p[2];
end);

InstallMethod(Centre, [IsTriangulationFace], function(f)
    local p;
    p := CallFuncList(P1Circumcentre,List(Neighbours(f),FromPos));
    SetRadius(f, p[2]); # we computed it for free
    return p[1];
end);
    
InstallMethod(Pos, [IsTriangulationFace], f->P1Barycentre(List(Neighbours(f),FromPos)));
    
InstallMethod(Neighbours, [IsTriangulationFace], f->Neighbours(f,Neighbour(f)));
InstallMethod(Neighbours, [IsTriangulationFace, IsTriangulationEdge], function(f,e)
    local n;
    n := [];
    repeat
        Add(n,e);
        e := Next(e);
    until e=n[1];
    return n;
end);

InstallMethod(Valency, [IsTriangulationFace], function(f)
    local n, e0, e;
    n := 0;
    e0 := Neighbour(f);
    e := e0;
    repeat
        n := n+1;
        e := Next(e);
    until e=e0;
    while n<>3 do Error("I found a face with ",n,"<>3 adjacent edges."); od;
    
    return n;
end);

################################################################

InstallMethod(ViewString, "(IMG) for a triangulation",
        [IsSphereTriangulation],
        t->CONCAT@FR("<triangulation with ",Length(t!.v)," vertices, ",Length(t!.e)," edges and ",Length(t!.f)," faces>"));

InstallMethod(String, "(IMG) for a triangulation",
        [IsSphereTriangulation],
        t->"DelaunayTriangulation(...)");

InstallMethod(DisplayString, "(IMG) for a triangulation",
        [IsSphereTriangulation],
        function(t)
    local i, j, s;
    s := "   vertex | position                                 | neighbours\n";
    Append(s,"----------+------------------------------------------+-----------------\n");
    for i in t!.v do
        Append(s,String(CONCAT@FR("Vertex ",i!.index),9));
        Append(s," | ");
        Append(s,String(Pos(i),-40));
        Append(s," |");
        for j in Neighbours(i) do APPEND@FR(s," ",j!.index); od;
        Append(s,"\n");
    od;
    Append(s,"----------+------------------------------------------+-----------------\n");
    Append(s,"     edge | position                                 |frm to lt rt rev\n");
    Append(s,"----------+------------------------------------------+-----------------\n");
    for i in t!.e do
        Append(s,String(CONCAT@FR("Edge ",i!.index),9));
        Append(s," | ");
        Append(s,String(Pos(i),-40));
        Append(s," |");
        for j in [From,To,Left,Right,Opposite] do Append(s,String(j(i)!.index,3)); od;
        Append(s,"\n");
    od;
    Append(s,"----------+------------------------------------------+----------v-----------\n");
    Append(s,"     face | position                                 | radius   | neighbours\n");
    Append(s,"----------+------------------------------------------+----------+-----------\n");
    for i in t!.f do
        Append(s,String(CONCAT@FR("Face ",i!.index),9));
        Append(s," | ");
        Append(s,String(Pos(i),-40));
        Append(s," |");
        Append(s,String(Radius(i),-9));
        Append(s," |");
        for j in Neighbours(i) do Append(s," "); Append(s,String(j!.index)); od;
        Append(s,"\n");
    od;
    Append(s,"----------+------------------------------------------+----------+-----------\n");
    return s;
end);
INSTALLPRINTERS@(IsSphereTriangulation);

BindGlobal("LOCATE@", function(t,f0,p)
    # for an initial face f0 and a P1Point p
    # f0 is allowed to be <fail>, in which case the first face is chosen
    # returns either [face,barycentric_coords],
    #             or [face,edge,edge_coord],
    #             or [face,edge_in,edge_out,vertex]
    local y, yc, baryc, n, ymin, e, emin, i, seen;
    
    if f0=fail then f0 := t!.f[1]; fi;
    # bad, this can cost linear time, and not logarithmic.
    # we should use a "rho" method to detect loops
    seen := BlistList([1..Length(t!.f)],[]);
    repeat
        while seen[f0!.index] do
            Error("We're stuck in a loop trying to locate a face. Repent.");
        od;
        baryc := [];
        yc := [];
        ymin := 2*@.ro;
        n := Neighbours(f0);
        for e in n do
            y := P1Image(InverseP1Map(Map(e)),p);
            Add(baryc,y);
            y := SphereP1Y(y);
            Add(yc,y);
            if y<ymin then ymin := y; emin := e; fi;
        od;
        if ymin>@.reps then # inside face
            return [f0,yc];
        fi;
        if ymin>-@.reps then # on edge or vertex
            y := Filtered([1..3],i->yc[i]<@.reps);
            if Size(y)=1 then # on edge
                return [f0,emin,RealPart(P1Coordinate(baryc[y[1]]))];
            elif Size(y)=2 then # at vertex
                y := First(y,i->(i+1) mod 3 in y);
                return [f0,n[1+(y+1) mod 3],n[y],From(n[y])];
            else
                Error("There is probably a triangle with a flat angle. I'm stuck");
            fi;
        fi;
        seen[f0!.index] := true;
        f0 := Right(emin);
    until false;
end);
InstallMethod(LocateInTriangulation, "(IMG) for a triangulation and point",
        [IsSphereTriangulation, IsP1Point],
        function(t,p)
    return LOCATE@(t,fail,p)[1];
end);
InstallMethod(LocateInTriangulation, "(IMG) for a triangulation, vertex and point",
        [IsSphereTriangulation, IsTriangulationVertex, IsP1Point],
        function(t,v,p)
    return LOCATE@(t,Left(Neighbour(v)),p)[1];
end);
InstallMethod(LocateInTriangulation, "(IMG) for a triangulation, edge and point",
        [IsSphereTriangulation, IsTriangulationEdge, IsP1Point],
        function(t,e,p)
    return LOCATE@(t,Left(e),p)[1];
end);
InstallMethod(LocateInTriangulation, "(IMG) for a triangulation, face and point",
        [IsSphereTriangulation, IsTriangulationFace, IsP1Point],
        function(t,f,p)
    return LOCATE@(t,f,p)[1];
end);

BindGlobal("YRATIO@", function(a,b,c,d)
    return SphereP1Y(P1XRatio(a,b,c,d));
end);

BindGlobal("RESETVERTEX@", function(v)
#    ResetFilterObj(e, HasValency);
end);

BindGlobal("RESETEDGE@", function(e)
    ResetFilterObj(e, HasMap);
    ResetFilterObj(e, HasPos);
    ResetFilterObj(e, HasFromPos);
    ResetFilterObj(e, HasToPos);
    ResetFilterObj(e, HasLength);
end);

BindGlobal("RESETFACE@", function(f)
    ResetFilterObj(f, HasCentre);
    ResetFilterObj(f, HasPos);
    ResetFilterObj(f, HasRadius);
end);

DeclareGlobalFunction("FLIPEDGE@"); # since it calls itself recursively, first declare it
InstallGlobalFunction(FLIPEDGE@, function(e,multi)
    # flip the edge e
    # if multi, then do as many (>= 0) flips till the triangulation is Delaunay; otherwise,
    # do exactly one flip
    local a, b, p, q, bp, pa, aq, qb, f, paq, qbp, opa, oqb;

    f := Opposite(e);
    a := From(e); b := From(f);
    bp := Next(e); pa := Next(bp);
    aq := Next(f); qb := Next(aq);
    p := From(pa); q := From(qb);

    if not multi or YRATIO@(Pos(p),Pos(q),Pos(a),Pos(b))>@.rz then
        paq := Left(f); RESETFACE@(paq);
        qbp := Left(e); RESETFACE@(qbp);
        opa := Opposite(pa); oqb := Opposite(qb);
        e!.From := p; e!.To := q; e!.Next := qb; e!.Prevopp := Opposite(bp); RESETEDGE@(e);
        f!.From := q; f!.To := p; f!.Next := pa; f!.Prevopp := Opposite(aq); RESETEDGE@(f);
        qbp!.Neighbour := e;
        paq!.Neighbour := f;
        a!.Neighbour := aq;
        b!.Neighbour := bp;
        pa!.Prevopp := e; pa!.Next := aq; pa!.Left := paq; opa!.Right := paq;
        qb!.Prevopp := f; qb!.Next := bp; qb!.Left := qbp; oqb!.Right := qbp;
        aq!.Prevopp := opa; aq!.Next := f;
        bp!.Prevopp := oqb; bp!.Next := e;

        if HasGroupElement(e) then
            pa!.GroupElement := GroupElement(e)^-1*GroupElement(pa);
            Opposite(pa)!.GroupElement := pa!.GroupElement^-1;
            qb!.GroupElement := GroupElement(e)*GroupElement(qb);
            Opposite(qb)!.GroupElement := qb!.GroupElement^-1;
        fi;

        if multi then
            FLIPEDGE@(aq,true);
            FLIPEDGE@(qb,true);
        fi;
    fi;
end);

BindGlobal("CHECKTRIANGULATION@", function(t)
    local x;
    x := Filtered(t!.v,v->not IsIdenticalObj(From(Neighbour(v)),v));
    if x<>[] then return ["From(Neighbour(v)) <> v: ",x]; fi;
    x := Filtered(t!.e,e->not IsIdenticalObj(Opposite(Opposite(e)),e));
    if x<>[] then return ["Opposite(Opposite(e)) <> e: ",x]; fi;
    x := Filtered(t!.e,e->not IsIdenticalObj(Opposite(e),Prevopp(Next(e))));
    if x<>[] then return ["Opposite(e) <> Prevopp(Next(e)): ",x]; fi;
    x := Filtered(t!.e,e->not IsIdenticalObj(e,Next(Next(Next(e)))));
    if x<>[] then return ["Next(Next(Next(e))) <> e: ",x]; fi;
    x := Filtered(t!.e,e->not IsIdenticalObj(From(e),From(Prevopp(e))));
    if x<>[] then return ["From(e) <> From(Prevopp(e)): ",x]; fi;
    x := Filtered(t!.e,e->not IsIdenticalObj(Left(e),Left(Next(e))));
    if x<>[] then return ["Left(e) <> Left(Next(e)): ",x]; fi;
    x := Filtered(t!.e,e->not IsIdenticalObj(Left(e),Right(Opposite(e))));
    if x<>[] then return ["Left(e) <> Right(Opposite(e)): ",x]; fi;
    x := Filtered(t!.e,e->not IsIdenticalObj(From(e),To(Opposite(e))));
    if x<>[] then return ["From(e) <> To(Opposite(e)): ",x]; fi;
    x := Filtered(t!.f,f->not IsIdenticalObj(Left(Neighbour(f)),f));
    if x<>[] then return ["Left(Neighbour(f)) <> f: ",x]; fi;
    if ValueOption("nodelaunay")=fail then
        x := Filtered(t!.e,e->YRATIO@(ToPos(Next(e)),ToPos(Next(Opposite(e))),FromPos(e),ToPos(e))>@.rz);
        if x<>[] then return ["Delaunay condition fails: ",x]; fi;
        x := Filtered(t!.f,f->not IsIdenticalObj(LocateInTriangulation(t,f,Pos(f)),f));
        if x<>[] then return ["Locate(f,Pos(f)) <> f: ",x]; fi;
    fi;
    return true;
end);

BindGlobal("FIXDELAUNAY@", function(t)
    local idle, e;
    repeat
        idle := true;
        for e in t!.e do
            if YRATIO@(ToPos(Next(e)),ToPos(Next(Opposite(e))),FromPos(e),ToPos(e))>@.rz then
                FLIPEDGE@(e,true);
                idle := false;
            fi;
        od;
    until idle;
end);

BindGlobal("ADDTOTRIANGULATION@", function(t,f0,p,delaunay)
    # adds point p, in face f0, to triangulation t
    local v, e0, e1, e2, f, oldn;
    
    v := Objectify(TYPE_VERTEX, rec(index := Length(t!.v)+1)); Add(t!.v,v);
    SetPos(v, p);
    oldn := [];
    for e0 in Neighbours(f0) do
        Add(oldn,e0);
        f := Objectify(TYPE_FACE, rec(index := Length(t!.f)+1)); Add(t!.f,f);
        e1 := Objectify(TYPE_EDGE, rec(index := Length(t!.e)+1)); Add(t!.e,e1);
        e2 := Objectify(TYPE_EDGE, rec(index := Length(t!.e)+1)); Add(t!.e,e2);
        SetNeighbour(f, e0);
        SetFrom(e1, To(e0));
        SetFrom(e2, v);
        SetLeft(e1, f);
        SetLeft(e2, f);
        e0!.Left := f;
        Opposite(e0)!.Right := f;
        SetPrevopp(e1, Opposite(e0));
        Next(e0)!.Prevopp := e1;
        SetNext(e1, e2);
        SetNext(e2, e0);
        e0!.Next := e1;
        
        if HasGroupElement(e0) then
            SetGroupElement(e1, One(GroupElement(e0)));
            SetGroupElement(e2, One(GroupElement(e0)));
        fi;
    od;
    SetNeighbour(v, e2);
    
    for e0 in oldn do
        SetOpposite(Next(Next(e0)),Prevopp(e0));
        SetOpposite(Prevopp(e0),Next(Next(e0)));
    od;
    for e0 in oldn do
        SetPrevopp(Next(Next(e0)),Opposite(Next(e0)));
    od;
    
    t!.f[f0!.index] := Remove(t!.f);
    t!.f[f0!.index]!.index := f0!.index;
    
    # flip diagonals if needed, to preserve Delaunay condition
    if delaunay then
        for e0 in oldn do FLIPEDGE@(e0,true); od;
    fi;
    return v;
end);

InstallMethod(AddToTriangulation, [IsSphereTriangulation, IsP1Point],
        function(t,p)
    return ADDTOTRIANGULATION@(t,LOCATE@(t,fail,p)[1],p,true);
end);
InstallMethod(AddToTriangulation, [IsSphereTriangulation, IsP1Point, IsBool],
        function(t,p,delaunay)
    return ADDTOTRIANGULATION@(t,LOCATE@(t,fail,p)[1],p,delaunay);
end);
InstallMethod(AddToTriangulation, [IsSphereTriangulation, IsTriangulationFace, IsP1Point],
        function(t,f,p)
    return ADDTOTRIANGULATION@(t,f,p,true);
end);
InstallMethod(AddToTriangulation, [IsSphereTriangulation, IsTriangulationFace, IsP1Point, IsBool],
        ADDTOTRIANGULATION@);

InstallMethod(RemoveFromTriangulation, [IsSphereTriangulation, IsTriangulationVertex],
        function(t,v)
    # remove vertex v from triangulation t. flip edges as needed till v becomes
    # trivalent, then zap it.
    local y, miny, e0, e1, e2, e, f, j;
    
    while Valency(v)>3 do
        e0 := Neighbour(v);
        e1 := Prevopp(e0);
        e2 := Prevopp(e1);
        miny := @.ro/@.rz; # positive infinity
        # compute powers of circle on triplet-of-neighbours wrt v, find smallest
        repeat
            y := P1Circumcentre(ToPos(e0),ToPos(e1),ToPos(e2));
            y := y[2]^2-P1Distance(y[1],Pos(v))^2;
            if y<miny then miny := y; e := e1; fi;
            e0 := e1;
            e1 := e2;
            e2 := Prevopp(e2);
        until e0=Neighbour(v);
        FLIPEDGE@(e,false);
    od;

    # remove two faces and 6 edges, recycle face Left(e0).
    e0 := Neighbour(v); e1 := Prevopp(e0); e2 := Prevopp(e1);
    
    # deallocate faces Left(e1) and Left(e2)
    for e in [e1,e2] do
        j := Left(e)!.index;
        f := Remove(t!.f);
        if j<=Length(t!.f) then t!.f[j] := f; t!.f[j]!.index := j; fi;
    od;

    # deallocate edges in and out of v
    for e in [e0,e1,e2,Opposite(e0),Opposite(e1),Opposite(e2)] do
        j := e!.index;
        e := Remove(t!.e);
        if j<=Length(t!.e) then t!.e[j] := e; t!.e[j]!.index := j; fi;
    od;

    # deallocate vertex v
    j := v!.index;
    v := Remove(t!.v);
    if j<= Length(t!.v) then t!.v[j] := v; t!.v[j]!.index := j; fi;
    
    e0 := Next(e0); e1 := Next(e1); e2 := Next(e2);
    
    f := Left(e0); RESETFACE@(f);
    f!.Neighbour := e0;
    From(e0)!.Neighbour := e0; e0!.Left := f; Opposite(e0)!.Right := f; e0!.Next := e1; e0!.Prevopp := Opposite(e2);
    From(e1)!.Neighbour := e1; e1!.Left := f; Opposite(e1)!.Right := f; e1!.Next := e2; e1!.Prevopp := Opposite(e0);
    From(e2)!.Neighbour := e2; e2!.Left := f; Opposite(e2)!.Right := f; e2!.Next := e0; e2!.Prevopp := Opposite(e1);
end);

InstallMethod(DelaunayTriangulation, "(IMG) for a list of points and a quality",
        [IsList, IsFloat],
        function(points,quality)
    local t, i, order, n, im, p, d, idle, print;
    
    while not ForAll(points,IsP1Point) do
        Error("DelaunayTriangution: argument should be a list of points on P1");
    od;
    
    n := Length(points);
    if n=0 then points := [P1infinity]; n := 1; fi;    
    d := List(points,x->P1Distance(points[n],x));
    order := [n]; # points[order[1]] is last point, presumably infinity
    im := List(d,v->AbsoluteValue(v-@.pi/2));
    i := POSITIONID@(im,MinimumList(im));

    if im[i]>=@.pi/6 then # all points are more or less aligned to points[order[1]]
        points := ShallowCopy(points);
        i := POSITIONID@(d,MaximumList(d));
        if d[i]<@.pi/2 then # actually all points are close to points[order[1]]
            Add(points,P1Antipode(points[n]));
            Add(order,Length(points));
        else
            Add(order,i);
        fi;
        for p in [P1one,P1Point(@.i),P1Point(-@.o),P1Point(-@.i)] do
            t := MoebiusMap(points[n],points[order[2]]);
            Add(points,P1Image(t,p));
            Add(order,Length(points));
        od;
    else # points[i] is roughly at 90 degrees from points[n]
        t := MoebiusMap(points[n],points[i],P1Antipode(points[n]));
        # so t(0)=points[n], t(1)=points[i]. Try to find points close to t^-1(infty,-1,i,-i).
        i := InverseP1Map(t); im := List(points,x->P1Image(i,x));
        for p in [P1infinity,P1one,P1Point(@.i),P1Point(-@.o),P1Point(-@.i)] do
            d := List(im,x->P1Distance(x,p));
            i := POSITIONID@(d,MinimumList(d));
            if d[i]>=@.pi/6 then
                if Length(points)=n then points := ShallowCopy(points); fi;
                Add(points,P1Image(t,p));
                Add(order,Length(points));
            else
                Add(order,i);
            fi;
        od;
    fi;
    Assert(1,IsDuplicateFreeList(order),"DelaunayTriangulation couldn't create octahedron");

    Append(order,Difference([1..n],order)); # so now order[1..6] is roughly an octahedron:
    # points{order{[1..6]}} = [0,infty,1,i,-1,-i]

    # create the octahedron
    t := rec(v := List([1..6],i->Objectify(TYPE_VERTEX, rec(index := order[i]))),
             e := List([1..24],i->Objectify(TYPE_EDGE, rec(index := i))),
             f := List([1..8],i->Objectify(TYPE_FACE, rec(index := i))));
    for i in [1..6] do SetPos(t.v[i], points[order[i]]); od;
    for i in [1..24] do SetOpposite(t.e[i], t.e[i-(-1)^i]); od;
    for i in [[1,2,1,3],[2,6,4,13],[3,6,1,23],[4,5,7,17],[5,5,8,19],[6,1,7,4],
            [7,1,5,9],[8,3,6,18],[9,3,5,20],[10,4,3,14],[11,4,2,24],[12,2,3,10],
            [13,2,4,15],[14,3,3,12],[15,3,4,2],[16,6,6,8],[17,6,7,6],[18,1,6,16],
            [19,1,8,21],[20,4,5,7],[21,4,8,5],[22,5,2,11],[23,5,1,1],[24,2,2,22]] do
        SetFrom(t.e[i[1]], t.v[i[2]]);
        SetNeighbour(t.v[i[2]],t.e[i[1]]);
        SetLeft(t.e[i[1]], t.f[i[3]]);
        SetNeighbour(t.f[i[3]],t.e[i[1]]);
        SetNext(t.e[i[1]], t.e[i[4]]);
        SetPrevopp(t.e[i[4]], Opposite(t.e[i[1]]));
    od;

    # now add the other points
    for i in [7..Length(points)] do
        p := points[order[i]];
        im := LOCATE@(t,fail,points[order[i]]);
        while Length(im)=4 do # vertex
            Error("Two vertices coincide: ",p," and ",im[4]);
        od;
        ADDTOTRIANGULATION@(t,im[1],p,true);
        t.v[i]!.index := order[i];
    od;
    
    t.v{order} := ShallowCopy(t.v); # reorder the points as they were before
    
    repeat
        idle := true;
        for i in t.f do
            if HasRadius(i) then continue; fi;
            p := Radius(i) / MinimumList(List(Neighbours(i),Length));            
            if p > quality then
                ADDTOTRIANGULATION@(t,i,Centre(i),true);
                idle := false;
            fi;
        od;
    until idle;

    for i in [n+1..Length(t.v)] do # remember these are added vertices
        SetIsFake(t.v[i], true);
    od;
    
    t := Objectify(TYPE_TRIANGULATION,t);
    return t;
end);
InstallMethod(DelaunayTriangulation, "(IMG) for a list of points",
        [IsList], points->DelaunayTriangulation(points,@.rinf));

InstallGlobalFunction(EquidistributedP1Points, function(N)
    # creates a list of N points equidistributed on the sphere
    local t, x, p, r;

    p := [];

    while Length(p)<Minimum(N,10) do # add a little randomness
        x := List([1..3],i->Random([-10^5..10^5]));
        if x<>[0,0,0] then # that would be VERY unlucky
            Add(p,P1Sphere(@.ro*x));
        fi;
    od;
    if Length(p)=N then return p; fi;
    
    t := DelaunayTriangulation(p);
    r := @.pi/2;
    while Length(t!.v)<N do
        p := First(t!.f,x->Radius(x)>=r);
        if p=fail then r := r*3/4; continue; fi;
        AddToTriangulation(t,p,Centre(p));
    od;
    return List(t!.v,Pos);
end);

InstallMethod(WiggledTriangulation, [IsSphereTriangulation,IsObject],
        function(t,movement)
    # movement is either a list of new positions for the vertices
    # (in which case the vertices, edges etc. should be wiggled to
    # their new positions), or a Möbius transformation, or true.
    local r, i, j;
    r := rec(v := List([1..Length(t!.v)],i->Objectify(TYPE_VERTEX, rec(index := i))),
             e := List([1..Length(t!.e)],i->Objectify(TYPE_EDGE, rec(index := i))),
             f := List([1..Length(t!.f)],i->Objectify(TYPE_FACE, rec(index := i))));
    for i in [1..Length(t!.v)] do
        SetNeighbour(r.v[i], r.e[Neighbour(t!.v[i])!.index]);
        SetIsFake(r.v[i],IsFake(t!.v[i]));
    od;
    for i in [1..Length(t!.e)] do
        SetFrom(r.e[i], r.v[From(t!.e[i])!.index]);
        SetLeft(r.e[i], r.f[Left(t!.e[i])!.index]);
        for j in [Next,Prevopp,Opposite] do
            Setter(j)(r.e[i], r.e[j(t!.e[i])!.index]);
        od;
        if HasGroupElement(t!.e[i]) then
            SetGroupElement(r.e[i], GroupElement(t!.e[i]));
        fi;
    od;
    for i in [1..Length(t!.f)] do
        SetNeighbour(r.f[i], r.e[Neighbour(t!.f[i])!.index]);
    od;
    if IsList(movement) then
        r.wiggled := @.rz;
        for i in [1..Length(r.v)] do
            if IsFake(r.v[i]) then
                SetPos(r.v[i], Pos(t!.v[i]));
            else
                r.wiggled := r.wiggled + P1Distance(Pos(t!.v[i]), movement[i]);
                SetPos(r.v[i], movement[i]);
            fi;
        od;
    elif IsP1Map(movement) then
        for i in [1..Length(r.v)] do
            SetPos(r.v[i], P1Image(movement,Pos(t!.v[i])));
        od;
        for i in [1..Length(r.e)] do
            SetPos(r.e[i], P1Image(movement,Pos(t!.e[i])));
            SetMap(r.e[i], CompositionP1Map(movement,Map(t!.e[i])));
            if HasGroupElement(t!.e[i]) then SetGroupElement(r.e[i], GroupElement(t!.e[i])); fi;
        od;
        for i in [1..Length(r.f)] do SetPos(r.f[i], P1Image(movement,Pos(t!.f[i]))); od;
    fi;

    return Objectify(TYPE_TRIANGULATION, r);
end);

InstallMethod(ShallowCopy, [IsSphereTriangulation],
        t->WiggledTriangulation(t,fail));

InstallMethod(ClosestFaces, [IsTriangulationVertex], x->List(Neighbours(x),Left));
InstallMethod(ClosestFaces, [IsTriangulationEdge], x->[Left(x),Right(x)]);
InstallMethod(ClosestFaces, [IsTriangulationFace], x->[x]);

InstallMethod(ClosestVertices, [IsTriangulationVertex], x->[x]);
InstallMethod(ClosestVertices, [IsTriangulationEdge], x->[From(x),To(x)]);
InstallMethod(ClosestVertices, [IsTriangulationFace], x->List(Neighbours(x),From));

BindGlobal("INTERPOLATE_ARC@", function(l)
    # interpolate along points of l
    local r, i, p;
    r := ShallowCopy(l);
    i := 1;
    while i<Length(r) do
        if P1Distance(r[i],r[i+1])>@.pi/12 then
            Add(r,P1Barycentre(r[i],r[i+1]),i+1);
        else
            i := i+1;
        fi;
    od;
    return r;
end);

BindGlobal("PRINTPT@", function(f,p1p,sep,s)
    local p;
    p := sep*SphereP1(p1p);
    PrintTo(f, p[1], " ", p[2], " ", p[3], s, "\n");
end);

BindGlobal("PRINTARC@", function(f,a,col,sep)
    local j;
    a := INTERPOLATE_ARC@(a);
    PrintTo(f, "ARC ",Length(a)," ",String(col[1])," ",String(col[2])," ",String(col[3]),"\n");
    for j in a do
        PRINTPT@(f, j, sep, "");
    od;
end);

BindGlobal("PRINTEDGE@", function(f,e,col,sep)
    local a, j, t, delta, next;
    t := 0.;
    delta := 1.;
    a := [FromPos(e)];
    while t < 1. do
        next := P1Image(Map(e),P1Point(t+delta));
        if P1Distance(a[Length(a)],next) > @.pi/12 then
            delta := delta / 2.;
        elif delta < 0.5 and P1Distance(a[Length(a)],next) < @.pi/50 then
            delta := 2. * delta;
        else
            Add(a,next);
            t := t+delta;
        fi;
    od;
    a[Length(a)] := ToPos(e);
    
    PrintTo(f, "ARC ",Length(a)," ",String(col[1])," ",String(col[2])," ",String(col[3]),"\n");
    for j in a do
        PRINTPT@(f, j, sep, "");
    od;    
end);

BindGlobal("PRINTPOINTS@", function(f,t,extrapt)
    local i, x, n, arcs, labels;
    
    arcs := ValueOption("noarcs")=fail;
    labels := ValueOption("nolabels")=fail;
    
    if arcs then
        n := Length(t!.v)+Length(t!.f);
    else
        n := Number(t!.v,v->not IsFake(v));
    fi;
    PrintTo(f, "POINTS ",n+Length(extrapt),"\n");
    for i in t!.v do
        if IsFake(i) and arcs then
            PRINTPT@(f, Pos(i), @.ro, " 0.5");
        elif not IsFake(i) then
            if labels then
                x := ViewString(CleanedP1Point(Pos(i),@.p1eps));
                RemoveCharacters(x,"<>");
                x := Concatenation(" 2.0 ",x);
            else
                x := " 2.0";
            fi;
            PRINTPT@(f, Pos(i), @.ro, x);
        fi;
    od;
    if arcs then
        for i in t!.f do PRINTPT@(f, Pos(i), @.ro, " 1.0"); od;
    fi;
    for i in extrapt do PRINTPT@(f, Pos(i), @.ro, " 0.5"); od;
end);

BindGlobal("PRINTARCS@", function(f, edges, arcs, radius)
    local a, e, j, k;
    if ValueOption("noarcs")<>fail then
        PrintTo(f, "ARCS 0\n");
    else
        PrintTo(f, "ARCS ", Length(edges)+Length(arcs),"\n");
        for e in edges do
            if From(e)!.index>To(e)!.index then # print only in 1 direction
                continue;
            fi;
            j := [128,64,64];
            k := [64,128,64];
            if not (HasGroupElement(e) and IsOne(GroupElement(e))) then
                j := [255,64,64];
            else
                k := [64,255,64];
            fi;
            PRINTEDGE@(f, e, j, radius);
            PRINTARC@(f, [Pos(Left(e)),Pos(e),Pos(Right(e))], k, radius);
        od;
        for a in arcs do PRINTARC@(f, a[3], a[1], a[2], radius); od;
    fi;
end);

InstallMethod(Draw, "(IMG) for a triangulation",
        [IsSphereTriangulation],
        function(t)
    local s, f;
    s := ""; f := OUTPUTTEXTSTRING@FR(s);
    
    if ValueOption("upper")<>fail then
        PrintTo(f,"UPPER");
    fi;
    if ValueOption("lower")<>fail then
        PrintTo(f,"LOWER");
    fi;
    
    PRINTPOINTS@(f,t,[]);
    PRINTARCS@(f,t!.e,[],@.ro);
    
    Info(InfoIMG,3,"calling javaplot with:\n",s);
    JAVAPLOT@(InputTextString(s));
end);
##############################################################################

#E triangulations.gi . . . . . . . . . . . . . . . . . . . . . . . .ends here
