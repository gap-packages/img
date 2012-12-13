m := PolynomialIMGMachine(2,[1/7]);
f := StateSet(m);
i := GroupHomomorphismByImages(f,f,[f.1^(f.2*f.1),f.2^(f.1),f.3,f.4]);
SetInfoLevel(InfoFR,3);
ProfileGlobalFunctions(true);
ProfileOperationsAndMethodsOn();
RationalFunction(m^(i^15));

bg := function(f,rel)
    local n, i, g, a, s, u, v;

    g := GeneratorsOfGroup(f);
    n := Length(g);
    a := [];
    for i in [1..n] do
        s := ShallowCopy(g);
        u := Subword(rel,i,i);
        v := Subword(rel,i mod n+1,i mod n+1);
        s[Position(g,u)] := v;
        s[Position(g,v)] := u^v;
	s := GroupHomomorphismByImages(f,f,g,s);
        SetName(s,Concatenation("(",String(u),String(v),")"));
        Add(a,s);
	s := s^-1;
        SetName(s,Concatenation("(",String(v),String(u),")"));
        Add(a,s);
    od;
    return a;
end;

simpl := function(m,rel)
    local bgen, g, i, len, newm, newlen, twist, idle;

    bgen := bg(Source(m),rel);
    g := GeneratorsOfGroup(Source(m));
    len := infinity;
    twist := [];
    repeat
        idle := true;
        newm := List(bgen,t->t*m);
        newlen := List(newm,t->Sum(g,x->Length(x^t)));
        i := Position(newlen,Minimum(newlen));
        if newlen[i] < len then
            len := newlen[i];
            m := newm[i];
            Add(twist,bgen[i]);
            idle := false;
        fi;
    until idle;
    return List(Reversed(twist),Inverse);
end;
