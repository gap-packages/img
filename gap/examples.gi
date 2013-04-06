#############################################################################
##
#W examples.gi                                              Laurent Bartholdi
##
#Y Copyright (C) 2012, Laurent Bartholdi
##
#############################################################################
##
##  All interesting examples of IMG's I came through
##
#############################################################################

#############################################################################
##
#F DBRationalIMGGroup
##
BindGlobal("IMGDB@", []);
Add(IMGDB@, function()
    local e, f, z, h, H, Z;

    f := function(arg)
        local G;
        G := CallFuncList(FRGroup,arg{[4..Length(arg)]});
        SetName(G,Concatenation("IMG(",arg[2],")"));
        G!.Correspondence := arg[3];
        Add(IMGDB@,[arg[1],arg[3],G]);
    end;
    e := function(p)
        local k;
        k := AlgebraicExtension(Rationals,p);
        H := RootOfDefiningPolynomial(k);
        Z := Indeterminate(k,"z":new);
    end;

    z := Indeterminate(Rationals,"z":new);
    h := Indeterminate(Rationals,"h":new);
    Remove(IMGDB@);

    f([2,1],"z^2",z^2,
      "a=<,a>(1,2)");
    f([2,2],"z^-2",z^-2,
      "a=<,a>(1,2)");
    f([3,1],"(z-1)^2",(z-1)^2,
      "a=<,b>(1,2)","b=<,b*a/b>");
    f([3,2],"(2z-1)^2",(2*z-1)^2,
      "a=(1,2)","b=<a,b>");
    f([3,3],"(2z-1)^-2",(2*z-1)^-2,
      "a=<,a^-1/b>(1,2)","b=<a,b>");
    f([3,4],"((z-1)/(z+1))^2",((z-1)/(z+1))^2,
      "a=<,b>(1,2)","b=<a,a^-1/b>");
    f([3,5],"((z-1)/z)^2",((z-1)/z)^2,
      "a=<,b>(1,2)","b=<,a^-1/b>");
    e(2*h^3+2*h^2+2*h+1);
    f([4,1,1],"(z/h+1)^2|2h^3+2h^2+2h+1=0,h~-0.64",(Z/H+1)^2,
      "a=<,a>(1,2)","b=<,a/c>","c=<b,a/b>");
    f([4,1,2],"((z/h+1)^2|2h^3+2h^2+2h+1=0,h~-0.17-0.86i",(Z/H+1)^2,
      "a=(1,2)","b=<a,>","c=<b,c>");
    f([4,2],"((z/i+1)^2",(z/E(4)+1)^2,
      "a=(1,2)","b=<c,c*a/c>","c=<,c*b/c>");
    e(h^3+h^2+2*h+1);
    f([4,3,1],"(z/h+1)^2|h^3+h^2+2h+1=0,h~-0.56",(Z/H+1)^2,
      "a=<,a>(1,2)","b=<,a/c>","c=<b,a/c>");
    f([4,3,2],"(z/h+1)^2|h^3+h^2+2h+1=0,h~-0.21-1.30i",(Z/H+1)^2,
      "a=<,a*c/a>(1,2)","b=<a,>","c=<a*b/a,>");
    H := (1+Sqrt(-7))/4;
    f([4,4],"((z+h)/((-1-2h)z+h))^2|2h^2-h+1=0",((z+H)/((-1-2*H)*z+H))^2,
      "a=(1,2)","b=<c^-1*b*c,a>","c=<c^-1/b/a,c>");
    f([4,5],"((z-1)/(-2z-1))^2",((z-1)/(-2*z-1))^2,
      "a=<,a*c*b/c/a>(1,2)","b=<,a*c*b*a/b/c/a>","c=<a*c/a,b^-1/c/a>");
    f([4,7,1],"((z+i)/(z-i))^2",((z+E(4))/(z-E(4)))^2,
      "a=(1,2)","b=<c*a*b,b*a*b>","c=<b,b*a*c*a*b>");
    f([4,7,2],"((z+1+sqrt2)/(z-1-sqrt2))^2",((z+1+Sqrt(2))/(z-1-Sqrt(2)))^2,
      "a=(1,2)","b=<a,a*b*c>","c=<b,a*b*c*b*a>");
    f([4,7,3],"((z+1-sqrt2)/(z-1+sqrt2))^2",((z+1-Sqrt(2))/(z-1+Sqrt(2)))^2,
      "a=(1,2)","b=<b*a*c,a>","c=<b,c>");
    e(h^3+4*h^2+6*h+2);
    f([4,8,1],"(hz+1)^-2|h^3+4h^2+6h+2=0,h~-0.45",(H*Z+1)^-2,
      "a=<,a^-1/c/b>(1,2)","b=<a,>","c=<a^-1*b*a,a^-1*c*a>");
    f([4,8,2],"(hz+1)^-2|h^3+4h^2+6h+2=0,h~-1.77+1.11i",(H*Z+1)^-2,
      "a=<,c^-1/a/b>(1,2)","b=<,c^-1/a*c>","c=<b,c>");
    H := (1+Sqrt(-7))/4;
    f([4,9],"((z-h/(1+2h))/(z+h))^2|2h^2-h+1=0",((z-H/(1+2*H))/(z+H))^2,
      "a=(1,2)","b=<a*c*a,c^-1/b/a>","c=<a,a*b*a>");
    f([4,10],"(-z/2+1)^-2",(-z/2+1)^-2,
      "a=<,b^-1>(1,2)","b=<c,b^-1*a>","c=<,b>");
    f([4,11],"((z-1)/(z+1+i))^2",((z-1)/(z+1+E(4)))^2,
      "a=<,c^-1*a*b>(1,2)","b=<b^-1/a,c>","c=<,c^-1*a>");
    e(h^3-h^2+3*h+1);
    f([4,12,1],"((z+h)/(z-h))^2|h^3-h^2+3h+1=0,h~-0.29",((Z+H)/(Z-H))^2,
      "a=<,c>(1,2)","b=<a,c/b/c/a/c>","c=<,c*b/c>");
    f([4,12,2],"((z+h)/(z-h))^2|h^3-h^2+3h+1=0,h~-0.64-1.72i",((Z+H)/(Z-H))^2,
      "a=<,b*a*c/a/b>(1,2)","b=<c^-1/a/b,b*a*c*a/c/a/b>","c=<,b*a*c*b/c/a/b>");
    f([4,13],"((z-1)/(z+1/2+i/2))^2",((z-1)/(z+1/2+E(4)/2))^2,
      "a=<,b>(1,2)","b=<,a^-1/c/b>","c=<a,b*c/b>");
    f([4,14],"((z-1)/(-z/2-1))^2",((z-1)/(-z/2-1))^2,
      "a=<,b>(1,2)","b=<a,b*c/b>","c=<,b/c/b/a/b>");
    f([4,15,1],"((-3/2+sqrt5/2)z+1)^-2",((-3/2+Sqrt(5)/2)*z+1)^-2,
      "a=<,c^-1>(1,2)","b=<,c^-1*a>","c=<,c^-1*b>");
    f([4,15,2],"((-3/2-sqrt5/2)z+1)^-2",((-3/2-Sqrt(5)/2)*z+1)^-2,
      "a=<,a^-1>(1,2)","b=<,c^-1*a>","c=<b,>");
    f([4,16],"((z-1)/(z+1/2+sqrt3i/2))^2",((z-1)/(z+E(6)))^2,
      "a=<,a*b/a>(1,2)","b=<,c^-1/b/a>","c=<a,>");
end);

BindGlobal("VALUERATIONAL@", function(f,x)
    local n, d, finv;
    if x=infinity then
        finv := Value(f,1/IndeterminateOfUnivariateRationalFunction(f));
        n := Value(NumeratorOfRationalFunction(finv),0);
        d := Value(DenominatorOfRationalFunction(finv),0);
        if d=0 then return infinity; else return n/d; fi;
    else
        n := Value(NumeratorOfRationalFunction(f),x);
        d := Value(DenominatorOfRationalFunction(f),x);
        if d=0 then return infinity; else return n/d; fi;
    fi;
end);

BindGlobal("CVQUADRATICRATIONAL@", function(f)
    local roots, m, z;
    z := IndeterminateOfUnivariateRationalFunction(f);
    roots := [];
    m := function(func,exp)
        local d, e, r, k;
        d := CoefficientsOfUnivariatePolynomial(NumeratorOfRationalFunction(Derivative(func)));
        if Length(d)=2 then
            r := [-d[1]/d[2]];
        elif Length(d)=3 then
            e := d[2]^2-4*d[3]*d[1];
            if IsRat(e) or @.isc(e) then
                e := Sqrt(e);
            else
                k := AlgebraicExtension(z^2-e);
                e := RootOfDefiningPolynomial(k);
                z := Indeterminate(k,String(z):new);
                f := Value(f,z);
            fi;
            r := [(-d[2]-e)/2/d[3],(-d[2]+e)/2/d[3]];
        else return; fi;
        if exp=1 then
            UniteSet(roots,r);
        else
            UniteSet(roots,List(r,x->VALUERATIONAL@(1/z,x)));
        fi;
    end;
    m(f,1);
    m(1/Value(f,1/z),-1);
    if Length(roots)=0 then
        roots := [0,infinity];
    fi;
    if Length(roots)>=2 and not infinity in roots and AbsoluteValue(roots[1]-roots[2])<Sqrt(@.reps) then
        Remove(roots,1);
    fi;
    if Length(roots)=1 then
        Add(roots,infinity);
    fi;
    return Set(roots,x->VALUERATIONAL@(f,x));
end);

BindGlobal("CANONICALQUADRATICRATIONAL@", function(f)
    local d, m, z;
    if DegreeOfLaurentPolynomial(NumeratorOfRationalFunction(f))>2 or
       DegreeOfLaurentPolynomial(DenominatorOfRationalFunction(f))>2 then
        return fail;
    fi;
    z := IndeterminateOfUnivariateRationalFunction(f);
    m := [];
    for d in CVQUADRATICRATIONAL@(f) do
        if IsInfinity(d) then Append(m,[1,0]); else Append(m,[d,1]); fi;
    od;
    f := Value((m[2]*z-m[1])/(-m[4]*z+m[3]),Value(f,(m[3]*z+m[1])/(m[4]*z+m[2])));
    m := VALUERATIONAL@(f,infinity);
    if IsZero(m) or IsInfinity(m) then m := VALUERATIONAL@(f,0); fi;
    if IsZero(m) then
        return [z^2];
    elif IsInfinity(m) then
        return [z^-2];
    fi;
    return Set([Value(f,z*m)/m,m/Value(f,1/m/z)]);
end);

InstallGlobalFunction(DBRationalIMGGroup, function(arg)
    local i, f;
    if IsFunction(IMGDB@[1]) then IMGDB@[1](); fi; # bootstrap
    if arg=[] then
        return  IMGDB@;
    elif IsCollection(arg) and IsInt(arg[1]) then
        for i in IMGDB@ do
            if arg=i[1] then return i[3]; fi;
        od;
        return fail;
    elif IsUnivariateRationalFunction(arg[1]) then
        f := CANONICALQUADRATICRATIONAL@(arg[1]);
        for i in IMGDB@ do
            if i[2] in f then return i[3]; fi;
        od;
        return fail;
    else
        Error("Argument should be indices in Dau's table, or a rational function\n");
    fi;
end);

################################################################

InstallGlobalFunction(PoirierExamples, function(arg)
    if arg=[1] then
        return PolynomialSphereMachine(2,[1/7],[]);
    elif arg=[2] then
        return PolynomialSphereMachine(2,[],[1/2]);
    elif arg=[3,1] then
        return PolynomialSphereMachine(2,[],[5/12]);
    elif arg=[3,2] then
        return PolynomialSphereMachine(2,[],[7/12]);
    elif arg=[4,1] then
        return PolynomialSphereMachine(3,[[3/4,1/12],[1/4,7/12]],[]);
    elif arg=[4,2] then
        return PolynomialSphereMachine(3,[[7/8,5/24],[5/8,7/24]],[]);
    elif arg=[4,3] then
        return PolynomialSphereMachine(3,[[1/8,19/24],[3/8,17/24]],[]);
    elif arg=[5] then
        return PolynomialSphereMachine(3,[[3/4,1/12],[3/8,17/24]],[]);
    elif arg=[6,1] then
        return PolynomialSphereMachine(4,[],[[1/4,3/4],[1/16,13/16],[5/16,9/16]]);
    elif arg=[6,2] then
        return PolynomialSphereMachine(4,[],[[1/4,3/4],[3/16,15/16],[7/16,11/16]]);
    elif arg=[7] then
        return PolynomialSphereMachine(5,[[0,4/5],[1/5,2/5,3/5]],[[1/5,4/5]]);
    elif arg=[9,1] then
        return PolynomialSphereMachine(3,[[0,1/3],[5/9,8/9]],[]);
    elif arg=[9,2] then
        return PolynomialSphereMachine(3,[[0,1/3]],[[5/9,8/9]]);
    fi;
end);
#############################################################################

#E examples.gi. . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
