#############################################################################
##
#W p1_ieee754.gi                                            Laurent Bartholdi
##
#Y Copyright (C) 2012-2013, Laurent Bartholdi
##
#############################################################################
##
##  This file contains helper code for P1 points, in their hardware
##  representation
##
#############################################################################

################################################################
# p1 points
################################################################
InstallMethod(P1Point, [IsIEEE754FloatRep], z->C22P1POINT(NewFloat(IsPMComplex,z),1.0_z));
InstallMethod(P1Point, [IsIEEE754FloatRep,IsIEEE754FloatRep],
        function(r,i)
    return C22P1POINT(NewFloat(IsPMComplex,r,i),1.0_z);
end);
InstallMethod(P1Point, [IsPMComplex], z->C22P1POINT(z,1.0_z));
InstallMethod(P1Point, [IsPMComplex,IsPMComplex], C22P1POINT);

InstallOtherMethod(P1Point, [IsString,IsString], 1, function(re,im)
    # try first if IEEE754 P1s
    if @.isc=IsPMComplex then return STRINGS2P1POINT(re,im); fi;
    TryNextMethod();
end);

InstallMethod(P1INFINITY@, [IsIEEE754P1Point], p->C22P1POINT(1.0_z,0.0_z));

InstallMethod(P1Coordinate, [IsIEEE754P1Point], function(p)
    p := P1POINT2C2(p);
    return p[1]/p[2];
end);

InstallMethod(DisplayString, [IsIEEE754P1Point], x->P1POINT2STRING(20,x));
InstallMethod(ViewString, [IsIEEE754P1Point], x->Concatenation("<",P1POINT2STRING(5,x),">"));
InstallMethod(String, [IsIEEE754P1Point], x->Concatenation("P1Point(",P1POINT2STRING(20,x),")"));

InstallMethod(EQ, [IsIEEE754P1Point,IsIEEE754P1Point], EQ_P1POINT);
InstallMethod(LT, [IsIEEE754P1Point,IsIEEE754P1Point], LT_P1POINT);

InstallMethod(P1Sphere, [IsList], 1, # higher method than generic
        function(v)
    if ForAll(v,IsIEEE754FloatRep) then return P1SPHERE(v); fi;
    TryNextMethod();
end);

InstallMethod(SphereP1, [IsIEEE754P1Point], SPHEREP1);
InstallMethod(SphereP1Y, [IsIEEE754P1Point], SPHEREP1Y);
InstallMethod(P1Antipode, [IsIEEE754P1Point], P1ANTIPODE);
InstallMethod(CleanedP1Point, [IsIEEE754P1Point,IsFloat], CLEANEDP1POINT);

InstallMethod(P1Barycentre, [IsP1PointList], 1, # higher method than generic
        function(list)
    if ForAll(list,IsIEEE754P1Point) then
        return P1BARYCENTRE(list);
    fi;
    TryNextMethod();
end);
InstallMethod(P1Barycentre, [IsIEEE754P1Point,IsIEEE754P1Point,IsIEEE754P1Point],
        function(arg) return P1BARYCENTRE(arg); end);
        
InstallMethod(P1Midpoint, [IsIEEE754P1Point,IsIEEE754P1Point], P1MIDPOINT);

InstallMethod(P1Distance, [IsIEEE754P1Point,IsIEEE754P1Point], P1DISTANCE);

InstallMethod(XRatio, [IsIEEE754P1Point,IsIEEE754P1Point,IsIEEE754P1Point,IsIEEE754P1Point], XRATIO);

InstallMethod(P1XRatio, [IsIEEE754P1Point,IsIEEE754P1Point,IsIEEE754P1Point,IsIEEE754P1Point], P1XRATIO);

InstallMethod(P1Circumcentre, [IsIEEE754P1Point,IsIEEE754P1Point,IsIEEE754P1Point], P1CIRCUMCENTRE);

###############################################################################
### P1 maps
###############################################################################
BindGlobal("ViewStringUnivariateLaurent", function ( fam, cofs, val, name )
    local str, zero, one, mone, i, c, lc, s;
    str := "";
    zero := fam!.zeroCoefficient;
    one := fam!.oneCoefficient;
    mone := - one;
    if IsInt( name ) then
        if HasIndeterminateName( fam, name ) then
            name := IndeterminateName( fam, name );
        else
            name := Concatenation( "x_", String( name ) );
        fi;
    fi;
    if Length( cofs ) = 0 then
        return ViewString( zero );
    fi;
    lc := Length( cofs );
    if cofs[lc] = zero then
        repeat
            lc := lc - 1;
        until cofs[lc] <> zero;
    fi;
    for i in [ lc, lc - 1 .. 1 ] do
        if cofs[i] <> zero then
            c := "*";
            if i < lc then
                if IsRat( cofs[i] ) then
                    if cofs[i] = one then
                        Append( str, "+" );
                        c := "";
                    elif cofs[i] > 0 then
                        Append( str, "+" );
                        Append( str, ViewString( cofs[i] ) );
                    elif cofs[i] = mone then
                        Append( str, "-" );
                        c := "";
                    else
                        Append( str, ViewString( cofs[i] ) );
                    fi;
                elif cofs[i] = one then
                    Append( str, "+" );
                    c := "";
                elif cofs[i] = mone then
                    Append( str, "-" );
                    c := "";
                else
                    Append( str, "+" );
                    s := ViewString( cofs[i] );
                    if '+' in s or '-' in s then
                        s := Concatenation( "(", s, ")" );
                    fi;
                    Append( str, s );
                fi;
            elif cofs[i] = one then
                c := "";
            elif cofs[i] = mone then
                Append( str, "-" );
                c := "";
            else
                s := ViewString( cofs[i] );
                if not IsRat( cofs[i] ) and ('+' in s or '-' in s) then
                    s := Concatenation( "(", s, ")" );
                fi;
                Append( str, s );
            fi;
            if i + val <> 1 then
                Append( str, c );
                Append( str, name );
                if i + val <> 2 then
                    Append( str, "^" );
                    Append( str, String( i + val - 1 ) );
                fi;
            elif cofs[i] = one then
                Append( str, ViewString( one ) );
            elif cofs[i] = mone then
                Append( str, ViewString( one ) );
            fi;
        fi;
    od;
    return str;
end);

InstallMethod(ViewString, [IsIEEE754P1Map], function(map)
    local c, ind, fam;
    ind := IndeterminateNumberOfLaurentPolynomial(map);
    fam := FamilyObj(map);
   
    IsPolynomial(map);
    if IsLaurentPolynomial(map) then
        c := CoefficientsOfLaurentPolynomial(map);
        return Concatenation("<", ViewStringUnivariateLaurent(fam, c[1], c[2], ind), ">");
    else
        c := CoefficientsOfUnivariateRationalFunction(map);
        return Concatenation("<(", ViewStringUnivariateLaurent(fam, c[1], Maximum(c[3],0), ind), ")/(", ViewStringUnivariateLaurent(fam, c[2], Maximum(-c[3],0), ind), ")>");
    fi;
end);
                   
InstallMethod(P1MAPBYCOEFFICIENTS2@, [IsPMComplex,IsList,IsList],
        function(dummy,numer,denom)
    return MAT2P1MAP([[1,0],[0,1]]*1.0_z*[numer,denom]);
end);

InstallMethod(AsP1Map, [IsIEEE754FloatRep], x->AsP1Map(NewFloat(IsPMComplex,x)));

InstallMethod(AsP1Map, [IsPMComplex], x->MAT2P1MAP([[x],[One(x)]]));

InstallMethod(AsP1Map, [IsIEEE754P1Map], IdFunc);

InstallMethod(CoefficientsOfP1Map, [IsIEEE754P1Map], P1MAP2MAT);

InstallOtherMethod(CoefficientsOfUnivariateLaurentPolynomial, [IsIEEE754P1Map],
        function(map)
    local s, t, i, deg;
    map := CoefficientsOfP1Map(map);
    deg := Length(map[1]);
    s := []; t := [];
    for i in [1..2] do
        s[i] := 1; while  s[i]<=deg and IsZero(map[i][s[i]]) do
            s[i] := s[i]+1;
        od;
        t[i] := deg; while s[i]<=t[i] and IsZero(map[i][t[i]]) do
            t[i] := t[i]-1;
        od;
    od;
    while s[2]<>t[2] do
        Error(map," is not a Laurent polynomial");
    od;
    if s[1]>t[1] then
        return [[],0];
    else
        return [map[1]{[s[1]..t[1]]}/map[2][s[2]],s[1]-s[2]];
    fi;
end);

InstallMethod(CoefficientsOfUnivariateRationalFunction, [IsIEEE754P1Map],
        function(map)
    local s, t, i, deg;
    map := CoefficientsOfP1Map(map);
    deg := Length(map[1]);
    s := []; t := [];
    for i in [1..2] do
        s[i] := 1; while  s[i]<=deg and IsZero(map[i][s[i]]) do
            s[i] := s[i]+1;
        od;
        t[i] := deg; while s[i]<=t[i] and IsZero(map[i][t[i]]) do
            t[i] := t[i]-1;
        od;
    od;
    if s[1]>t[1] then
        return [[],[One(map[2][s[2]])],0];
    else
        return [map[1]{[s[1]..t[1]]},map[2]{[s[2]..t[2]]},s[1]-s[2]];
    fi;
end);

InstallMethod(IndeterminateNumberOfUnivariateRationalFunction, [IsIEEE754P1Map],
        map->@.inum);

InstallMethod(IndeterminateOfUnivariateRationalFunction, [IsIEEE754P1Map],
        map->P1z);

InstallMethod(MoebiusMap, "(IMG) for images of 0,1,infinity",
        [IsIEEE754P1Point,IsIEEE754P1Point,IsIEEE754P1Point],
        P1MAP3);

InstallMethod(MoebiusMap, "(IMG) for images of 0,infinity",
        [IsIEEE754P1Point,IsIEEE754P1Point],
        P1MAP2);

InstallMethod(MoebiusMap, "(IMG) for image of infinity",
        [IsIEEE754P1Point],
        P1MAP1);

InstallMethod(P1Path, "(IMG) Möbius transformation 0->p, 1->q, infty->P1Antipode(p)",
        [IsIEEE754P1Point,IsIEEE754P1Point],
        P1PATH);

InstallMethod(IsPolynomial, "(IMG) for a P1 map",
        [IsIEEE754P1Map],
        P1ISPOLYNOMIAL);

InstallMethod(DegreeOfP1Map, [IsIEEE754P1Map], DEGREEOFP1MAP);
        
InstallMethod(CleanedP1Map, [IsIEEE754P1Map,IsIEEE754FloatRep], CLEANEDP1MAP);

InstallOtherMethod(InverseP1Map, [IsIEEE754P1Map], INVERSEP1MAP);

InstallOtherMethod(CompositionP1Map, [IsIEEE754P1Map,IsIEEE754P1Map],
        COMPOSEP1MAP);

InstallOtherMethod(CompositionP1Map, [IsP1Map, IsP1Map, IsP1Map],
        function(map3,map2,map1)
    return CompositionP1Map(map3,CompositionP1Map(map2,map1));
end);

InstallMethod(P1Image, [IsIEEE754P1Map,IsIEEE754P1Point], P1IMAGE);

InstallMethod(P1PreImage, [IsIEEE754P1Map,IsIEEE754P1Point], P1PREIMAGE);

InstallMethod(P1PreImages, [IsIEEE754P1Map,IsIEEE754P1Point], P1PREIMAGES);

InstallMethod(CriticalPointsOfP1Map, [IsIEEE754P1Map], P1MAPCRITICALPOINTS);

InstallMethod(P1MapByZerosPoles, [IsP1PointList,IsP1PointList,IsIEEE754P1Point,IsIEEE754P1Point],
        P1MAPBYZEROSPOLES);

InstallMethod(ComplexConjugate, "(IMG) for a P1 map", [IsIEEE754P1Map], P1MAPCONJUGATE);

InstallMethod(Primitive, "(IMG) for a P1 map", [IsIEEE754P1Map], P1MAPPRIMITIVE);

InstallMethod(Derivative, "(IMG) for a P1 map", [IsIEEE754P1Map], 2, P1MAPDERIVATIVE);

InstallMethod(P1MapScaling, "(IMG) for a P1 map and a point", [IsIEEE754P1Map,IsIEEE754P1Point], P1MAPSCALING);

InstallMethod(NumeratorP1Map, "(IMG) for a P1 map", [IsIEEE754P1Map], P1MAPNUMERATOR);

InstallMethod(DenominatorP1Map, "(IMG) for a P1 map", [IsIEEE754P1Map], P1MAPDENOMINATOR);

InstallMethod(One, [IsIEEE754P1Map], f->P1MapByCoefficients([1.0_z],[1.0_z]));
InstallMethod(Zero, [IsIEEE754P1Map], f->P1MapByCoefficients([0.0_z],[1.0_z]));
        
InstallMethod(SUM, IsIdenticalObj, [IsIEEE754P1Map,IsIEEE754P1Map], P1MAP_SUM);
InstallMethod(DIFF, IsIdenticalObj, [IsIEEE754P1Map,IsIEEE754P1Map], P1MAP_DIFF);
InstallMethod(PROD, IsIdenticalObj, [IsIEEE754P1Map,IsIEEE754P1Map], P1MAP_PROD);
InstallMethod(QUO, IsIdenticalObj, [IsIEEE754P1Map,IsIEEE754P1Map], P1MAP_QUO);
InstallMethod(INV, [IsIEEE754P1Map], P1MAP_INV);
InstallMethod(AINV, [IsIEEE754P1Map], P1MAP_AINV);
Apply([SUM,DIFF,PROD,QUO],function(op)
    local make;
    make := function(x)
        if IsRat(x) and not IsInt(x) then
	    return make(NumeratorRat(x))/make(DenominatorRat(x));
        elif not IsPMComplex(x) then
            x := NewFloat(IsPMComplex,x);
        fi;
        return MAT2P1MAP([[x],[One(x)]]);
    end;
    InstallMethod(op, [IsIEEE754P1Map,IsScalar], function(x,y) return op(x,make(y)); end);
    InstallMethod(op, [IsScalar,IsIEEE754P1Map], function(x,y) return op(make(x),y); end);
    return true;
end);

InstallMethod(P1INTERSECT@, [IsIEEE754P1Map,IsIEEE754P1Map,IsIEEE754P1Map],
        P1INTERSECT_IEEE754);
    
InstallMethod(P1ROTATION@, [IsIEEE754P1Point,IsP1PointList,IsP1PointList],
        function(dummy,points,extra)
    return P1ROTATION_IEEE754(points,extra);
end);

###########################################################################

#E p1.gi . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .ends here
