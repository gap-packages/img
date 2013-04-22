F := function(angle)
    local gens, airplane, mating, pol, admaddresses, count, hom, basis,
    v, f, a, b, c, x, y, z, n, i, t, xx, yy ,zz, image;

    if angle <1/7 or angle >=2/7 then return fail; fi;
    airplane := PolynomialIMGMachine(2,[3/7],[]);
    if IsOddInt(DenominatorRat(angle)) then
        pol := ChangeFRMachineBasis(PolynomialIMGMachine(2,[angle],[]),(1,2));
    else
        pol := ChangeFRMachineBasis(PolynomialIMGMachine(2,[],[angle]),(1,2));
    fi;
    mating := Mating(airplane,pol);
    gens := GeneratorsOfFRMachine(mating);
    n := Length(gens) -3;
    f := StateSet(mating);
    a := f.1; b := f.2; c := f.3;
    xx:=[1];
    yy:=[]; 
    zz:=[];
    x :=One(f.1); y:=One(f.1); z:=One(f.1);
    t:=ADMADDRESSES@fr(pol);  
    if t = [fail,fail] then return fail; fi;
    admaddresses := t[2];
    image := t[1];

    t := angle;
    count :=[0,0,0,0];
    count[2]:= count[2]+1;
    for i in [2..n] do
        t:=2*t;
        if t>1 then t:=t-1; fi;  
        if t<1/7 then count[1]:= count[1]+1; continue; fi; 
        if t<2/7 then count[2]:= count[2]+1; continue; fi; 
        if t<4/7 then count[3]:= count[3]+1; continue; fi; 
        count[4]:= count[4]+1; 
    od; 

#    v := ExtRepOfObj(IMGRelator(pol));#???
#    v := v{[1,3..Length(v)-1]};#???
    for i in [1..count[4]] do
        z := z*gens[n+4-i];
    od;
    for i in [1+count[4]..count[4]+count[3]] do
        y := y*gens[n+4-i];
    od;
    for i in [1+count[4]+count[3]..count[4]+count[3]+count[2]] do
        x := x*gens[n+4-i];
    od;
    for i in [1+count[4]+count[3]+count[2]..count[4]+count[3]+count[2]+count[1]] do
        z := z*y*x*gens[n+4-i]*x^-1*y^-1;
    od;

    basis:=ShallowCopy(gens);
    basis[2]:=b^(x*c^-1*x^-1*y^-1*a^-1);
    basis[3]:=c^(x^-1*y^-1*a^-1);
    for i in [1..count[4]] do
        basis[n+4-i]:=basis[n+4-i]^(a^-1);
    od;
    for i in [1+count[4]..count[4]+count[3]] do
        basis[n+4-i]:=basis[n+4-i];
    od;
    for i in [1+count[4]+count[3]..count[4]+count[3]+count[2]] do
        basis[n+4-i]:=basis[n+4-i]^(c^-1*x^-1*y^-1*a^-1);
    od;
    for i in [1+count[4]+count[3]+count[2]..count[4]+count[3]+count[2]+count[1]] do
        basis[n+4-i]:=basis[n+4-i]^(x^-1*y^-1*a^-1);
    od;

    mating := ChangeFRMachineBasis(mating,[Inverse(y*a),Inverse(a*y*x*c*x^-1)]);
    hom := GroupHomomorphismByImages(f,f,basis,gens);
    mating := mating^hom;

    x := FreeGroup(n);
    hom := GroupHomomorphismByImages(x,f,GeneratorsOfGroup(x),gens{[4..n+3]});
    x := SubFRMachine(mating,hom);
    if x=fail then
        return mating; # could not find submachine
    fi;
    y := SupportingRays(x);
    if not IsList(y) then # could not find supporting ray
        return y;
    fi;
    y := 2-2*y[2][1][1];
    return y-Int(y);
end;
