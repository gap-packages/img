a := PSZAlgebra(2);
d := a.1*a.2-a.2*a.1;
v := a.2;
maxn := 20;
k := GF(2);
s := [[VectorSpace(k,[],Zero(a)),VectorSpace(k,[v])],[VectorSpace(k,[d])]];
for n in [4..2*maxn] do
  for j in [1..n-1] do
    i := n-j;
    if not IsBound(s[i]) then s[i] := []; fi;
    b := [];
    if i>1 then Append(b,List(Basis(s[i-1][j]),x->x*d-d*x)); fi;
    if j>1 then Append(b,List(Basis(s[i][j-1]),x->x*v-v*x)); fi;
    if IsOddInt(i) and IsOddInt(j) then
      Append(b,List(Basis(s[(i+1)/2][(j+1)/2]),x->x^2));
    fi;
    s[i][j] := VectorSpace(k,b,Zero(a));
  od;
od;
x := Indeterminate(Rationals,"x");
t := Indeterminate(Rationals,"t");
Sum([1..Length(s)],i->Sum([1..Length(s[i])],j->Dimension(s[i][j])*t^(i-1)*x^(j-1)));

# 2*x^24*t^14+x^23*t^15+2*x^23*t^14+x^23*t^13+x^22*t^14+x^22*t^13+x^21*t^13+x^21*t^12+x^20*t^13+2*x^20*t^12+x^20*t^11+2*x^19*t^12+2*x^19*t^11+x^18*t^12+3*x^18*t^11+2*x^18*t^10+2*x^17*t^11+3*x^17*t^10+x^17*t^9+3*x^16*t^10+2*x^16*t^9+x^15*t^10+2*x^15*t^9+x^15*t^8+x^14*t^9+x^14*t^8+x^13*t^8+x^13*t^7+x^12*t^8+2*x^12*t^7+x^12*t^6+2*x^11*t^7+2*x^11*t^6+x^10*t^7+3*x^10*t^6+x^10*t^5+x^9*t^6+x^9*t^5+x^8*t^5+x^8*t^4+x^7*t^5+2*x^7*t^4+x^7*t^3+2*x^6*t^4+x^6*t^3+x^5*t^3+x^5*t^2+x^4*t^3+2*x^4*t^2+x^3*t^2+x^3*t+x^2*t^2+x^2*t+x^2+x*t+x+t
