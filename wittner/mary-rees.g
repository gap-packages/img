x := n->4/3*(13*2^(n-3)-1)/(2^(n+3)-1);
y := n->2^(n+1)/(2^(n+3)-1);
z := n->2/7*(6*2^n+1)/(2^(n+3)-1);
Read("F.g");

for i in [3,6..30] do
    t := Runtime(); j := F(z(i)); Print("n=",i," z=",z(i), " F(z)=",j, " time=",Runtime()-t,"\n");
od;
#for i in [1..15] do
#    t := Runtime(); j := F(y(i)); Print("n=",i," y=",y(i), " F(y)=",j, " time=",Runtime()-t,"\n");
#od;
a := function(n)
    if IsOddInt(n) then
        return (9*2^(n-1)-7)/(8*2^n-1);
    else
        return (9*2^(2*n+2)+51*2^(n-1)-7)/(2^(2*n+6)-1);
    fi;
end;
b := function(n)
    if IsOddInt(n) then
        return 2/3*(5*2^n-1)/(8*2^n-1);
    else
        return 1/3*(10*2^n-1)/(8*2^n-1);
    fi;
end;
