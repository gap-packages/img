#!/bin/sh
tail -n +4 $0 | pargap -r -q > log.$$ 2>&1
exit
################################################################
# Compute images, in parameter space, of Misiurewicz points,
# or of matings of Misiurewicz polynomials with rabbit/corabbit/airplane
#
mindenom := 8; # minimal denominator; all i/mindenom will be computed
maxdenom := 2^14; # maximal denominator
mindist := 1/10; # subdivide as long as denominator is small enough and
                 # distance between neighbouring points is >mindist
type := "airplane";

maxpcset := 16; # maximal number of post-critical points
################################################################

#ParReset();
ParEval("LoadPackage(\"fr\")");
#ParEval("SetInfoLevel(InfoFR,2)");
ParEval("EPS@fr.maxratio := MacFloat(16/10)");

################################################################
ParInstallTOPCGlobalFunction("makemeone", function(mindenom,maxdenom,mindist,maxpcset,type)
    local points, i, j, idle, c2i, i2c, obstructed, task, angle2, job;
    
    c2i := function(c)
        if IsInt(c) then return c; fi;
        return [Int(10^10*RealPart(c)),Int(10^10*ImaginaryPart(c))];
    end;
    i2c := function(i)
        if IsInt(i) then return i; fi;
        return Complex(i[1]/10^10,i[2]/10^10);
    end;
    MakeReadWriteGlobal("ErrorInner");
    ErrorInner := function(arg) JUMP_TO_CATCH(arg{[2..Length(arg)]}); end;
    if type="mandelbrot" then
	task := function(angle)
            local v;
            v := CALL_WITH_CATCH(RationalFunction,[PolynomialIMGMachine(2,[angle],false)]:param_unicritical);
	    if not v[1] then # gap error
		return 1;
	    elif IsRationalFunction(v[2]) then # z^2+c
		return c2i(Value(v[2],0));
	    elif IsRecord(v[2]) then # obstruction
		return 0;
	    else # fr error
		return 1;
	    fi;
	end;
    else # points in slice v3
        if type="rabbit" then
            angle2 := 1/7;
        elif type="airplane" then
            angle2 := 3/7;
        elif type="corabbit" then
            angle2 := 5/7;
        fi;
        obstructed := [1-angle2-1/7,1-angle2];
        task := function(angle)
            local v;
	    if angle >= obstructed[1] and angle <= obstructed[2] then
		return 0; # we know it's an obstruction
	    fi;
	    RUNTIME@fr := Runtime() + 3600*1000; # allow 1 hour
            v := CALL_WITH_CATCH(RationalFunction,[Mating(PolynomialIMGMachine(2,[angle],false),PolynomialIMGMachine(2,[angle2]))]:param_v:=3);
	    Info(InfoFR,1,"Spider converged to ",v," on ",MPI_Comm_rank());
            if not v[1] then
                return 1; # gap error
            elif IsRationalFunction(v[2]) then # 1 - (1+a)z^-1 + az^-2
                return c2i(CoefficientsOfUnivariateLaurentPolynomial(v[2])[1][1]);
            elif IsRecord(v[2]) then
		return 0; # obstruction
	    else
                return 1; # fr error
            fi;
        end;
    fi;

    points := [];

    job := [];
    # classical job
    for i in Combinations([0..maxpcset],2) do
	j := 2^i[2]-2^i[1];
	Append(job,[0..j-1]/j);
    od;
    j := AsSortedList(job);
    job := [];
    for i in [1..Length(j)] do
	if i=1 or j[i]<>j[i-1] then
	    Add(job,j[i]);
	fi;
    od;

    # Hamal Hubbard's question: only points in [2/7,1/3]
    if true then
	job := Filtered(job,angle->IsEvenInt(DenominatorRat(angle)) and angle >= 2/7 and angle <= 1/3);
    fi;

    # Xavier Buff's question: real polynomials with rabbit
    if false then
    fi;

    MasterSlave(function() # iterator
        local i, new;

	if IsBound(job) then
	    if job=[] then
		return NOTASK;
	    else
		i := Remove(job,1);
		Add(points,[i,fail]);
		return i;
	    fi;
	fi;

        i := Length(points);
        if i=0 then
            Add(points,[0,fail]);
            return 0;
        fi;
        if points[i][1]<1 and IsInt(mindenom*points[i][1]) then
            Add(points,[points[i][1]+1/mindenom,fail]);
            return points[i+1][1];
        fi;
	i := 2; while i <= Length(points) do
            if ForAll(points{[i-1,i]},p->DenominatorRat(p[1])<maxdenom and p[2]<>fail) then # something to subdivide
                if false and IS_COMPLEX(points[i-1][2]) and IS_COMPLEX(points[i][2]) and AbsoluteValue(points[i][2]-points[i-1][2])<mindist then
		    i := i+1;
		    continue;
		fi; # don't waste time here, points are very close
		new := (points[i-1][1]+points[i][1])/2;
                Add(points,[new,fail],i);
                return new; # new task
            fi;
	    i := i+1;
        od;
        return NOTASK; # done
    end,
    task, # task
      function(input,output) # check result, save locally
        points[PositionProperty(points,x->x[1]=input)][2] := i2c(output);
	Info(InfoFR,1,input," gives ",output," ",i2c(output));
        return NO_ACTION;
    end,
      Error); # update data

    return points;
end);

################################################################
points := makemeone(mindenom,maxdenom,mindist,maxpcset,type);

file := Concatenation(type,"-",String(maxpcset));
PrintTo(file,"# gnuplot data -- maxpcset=",maxpcset," type=",type,"\n");
#file := Concatenation(type,"-",String(maxdenom));
#PrintTo(file,"# gnuplot data -- maxdenom=",mindenom," maxdenom=",maxdenom," mindist=",mindist," type=",type,"\n");
lastinfinity := true;
for i in [1..Length(points)] do
    if IsInt(points[i][2]) then
	real := infinity;
	imag := infinity;
	lastinfinity := true;
    else
	if not lastinfinity and AbsoluteValue(points[i-1][2]-points[i][2])>10*mindist then
	    AppendTo(file,"infinity\t0\n"); # a jump in gnuplot
	fi;
	real := RealPart(points[i][2]);
	imag := ImaginaryPart(points[i][2]);
	lastinfinity := false;
    fi;
    AppendTo(file,real,"\t",imag,"\t",String(points[i][1]),"\t",STRING_DIGITS_MACFLOAT(6,MacFloat(points[i][1])),"\n");
od;
# hubbard.g . . . . . . . . . . . . . . . . . . . . . . . . . ends here
# recover angles:
# awk '$1=="master" {n=substr($3,1,length($3)-1); angle[n]=$4; split($4,a,"/"); if(length(a)==1)a[2]=1; angleval[n]=1.0*a[1]/a[2]} $3=="master:" {if(NF==7){printf "%.10g\t%.10g\t%s\t%g\n",substr($5,1,length($5)-1)/10000000000.0,$6/10000000000.0,angle[$1],angleval[$1]}else{print "infinity\tinfinity\t" angle[$1] "\t" angleval[$1]}}' < log.
# awk '{split($3,a,"/");if(a[2]==0)a[2]=1;b=a[1]*16384/a[2];seen[b]++} END{for(i=1;i<=11702;i++) if(seen[i]!=1) print i ",";for(i=14043;i<=16384;i++) if(seen[i]!=1) print i ","}' < rabbit-temp >

if false then

MakeReadWriteGlobal("ErrorInner");
ErrorInner := function(arg) JUMP_TO_CATCH(arg{[2..Length(arg)]}); end;

hard := [8199, 8850, 9349, 9457, 9785, 9800, 10508, 10628, 10822,
11279, 11308, 11573, 11618, 11690, 14082, 14139, 14211, 14383,
14457, 14685, 14779, 15085, 15700];

points := [];

for angle in angles2 do
  v := CALL_WITH_CATCH(RationalFunction,[Mating(PolynomialIMGMachine(2,[angle],false),PolynomialIMGMachine(2,[1/7]))]:param_v:=3);
  Info(InfoFR,1,"Angle ",angle,": spider converged to ",v);
  Add(points,[angle,v]);
od;

file := "xx";
PrintTo(file,"");
for i in [1..Length(points)] do
real := STRING_DIGITS_MACFLOAT(10,RealPart(points[i][2]));
imag := STRING_DIGITS_MACFLOAT(10,ImaginaryPart(points[i][2]));
AppendTo(file,real,"\t",imag,"\t",String(points[i][1]),"\t",STRING_DIGITS_MACFLOAT(6,MacFloat(points[i][1])),"\n");
od;

# plot [300:1200] [200:700] '< convert -negate -modulate 200 ~/math/GAP/fr/sandbox/v3.jpg avs:-' binary filetype=avs with rgbimage, '~/math/GAP/fr/sandbox/airplane-4096' using (-($1+7.15)*40+700):(-$2*160+450) with lines
fi;
a2c(x,y) = 2*(x+{0,1}*y)/(x+{0,1}*y+1)

plot [-0.7:3.75] [-1.98:1.98] '< convert -negate -colorspace Gray per3.jpg avs:-' binary filetype=avs origin=(-0.835,-1.995) dx=0.00445 dy=0.00445 with rgbimage,'rabbit-11-16384' using (real(a2c($1,$2))):(imag(a2c($1,$2))) with lines,'airplane-13'  using (real(a2c($1,$2))):(imag(a2c($1,$2))) with lines,'rabbit-11-16384'  using (real(a2c($1,-$2))):(imag(a2c($1,-$2))) with lines
set term pdfcairo size 29.7cm,21cm
set out "wittner.pdf"
replot
set term png size 1112,990
set out "wittner.png"
replot
plot [0.43:1.9] [0.5:1.98] '< convert -negate -colorspace Gray per3.jpg avs:-' binary filetype=avs origin=(-0.835,-1.995) dx=0.00445 dy=0.00445 with rgbimage,'rabbit-11-16384' using (real(a2c($1,$2))):($4 > 0.33333 && $4 < 0.666666 ? imag(a2c($1,$2)):1/0):(150+($4-0.333333)*150*3) with lines linew 2.0 palette,'airplane-13'  using (real(a2c($1,$2))):($4 > 0.142857 && $4 < 0.285715 ? imag(a2c($1,$2)):1/0):(30+($4-0.142857)*120*7) with lines linew 2.0 palette
