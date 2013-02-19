LoadPackage("fr");
LoadPackage("io");

parabolic := function(angle,z)
    local s, n;
    
    if IsEvenInt(DenominatorRat(angle)) then
        return z;
    fi;
    n := First(Integers,n->n>=1 and IsInt((2^n-1)*angle));
    s := IO_PipeThrough("./parabolic",[String(RealPart(z)),String(ImaginaryPart(z)),String(n),"1.0","0.0"],"");
    return Complex(s);
end;

points := [];

makegnuplot := function(outfile)
   local f, l, real, imag, z;

   f := IO_Popen2("/usr/bin/ssh",["-C","gauss04","grep -2 gives log.total"]);
   points := [];
   while true do
       l := IO_ReadLine(f.stdout);
       if l="" or l=fail then break; fi;
       z := SplitString(l,"",WHITESPACE);
       if Length(z)<3 or z[3]<>"gives" then continue; fi;
       if z[4]="[" and Length(z)<8 then
           Append(l,IO_ReadLine(f.stdout));
           z := SplitString(l,"",WHITESPACE);
       fi;
       if z[4]<>"[" then Add(points,[EvalString(z[2]),z[4]]); continue; fi;
       Add(points,[EvalString(z[2]),Complex(z[8])]);
   od;
   IO_Close(f.stdout);
   Sort(points);
   f := OutputTextFile(outfile,false);
   PrintTo(f, "# temporary gnuplot data\n");
   for l in points do
       if IsString(l[2]) then
           real := "inf"; imag := EvalString(l[2]);
       else
           z := parabolic(l[1],l[2]);
           z := 2*z/(z+1);
           real := STRING_DIGITS_MACFLOAT(10,RealPart(z));
           imag := STRING_DIGITS_MACFLOAT(10,ImaginaryPart(z));
       fi;
       PrintTo(f,real,"\t",imag,"\t",String(l[1]),"\t",STRING_DIGITS_MACFLOAT(8,MacFloat(l[1])),"\n");
   od;
   CloseStream(f);
end;

makegnuplot("rees-data");
Exec("grep -v 'inf.*1' rees-data > rees-temp");

maxpcset := 16;
job := [];
    # classical job
for i in Combinations([0..maxpcset],2) do
    j := 2^i[2]-2^i[1];
    UniteSet(job,Filtered([0..j-1]/j,angle->IsEvenInt(DenominatorRat(angle)) and angle >= 2/7 and angle <= 1/3));
od;
