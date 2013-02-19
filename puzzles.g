#  Puzzles_Fr(machine)  finds all puzzle's systems of level 0,
# note: machine must have f1*f2*..*fn as a relation 
# and fn must correspond to the infinity loop   

#  here puzzle system is a machine with extra markings:
#   return rec(
#      machine := machine,
#      group0 := group0,           a free group of rank n               
#      group1 := group1,           a free group of rank m>n
#      embedding := embedding,     a surjective homomorphism from group0 to group1    
#      covering := cover,          a covering homomorphism from group1 to group0  
#      homomorphism := hom,        a homeomorphism from group0 to StateSet(machine)
#      globallevel := 0,           
#      levels := List([1..Length(GeneratorsOfGroup(group0))], i->[0]),    list of levels of "puzzle pieces"
#      refinest := false,
#      coordinate0 :=coordinate0,     describe cyclic ordering of puzzle piece's sides
#      coordinate1 := coordinate1     describe cyclic ordering of puzzle piece's sides        
#   );

#  some properties:
# embedding("cyclic product of generators") = "cyclic product of generators"
# embedding(covering( "cyclic product of generators" )) = "cyclic product of generators"^degree;

# for any f in group0  
# homomorphism(f) is a pre-image of homomorphism(embedding(covering(f)))

# convention: if a is a pullback of b under  puzzle.covering,  
# then  puzzle.covering(a) =   (b^"degree")^(n*"cyclic product of generators")
# and n is minimal


# for any generator f of a group "group0" embedding(f) "contains" at most one critical generator of group1
# therefore it is possible to construct a minimal machine that satisfies a given puzzle system
# the last machine and th first machines are equal if a given puzzle system is complete  


# Other functions:

# IsRefinest_Fr := function(puzzle),  check if a puzzle syztem is refinest 

#PuzzleRenormalization_Fr := function(puzzle)
   # if puzzle.refinest = true, then
   # return a homomorphism that renormalize puzzle.machine

#PuzzlePullback_Fr := function(puzzle)
# if puzzle.refinest = false make a pullback of the puzzle system (puzzle.globallevel will increase by 1)




Abelianization_Fr := function(el)
    local v, i, count;
    i := 1;
    v := ExtRepOfObj(el);
    count := [];
    for i in [2,4..Length(v)] do
        if not IsBound(count[v[i-1]]) then count[v[i-1]] := 0; fi;
        count[v[i-1]] := count[v[i-1]] + v[i];
    od;
    for i in [1..Length(count)] do
        if IsBound(count[i]) and count[i]=0 then
           Unbind(count[i]);
        fi;
    od;
    return count;
end;


StateByMachine_Fr := function(machine, el, c, set) 
# find the preimage of an element "el" that is in the region separated by "c" and "set" 
   local v, i, g, degree, bool, j, b;
#   if WreathRecursion(a)(el)[2] = [1,2] then 
#     return WreathRecursion(a)(el)[1][1]*WreathRecursion(a)(el)[1][2];;
#   fi;

   degree := Length(AlphabetOfFRObject(machine));
   g := One(GeneratorsOfFRMachine(machine)[1]);
   for i in [1..degree] do
      v := ExtRepOfObj(WreathRecursion(machine)(el)[1][i]);
      v := v{[1,3..Length(v)-1]};
      if  not IsSubset(c,v) or WreathRecursion(machine)(el)[1][i]=One(GeneratorsOfFRMachine(machine)[1]) then    
        continue;
      fi;
      bool := true; 
      for j in [1..Length(set)] do
         b:=Intersection(v,set[j]);
         if b=[] then continue; fi;
         if not IsSubset(v,set[j]) then bool:=false; fi;    

      od;


      if bool then
#Print(v,bool,el,c,"\n",c,set,"\n");
#Print(WreathRecursion(machine)(el)[1][i],"\n\n");
         return WreathRecursion(machine)(el)[1][i];
      fi;
   od;
   return One(GeneratorsOfFRMachine(machine)[1]);
end;




LoopPreimagesByFrMachine_Fr := function(machine, el)
#find all preimages of "el", where "el" is in "StateSet(machine)"
# every preimage is given in the form: [preimage,degree,coordinate,[admissible coordinates]];
# coordinate = "admissible coordinates"[1]
  local a, b, c, s, k, gens, set, group;

  gens := GeneratorsOfFRMachine(machine); 
  group := StateSet(machine);
  set := [];
  for c in Cycles(PermList(Output(machine,el)), AlphabetOfFRObject(machine)) do
      if Length(c)>1 then

          a := Transition(machine,el,c[Length(c)]); 
          b := Length(c)-1;
          while b>0   do  #  ???and Length(Transition(machine,el,c[b])*a)<Length(a) do
              a := Transition(machine,el,c[b])*a;
              b := b-1;
          od;

#          if b>0 then
#             for k in [1..b] do
#                 a := a*Transition(machine,el,c[k]);  
#             od;
#          fi;    
          Add(set,[a,Length(c), c[b+1],c]);   
          continue;
      fi;
      Add(set,[Transition(machine,el,c[1]),1,c[1],c]);

  od;



  return set; 
end;

PuzzlePartialInclusion_Fr := function(puzzle) # is not working yet
  local i, j, set, s, gens0, v, k;
  
  gens0 := GeneratorsOfGroup(puzzle.group0);
  set :=[];
  for i in [1..Length(gens0)] do
     v := ExtRepOfObj(Image(puzzle.embedding,gens0[i]));
     v := Set(v{[1,3..Length(v)-1]});
     Add(set,[v,i]);
  od; 
  set := Set(set);
  s := [set[1]];
  j :=[s[1]];
  Print(j);
  for i in [2..Length(gens0)] do
    if IsSubset(j[Length(j)][1],set[i][1]) then
       Add(j[Length(j)],set[i]);
       Add(j,j[Length(j)][Length(j[Length(j)])]); 
       continue;
    fi; 
    v := false;
    while not v do
       Remove(j[Length(j)],1); 
       Remove(j);
       if IsSubset(j[Length(j)][1],set[i][1]) then
          Add(j[Length(j)],set[i]);
          Add(j,j[Length(j)][Length(j[Length(j)])]); 
          v := true;
       fi;
    od;  
Print(s,"\n\n");  

  od;
  while Length(j)> 0 do
     Remove(j[Length(j)],1);
     Remove(j);
  od; 
  Print(s);



  return fail;
end; 





IsGraterInAbelianization_Fr := function(el1, el2) 
    local v, i, count;
    i := 1;
    v := ExtRepOfObj(el1);
    count := [];
    for i in [2,4..Length(v)] do
        if not IsBound(count[v[i-1]]) then count[v[i-1]] := 0; fi;
        count[v[i-1]] := count[v[i-1]] + v[i];
    od;
    v := ExtRepOfObj(el2);
    for i in [2,4..Length(v)] do
        if not IsBound(count[v[i-1]]) then count[v[i-1]] := 0; fi;
        count[v[i-1]] := count[v[i-1]] - v[i];
    od;
    for i in count do
       if i < 0 then return false; fi;  
    od;
    return true;
end;






PullbackHomByMachine_Fr := function(machine, group0, hom) 
# for an surjective "hom" from "group" to "StateSet(machine)" 
# it constructs a "puzzle system", if it is possible 
 local i, degree, gens0, set, j, t, k, group, n, gens1, group1, s1, s0, q, cover, embedding, coordinate0, coordinate1, inclus;
 degree := Length(AlphabetOfFRObject(machine));
 gens0 := GeneratorsOfGroup(group0);  
 group := StateSet(machine);
 set:=[];
 coordinate0 := List([1..Length(gens0)],i-> i);
 Add(coordinate0, 1); 
 for i in [1..Length(gens0)] do
   t := LoopPreimagesByFrMachine_Fr(machine,Image(hom,gens0[i]));

   for j in t do
      if not j[1]=One(group) then
         for k in [1..Length(gens0)] do

            if IsGraterInAbelianization_Fr(Image(hom,gens0[k]),j[1]) then
                Add(set,[j[3],i,j[1],j[4],k]);
            fi;
         od;
         continue;
      fi;
      Add(set,[j[3],i,j[1],j[4],false]);
   od;


 od;

 set :=Set(set);
 n := Length(set);
 group1 := FreeGroup(n);
 gens1 := GeneratorsOfGroup(group1);
 coordinate1 := List([1..degree*Length(coordinate0)],i-> 0);
 for i in [1..Length(set)] do    for j in [1..Length(set)] do
     if not set[i][5] =false or Length(set[j][4])=1 or set[j][5]=1 
     then
        continue;
     fi;
     if  set[i][2] < set[j][4][Length(set[j][4])] and set[i][2] > set[j][4][1]
       or set[i][2] = set[j][4][Length(set[j][4])] and set[i][1]<set[j][1]
       or set[i][2] = set[j][4][1] and set[i][1]>set[j][1] 
     then
          set[i][5] := set[j][5];
     fi;

  od;   od;
  i := 1;
  while i < Length(set) and (set[i][5]= false or set[i][5]=1) do
     set[i][5] := 1;
     i := i+1;
  od;

  i := Length(set);
  while i >1 and (set[i][5]= false or set[i][5]=1) do
     set[i][5] := 1;
     i := i-1;
  od;
  q := false;
  j:=1;
  for i in [1..Length(set)] do
     if set[i][5]= false then
        set[i][5] := j;
        q := true; 
        continue;
     fi;
     if q then
        if not j=set[i][5] then 
            Print("Unable to construct puzzle\n");
            return fail;
        fi; 
     fi; 
     j := set[i][5];
     q := false;
  od;



  s0 := List([1..Length(gens0)], i->One(group1));
  s1 := List([1..Length(gens1)], i->One(group0));
  q := One(group0);
  for i in [1..Length(gens0)] do
    q:= gens0[i]*q;
  od;
  for i in [1..Length(gens1)] do
#Print(set[i],"\n");
     s1[i]:= q^(set[i][1]-1)*gens0[set[i][2]]^Length(set[i][4])*q^(1-set[i][1]);
     for j in set[i][4] do
        coordinate1[set[i][2]+Length(coordinate0)*(j-1)]:=i;
     od;
     if set[i][2]= 1 then
       for j in set[i][4] do
          if j >1 then
            coordinate1[Length(coordinate0)*(j-1)]:=i;
          else;
            coordinate1[Length(coordinate0)*degree]:=i;     
          fi; 
       od;
     fi;
  od;

  cover := GroupHomomorphismByImages(group1, group0, gens1, s1); 

  i := 1;
  while i <= Length(gens1) and set[i][5]=1 do
    s0[set[i][5]]:= gens1[i]*s0[set[i][5]];
    i:=i+1;
  od;
  while i <= Length(gens1) and set[i][5]>1 do
     s0[set[i][5]]:= gens1[i]*s0[set[i][5]];
     i:=i+1;
  od;  
  q := One(group1);
  for j in [2..Length(gens0)] do
     q:= s0[j]*q;
  od;
  while i <= Length(gens1)  do
     s0[1]:= q^-1*gens1[i]*q*s0[1];
     i:=i+1;
  od;  

  embedding := GroupHomomorphismByImages(group0, group1, gens0, s0); 





   return rec(
      machine := machine,
      group0 := group0,
      group1 := group1,
      embedding := embedding,
      covering := cover,
      homomorphism := hom,
      globallevel := 0,
      levels := List([1..Length(gens0)], i->[0]),
      refinest := false,
      coordinate0 :=coordinate0,
      coordinate1 := coordinate1  
   );

end;

PullbackHomOfPuzzle_Fr := function(puzzle) 
  # construct a homomorphism from puzzle.group1 to State StateSet(machine)
  # that is compatible with puzzle.homomorphism
  local gens0, gens1, set0, set1, i, j, k, t, set, hom;
  gens0 := GeneratorsOfGroup(puzzle.group0);
  gens1 := GeneratorsOfGroup(puzzle.group1);
  set0 := [];
  set1 := List([1..Length(gens1)],i->One(StateSet(puzzle.machine)));
  for i in [1..Length(gens1)] do
     j := Image(puzzle.homomorphism,Image(puzzle.covering,gens1[i]));
     set1[i] := WreathRecursion(puzzle.machine)(j)[1][1]*set1[i];
  od; 
  

 return  GroupHomomorphismByImages(puzzle.group1,StateSet(puzzle.machine),gens1,set1);
end;


IsRefinest_Fr := function(puzzle)
# check if a puzzle syztem is refinest


 local i, j, set, gens0, gens1, s, hom;
 gens0 := GeneratorsOfGroup(puzzle.group0);
 gens1 := GeneratorsOfGroup(puzzle.group1);
 hom := PullbackHomOfPuzzle_Fr(puzzle);
 set := List([1..Length(gens0)], i ->0);
 for i in [1..Length(gens1)] do
    if not Image(hom,gens1[i])=One(StateSet(puzzle.machine)) then
        j := Abelianization_Fr(Image(puzzle.covering,gens1[i])); 
        set[Position(j,Maximum(j))] := set[Position(j,Maximum(j))] + 1; 
    fi;
 od;

for i in set do
   if i>1 then
      puzzle.refinest := false; 
      return false; 
   fi;
od;
 
 puzzle.refinest := true; 
 return true; 
end;

CreatePuzzle_Fr := function(machine, set) # creates puzzle of level 0; 
  local  i, j, q, group, subgroup, group0, group1, gens, gens0, gens1, cover, embedding, hom,
        s, s0, s1, set2, v, n;
  set2 :=[];
#Print(set,"\n");
  for i in [1..Length(set)] do
     q :=[];
     for j in [1..Length(set[i])] do
         v := ExtRepOfObj(set[i][j]);
         v := Set(v{[1,3..Length(v)-1]});
         Add(q,v);
     od;
     set2:=Union(set2,q);
  od;

  if Length(set2)>1 then 
     Remove(set2,2);  # in this case set2[2] is the union of the rest
  fi; 


  gens := GeneratorsOfFRMachine(machine); 
  n := Length(gens)-1;
  group := StateSet(machine);
  s0 := [];
  s1 := [];
  s  := List([1..Length(set2)+2], i ->One(group));

  j := 0;

  for i in [1..n] do 
     if j<=Length(set2)-1 and not i in set2[j+1] then
         s[j+1]:=gens[i]*s[j+1]; 
         continue;
     fi;
     if j=Length(set2) and i in set2[j] then
         s[j+1]:=gens[i]*s[j+1]; 
         continue;
     fi;
     if j<=Length(set2) then j:=j+1; fi;
     s[j+1]:=gens[i]*s[j+1]; 
  od;
  q := One(group); 
  for j in [2..Length(set2)+1] do
    q:=s[j]*q;
  od;
  s[1] := q^-1*s[Length(s)]*q*s[1];
  Remove(s); 

  group0 := FreeGroup(Length(s));
  gens0 := GeneratorsOfGroup(group0); 

  hom := GroupHomomorphismByImages(group0,group,gens0,s);




  return PullbackHomByMachine_Fr(machine, group0, hom);
end;


Puzzles_Fr := function(machine) # we assume that the last generator is the adding element (a loop around infinity)


  local a, v, b, s, ss, c, t, n, image, admaddresses, criticalpoints, criticalvalues, gens, i, j, k, m, group, degree, puzzles;
  gens := GeneratorsOfFRMachine(machine); 
  group := StateSet(machine);
  t:=ADMADDRESSES_FR(machine);  
  if t = [fail,fail] then return fail; fi;
  admaddresses := t[2];
  image := t[1];  
  n:=Length(gens)-1;
  criticalpoints :=[];
  criticalvalues :=[];
  puzzles := [];

# looking for "critical points and values" 
  for s in [1..Length(gens)-1] do
        for c in Cycles(PermList(Output(machine,gens[s])), AlphabetOfFRObject(machine)) do
            if Length(c)>1 then
                ss:=[];
                t:=[]; 
                for k in [1..Length(c)] do
                   v := ExtRepOfObj(Transition(machine,gens[s],c[k]));
                   v := v{[1,3..Length(v)-1]};
                   if Length(v)<>Length(Set(v)) then return fail; Print("It is not an IMG machine!!"); fi;
                   Add(ss,v);
                   Add(t,Transition(machine,gens[s],c[k]));  
                od;    
                Add(criticalpoints,ss);
                Add(criticalvalues,t);   

            fi;

        od;
  od;


 # looking for \alpha fixed points 


 for i in [1..Length(criticalvalues)] do
   for k in [1..Length(criticalvalues[i])-1] do 
      if criticalvalues[i][k]= One(group) then continue; fi;
      ss:=[];
      t :=[criticalvalues[i][k]];
      for j in [1..Length(criticalvalues)] do
         if j=i then continue; fi;
         if IsSubset(criticalpoints[i][Length(criticalpoints[k])],criticalpoints[j][Length(criticalpoints[j])]) then
             Add(ss,criticalpoints[j][Length(criticalpoints[j])]);
             Add(t,criticalvalues[j][Length(criticalpoints[j])]);
         fi;
      od; 
      b:=[];
      for m in t do
         a :=[];
         s := StateByMachine_Fr(machine,m,criticalpoints[i][k],ss);
         while not s in a do
             Add(a,s);
             s:=  StateByMachine_Fr(machine,s,criticalpoints[i][k],ss);     
         od; 
         if s = One(group) then continue; fi;
         c:= Position(a,s);
         Add(b,a{[c,c+1..Length(a)]}); 
      od;

      if b = [] then continue; fi;

#Print(Preimage_Fr(machine,criticalvalues[i][k],criticalpoints[i][Length(criticalpoints[k])],ss));
    #  Print(i,t,"\n",Puzzle_Fr(machine,b),"\n\n");
       Add(puzzles, CreatePuzzle_Fr(machine,b));

   od;
  od;




#cover, embedding



#  Print(criticalpoints,criticalvalues);

 return puzzles;
end; 


PuzzleRenormalization_Fr := function(puzzle)
   # if puzzle.refinest = true, then
   # return a homomorphism that renormalize puzzle.machine

   local gens0, gens, gr;
  
   IsRefinest_Fr(puzzle);
   if not puzzle.refinest = true then
      return fail;
   fi;
   gens0 := GeneratorsOfGroup(puzzle.group0);
   gr := FreeGroup(Length(gens0));
  
   return  GroupHomomorphismByImages(gr,StateSet(puzzle.machine),GeneratorsOfGroup(gr),Image(puzzle.homomorphism,gens0));

end;

PartialPullbackOfOrdering_Fr := function(group0, subset, group1, embedding) 
# return a list of generators
# subset is in GeneratorsOfGroup(group0);
# embedding("cyclic product of generators") = "cyclic product of generators"


   local  i, j, sub0, sub1, gens0, gens1, n, m, gr, subgr, s, hom, v, order;
   sub0  := List([1..Length(subset)], i->ExtRepOfObj(subset[i])[1]);
   gens0 := GeneratorsOfGroup(group0);
   gens1 := GeneratorsOfGroup(group1); 
   sub1 :=  List([1..Length(subset)],i->[Image(embedding,subset[i]),i]); 
   sub1 := Set(sub1);
   order := List([1..Length(gens1)],i->gens1[i]);
   v := [];
   for s in sub1 do
       Add(order,s, Position(order,gens1[ExtRepOfObj(s[1])[1]]));  
       j :=  Abelianization_Fr(s[1]);
       for i in [1..Length(j)] do
         if not IsBound(j[i]) then continue; fi;
         Add(v,i);         
       od;
   od;
   for i in v do 
     Remove(order,Position(order,gens1[i]));
   od;   
########################################################################################### 
#searching for bugs:
   j := One(group1);
   for i in order do 
     if i in gens1 then
       j := i*j; 
     else
       j := i[1]*j;
     fi;
   od;
   for i in gens1 do 
     j := j*i^-1; 
   od;  
  if not  j = One(group1) then
      Print("there is a bug in PartialPullbackOfOrdering_Fr");
      return fail;
  fi;

###########################################################################################
  for i in [1..Length(order)] do 
     if not order[i] in gens1 then
        order[i] := subset[order[i][2]];
     fi;
   od;



return order;




########################################################################################### 
#Print(sub1,"\n\n",order,"\n\n",v); return fail;
#   j := One(group1);
#   for i in subset do
#       j := j* Image(embedding, i);
#   od;
#   j :=  Abelianization_Fr(j);
#  
#   for i in [1..Length(j)] do
#       if not IsBound(j[i]) then continue; fi;
#       Remove(order,Position(order,i));
#   od;
#   n := Length(sub0);
#   m := Length(sub1);

###########################################################################################
#   gr := FreeGroup(n+m);
#   s := [];
#   for i in sub0 do
#     Add(s,Image(embedding,gens0[i]));
#   od;
#   for i in sub1 do
#     Add(s,gens1[i]);
#   od;
#   subgr := Subgroup(group1,s);
#
#
#
#   j := One(group1);
#
#   for i in gens1 do
#       j := i*j;
#   od;
#Print(j, "\n",Length(s),"\n",subgr,"\n", j in subgr); return fail;


###################################################################################
#   i := ShortGroupWordInSet(subgr,j,Length(s)); 
###########################################################################################
   if Length(i)= 1 then return fail; fi;
   i := i[2];
   v := ExtRepOfObj(i);  
   v := v{[1,3..Length(v)-1]};
   j := [];
   for i in v do
     if i<=n then
        Add(j,gens0[sub0[i]],1);
     else
        Add(j,gens1[sub1[i-n]],1);
     fi;
   od;
   return j;
end;


PuzzlePullback_Fr := function(puzzle)
# if puzzle.refinest = false make a pullback of the puzzle system
 local a, b, q, i, j, v,  sets,  set1,  set,  gens0, gens1, newgens1, newgens0, order, ss, info, image, preim0, preimset0, preimset, 
 degree, s, s2, hom,  n, m, k, newgroup0, newhomomorphism, newcovering, newembedding, t, newgroup1, iteration, newlevels, inclus, bool;
 IsRefinest_Fr(puzzle);
 if puzzle.refinest then return false; fi;

 degree := Length(AlphabetOfFRObject(puzzle.machine)); 
 gens0 := GeneratorsOfGroup(puzzle.group0);
 gens1 := GeneratorsOfGroup(puzzle.group1);
 hom := PullbackHomOfPuzzle_Fr(puzzle);
 
 sets := List([1..Length(gens0)], i ->[]);
 set := List([1..Length(gens0)], i ->0);
 image := List([1..Length(gens1)], i ->0);

# searching for reducible puzzle pieces: 

 for i in [1..Length(gens0)] do     
    v := Abelianization_Fr(Image(puzzle.embedding,gens0[i]));
  
    for j in [1..Length(v)] do
       if not IsBound(v[j]) then continue; fi;
       s := Image(hom,gens1[j]);
       Add(sets[i], j);
       image[j] := i;
       if not s=One(StateSet(puzzle.machine)) then       
           set[i] := set[i] + 1; 
       fi;       
    od;
 od;

 s := [];

 for i in [1..Length(gens0)] do
    if set[i] = 1 then 
       Add(s,gens0[i]);
    fi;
  od;

  order := PartialPullbackOfOrdering_Fr(puzzle.group0, s, puzzle.group1, puzzle.embedding); 

  n := Length(order);  
  newgroup0 := FreeGroup(n);
  newgens0 := GeneratorsOfGroup(newgroup0);

# use puzzle.machine: 
 t :=[];
 m := 0; 
 set1 := []; 
 for i in order do
   if i in gens0 then
      Add(t,Image(puzzle.homomorphism, i));
      k := LoopPreimagesByFrMachine_Fr(puzzle.machine,t[Length(t)]);
      Add(set1,k);
      m := m+Length(k);
   else
      Add(t,Image(hom, i));
      k := LoopPreimagesByFrMachine_Fr(puzzle.machine,t[Length(t)]);
      Add(set1,k);
      m := m+Length(k);
   fi;
 od; 
preim0 := List([1..Length(gens0)],i->[]);
 for i in [1..Length(gens0)] do
    k := LoopPreimagesByFrMachine_Fr(puzzle.machine,Image(puzzle.homomorphism,gens0[i]));
    for j in k do
       Add(preim0[i],j[4]); 
    od;
    preim0[i] := Set(preim0[i]);
 od; 





#  create newcovering:

 newhomomorphism := GroupHomomorphismByImages(newgroup0,StateSet(puzzle.machine),GeneratorsOfGroup(newgroup0),t);

 newgroup1 := FreeGroup(m);
 newgens1 := GeneratorsOfGroup(newgroup1);
 
 q := One(newgroup0); 
 for j in [1..Length(newgens0)] do
   q:=newgens0[j]*q;
 od;

 j :=1;
 s :=[];
 info := []; 

 for t in [1..degree] do
    for i in [1..Length(set1)] do  #??? make a copy of set1?
      if Length(set1[i])=0 then continue; fi;
      if  not set1[i][1][3] = t then continue; fi;
      Add(s, q^(t-1)*newgens0[i]^set1[i][1][2]*q^(1-t));
      Add(info,[i,[t]]); #????Add(info,[i,set1[i][1][4]]);
      Remove(set1[i],1);
    od;
  od;

 newcovering := GroupHomomorphismByImages(newgroup1,newgroup0,newgens1,s);   


 newlevels :=[];
 q := [];
 for i in order do
    if i in gens0 then
       Add(q, Image(puzzle.embedding, gens0[Position(gens0,i)]));
       Add(newlevels,ShallowCopy(puzzle.levels[Position(gens0,i)]));
    else
       j := Abelianization_Fr(Image(puzzle.covering,gens1[Position(gens1,i)]));
       Add(q, gens1[Position(j, Maximum(j))]);
       Add(newlevels,ShallowCopy(puzzle.levels[Position(
             j, Maximum(j)
           )])); 
       Add(newlevels[Length(newlevels)],puzzle.globallevel+1);
    fi; 

 od;


a :=[];
for i in order do
   if i in puzzle.group1 then
     Add(a,i);
   else
     Add(a,Image(puzzle.embedding,i));
   fi;
od;
inclus := GroupHomomorphismByImages(newgroup0,puzzle.group1,newgens0,a);

 

  preimset0:= List([1..Length(gens0)],i->[]);
  for i in [1..Length(gens1)] do
    k := Abelianization_Fr(Image(puzzle.covering,gens1[i]));
    a := preimset0[Position(k,Maximum(k))];
    for j in [1..Maximum(k)] do 
       Add(a,i);
    od;
  od;

  preimset := List([1..Length(order)],i->[]);
  for i in [1..Length(order)] do
     if order[i] in gens0 then 
       preimset[i]:=ShallowCopy(preimset0[Position(gens0,order[i])]);
     else
       preimset[i]:=ShallowCopy(preimset0[image[Position(gens1,order[i])]]);
     fi;
  od;
  

  set := List([1..n], i ->One(newgroup1));
  set1 := List([1..n], i ->One(newgroup1)); 
  #set2 := One(newgens1); #set2 := List([1..n], i ->One(newgroup1)); 
  
  v:=List([1..n], i ->0);
  ss:= [1];   
  s2 := List([1..Length(gens1)], i ->0);
  k := 2;

  for i in [1..m] do

      j :=info[i][1];
      s := j; 
      q := 2;
      bool := false;
      while not bool do 
        for t in [1..Length(preimset[s])] do     
           if bool then continue; fi;
           j := preimset[s][t];          
           if gens1[j] in order then
              j := Position(order,gens1[j]);
           else
             j := Position(order,gens0[image[j]]);
           fi;
           if j=ss[Length(ss)] then 
               bool := true;
               a := t; 
           fi;
         od;
         if bool then
            Remove(preimset[s],a); 
         fi;
         if q=1 and bool then
           k:=k+1;
           continue;  
         fi;
         if bool then continue; fi;
         if q=2 and k<=n then
            Add(ss,k);
            q:=1;
            continue;
         fi;
         if q=1 then
            Remove(ss);
            Remove(ss);   
            q:=0;
            continue;
         fi; 
         Remove(ss);
         if Length(ss)=0 then
               Print("bug"); return fail;
         fi;
      od;
     
      if v[j]>0 and v[j]+1<i then
         a := List([1..i-v[j]-1], b->s2[v[j]+b]);   

         for q in [Minimum(a)..Maximum(a)] do
             set1[j]:=set[q]*set1[j];
         od;
      fi;
      v[j] := i;
      s2[i] := j; 
      set[j] := set1[j]^-1*newgens1[i]*set1[j]* set[j];
  od;



  newembedding := GroupHomomorphismByImages(newgroup0,newgroup1,newgens0,set);

###################################

#Print(inclus);
for i in newgens0 do

 #  j := LoopPreimagesByFrMachine_Fr(puzzle.machine,Image(newhomomorphism,Image(newcovering,Image(newembedding,i) )))[1][1];
 #  q := Image(newhomomorphism,i);
   j := Image(inclus,Image(newcovering,Image(newembedding,i) ));
   q := Image(puzzle.embedding,Image(puzzle.covering, Image(inclus,i)));
 
   if not j=q then
     Print(i,"\n\n",Abelianization_Fr(j),"\n",Abelianization_Fr(q));#
#  Print("\n\n",Image(puzzle.covering, Image(inclus,i)),"\n\n",Image(puzzle.embedding,Image(puzzle.covering, Image(inclus,i))),"\n");#
#     Print("\n\n\n\n",j,"\n\n",q,"\n\n\n",Image(newembedding,i),"\n\n",Image(newembedding,i),"mist\n\n\n");#Image(newembedding,i)
Print("\nthere is a bug somewhere\n");
   fi;
 
od;
#Print(Image(newembedding,newgens0[1]),"\n",Abelianization_Fr(Image(newembedding,newgens0[1])) );
##################################
  v := One(newgroup1);
   for i in [1..m] do
     v := newgens1[i]*v;
   od;
   q := One(newgroup0);
   for i in [1..n] do
     q := newgens0[i]*q;
   od;

#Print("\n\n\n", Image(newembedding,newgens0[14]),"\n",Image(newembedding,newgens0[15]));
#Print(ss);



##################????? search for bugs:
 
   if not  v^degree = Image(newembedding, Image(newcovering,v)) then
#Display(Image(newembedding, Image(newcovering,v))); 
     Print(Abelianization_Fr(Image(newembedding, Image(newcovering,v))),"there is a bug somewhere");
     return fail;
   fi;
##################?????  
  

 return  rec(
      machine := puzzle.machine,
      group0 := newgroup0,              
      group1 := newgroup1,           
      embedding := newembedding,    
      covering := newcovering,          
      homomorphism := newhomomorphism,       
      globallevel := puzzle.globallevel +1,           
      levels := newlevels,   
      refinest := false           
   ); 
end;


CreateIMGMachineFromPuzzle_Fr := function(puzzle)
# create a new img machine that is approximate puzzle.machine
# by the first return of critical points

 local i, j, v, q, k, a, b,  sets,  set1, set0, set2, set3, set4, set,  gens0, gens1, image, hom,  
 degree, s, inclus, bool;
 #IsRefinest_Fr(puzzle);
 #if puzzle.refinest then return false; fi;

 degree := Length(AlphabetOfFRObject(puzzle.machine)); 
 gens0 := GeneratorsOfGroup(puzzle.group0);
 gens1 := GeneratorsOfGroup(puzzle.group1);

 sets := List([1..Length(gens0)], i ->[]);
 set := [];
 inclus := List([1..Length(gens1)], i ->0);
 image := List([1..Length(gens1)], i ->0);
 hom := PullbackHomOfPuzzle_Fr(puzzle);

 for i in [1..Length(gens0)] do     
    v := Abelianization_Fr(Image(puzzle.embedding,gens0[i]));
  
    for j in [1..Length(v)] do
       if not IsBound(v[j]) then continue; fi;
       if not Image(hom,gens1[j])=One(StateSet(puzzle.machine)) then   
          Add(sets[i], j);          
       fi;
       inclus[j] := i; 
    od;
 od;

 set1 :=[];
 set0 :=[];

 for i in [1..Length(gens1)] do     
    v := Abelianization_Fr(Image(puzzle.covering,gens1[i]));
    j := Maximum(v);
    image[i]:= Position(v,j);
    if j > 1 then 
      Add(set,[i,j]);

      v := LoopPreimagesByFrMachine_Fr(puzzle.machine, Image(puzzle.homomorphism, gens0[image[i]]));
#Print(Int((Length(ExtRepOfObj(Image(puzzle.covering,gens1[i])))/2+1)/Length(gens0)/2));
      for k in v do
#Print(k[4],"\n");
          a := Int((Length(ExtRepOfObj(Image(puzzle.covering,gens1[i])))/2+1)/Length(gens0)/2) +1;
          if  a in k[4] 
          then
             Add(set1,k[4]-a);
          fi;
      od;
    fi;  
 od;
#Print(set1,"\n");
#Print(sets,"\nimage",image,"\ninclus",inclus,"\nset",set,"\n\n"); 
#Print(inclus,"\n",image,"\n",sets,"\n\n");
  set2 :=[]; 
  set3 :=[];
  set4 :=[];

  for i in set do 
     q := [sets[image[i[1]]]];
     j:= i[1];#Image(puzzle.embedding,gens0[image[i[1]]]);
     v := sets[image[j]];
     
     while not i[1] in v do
        a := []; 
        for j in v do
           a := Union(a,sets[image[j]]); 
        od;  
        Add(q,a,1);
        v := a;
     od;   
     Add(set0,q);
  od;

  for i in [1..Length(set0)] do
    q:=[];
    a := set[i][1];#Position(set1[i][1],Maximum(set1[i][1]));
    Add(q,a,1);
    for j in [2..Length(set0[i])] do
       bool := false;
       for k in [1..Length(set0[i][j])]do
           if a in sets[image[set0[i][j][k]]]  then
               if bool then 
                   Print("IMGMachine is not unique");
                   continue; 
               fi;  
               bool := true;
               a := set0[i][j][k];
               Add(q,a,1);
           fi; 

       od;

    od;
    Add(q,set[i][1],1);
    Add(set2,q);
  od;





  for i in [1..Length(set2)] do
     q:=[];

     for j in [1..Length(set2[i])-1] do
       Add(q,
              Int( (Length(ExtRepOfObj(Image(puzzle.covering,gens1[set2[i][j]])))/2+1)/Length(gens0)/2 ),
            1);
           
     od;
     Add(set3,q);
  od;

  for i in [1..Length(set2)] do

     q:=0;
     a:=1;
     for j in [1..Length(set3[i])] do       
       q := a*set3[i][j]+q;  
       a := a*degree;  
     od;
     Add(set4,[]);
     for j in set1[i] do
        Add(set4[i],q/(a-1)+j/degree);
     od;
  od;



  return [degree,set4,[]];
end;



CreateIMGMachineFromPuzzleOld_Fr := function(puzzle) #wrong algorithm
# create a new img machine that is approximate puzzle.machine
# by the first return of critical points

 local i, j, v, q, k, a, b,  sets,  set1, set2, set3, set4, set,  gens0, gens1, image,  
 degree, s, inclus, bool;
 #IsRefinest_Fr(puzzle);
 #if puzzle.refinest then return false; fi;

 degree := Length(AlphabetOfFRObject(puzzle.machine)); 
 gens0 := GeneratorsOfGroup(puzzle.group0);
 gens1 := GeneratorsOfGroup(puzzle.group1);

 sets := List([1..Length(gens0)], i ->[]);
 set := [];
 inclus := List([1..Length(gens1)], i ->0);
 image := List([1..Length(gens1)], i ->0);

# searching for reducible puzzle pieces: 

 for i in [1..Length(gens0)] do     
    v := Abelianization_Fr(Image(puzzle.embedding,gens0[i]));
  
    for j in [1..Length(v)] do
       if not IsBound(v[j]) then continue; fi;
       Add(sets[i], j);
       inclus[j] := i;
    od;
 od;


 for i in [1..Length(gens1)] do     
    v := Abelianization_Fr(Image(puzzle.covering,gens1[i]));
    j := Maximum(v);
    if j > 1 then 
      Add(set,[i,j]);
    fi;
    image[i]:= Position(v,j);
 od;

#Print(inclus,"\n",image,"\n",sets,"\n\n");
  set1 :=[];
  set2 :=[]; 
  set3 :=[];
  set4 :=[];
  for i in set do 
     q := [];
     j:= Image(puzzle.embedding,gens0[image[i[1]]]);
     v := Abelianization_Fr(j);
     
     while not IsBound(v[i[1]]) do
        Add(q,v,1);
        j := Image(puzzle.embedding,Image(puzzle.covering,j));
        v := Abelianization_Fr(j);

     od;
    
     Add(set1,q);
  od;
 # Print(set1,"\n\n");

  for i in [1..Length(set1)] do
    q:=[];
    a := set[i][1];#Position(set1[i][1],Maximum(set1[i][1]));
    Add(q,a,1);
    for j in [1..Length(set1[i])] do
       bool := false;
       for k in [1..Length(set1[i][j])]do
           if not IsBound(set1[i][j][k]) then continue; fi;
           if a in sets[image[k]] then
               if bool then 
                   Print("IMGMachine is not unique");
                   continue; 
               fi;  
               bool := true;
               a := k;
               Add(q,a,1);
           fi; 

       od;

    od;
    Add(set2,q);
  od;
 

  for i in [1..Length(set1)] do
     q:=[];

     for j in [1..Length(set2[i])] do

       Add(q,
              Int( (Length(ExtRepOfObj(Image(puzzle.covering,gens1[set2[i][j]])))/2+1)/Length(gens0)/2 ),
            1);
           
     od;
     Add(set3,q);
  od;


  for i in [1..Length(set1)] do

     q:=0;
     a:=1;
     for j in [1..Length(set3[i])] do       
       q := a*set3[i][j]+q;  
       a := a*degree;  
     od;
     Add(set4,q/(a-1));
  od;



#Print(set2,"\n\n");
#Print(set3,"\n\n",degree);
  return set4;
end;


CreateIMGMachineFromPuzzleOld2_Fr := function(puzzle) #wrong algorithm
# create a new img machine that is approximate puzzle.machine
# by the first return of critical points

 local i, j, v, q, k, a, b,  sets,  set1, set2, set3, set4, set,  gens0, gens1, image, hom,  
 degree, s, inclus, bool;
 #IsRefinest_Fr(puzzle);
 #if puzzle.refinest then return false; fi;

 degree := Length(AlphabetOfFRObject(puzzle.machine)); 
 gens0 := GeneratorsOfGroup(puzzle.group0);
 gens1 := GeneratorsOfGroup(puzzle.group1);

 sets := List([1..Length(gens0)], i ->[]);
 set := [];
 inclus := List([1..Length(gens1)], i ->0);
 image := List([1..Length(gens1)], i ->0);
 hom := PullbackHomOfPuzzle_Fr(puzzle);
# searching for reducible puzzle pieces: 

 for i in [1..Length(gens0)] do     
    v := Abelianization_Fr(Image(puzzle.embedding,gens0[i]));
  
    for j in [1..Length(v)] do
       if not IsBound(v[j]) then continue; fi;
       if not Image(hom,gens1[j])=One(StateSet(puzzle.machine)) then   
          Add(sets[i], j);          
       fi;
       inclus[j] := i; 
    od;
 od;

 set1 :=[];
 

 for i in [1..Length(gens1)] do     
    v := Abelianization_Fr(Image(puzzle.covering,gens1[i]));
    j := Maximum(v);
    image[i]:= Position(v,j);
    if j > 1 then 
      Add(set,[i,j]);

      v := LoopPreimagesByFrMachine_Fr(puzzle.machine, Image(puzzle.homomorphism, gens0[image[i]]));
#Print(Int((Length(ExtRepOfObj(Image(puzzle.covering,gens1[i])))/2+1)/Length(gens0)/2));
      for k in v do
#Print(k[4],"\n");
          a := Int((Length(ExtRepOfObj(Image(puzzle.covering,gens1[i])))/2+1)/Length(gens0)/2) +1;
          if  a in k[4] 
          then
             Add(set1,k[4]-a);
          fi;
      od;
    fi;  
 od;
#Print(set1,"\n");
#Print(sets,"\nimage",image,"\ninclus",inclus,"\nset",set); return fail;
#Print(inclus,"\n",image,"\n",sets,"\n\n");
  set2 :=[]; 
  set3 :=[];
  set4 :=[];
  for i in set do 
     q := [i[1]];
     j:= image[i[1]];#Image(puzzle.embedding,gens0[image[i[1]]]);
    # v := Abelianization_Fr(j);
     
     while not i[1] in sets[j] do
        #j := inclus[j]; 
        if Length(sets[j])> 1 then
           Print("IMGMachine is not unique");
           return fail; 
        fi;
        j := sets[j][1];
        Add(q,j);
        j:= image[j];
     od;
    
     Add(set2,q);
  od;
 # Print(set1,"\n\n");

#  for i in [1..Length(set1)] do
#    q:=[];
#    a := set[i][1];#Position(set1[i][1],Maximum(set1[i][1]));
#    Add(q,a,1);
#    for j in [1..Length(set1[i])] do
#       bool := false;
#       for k in [1..Length(set1[i][j])]do
#           if not IsBound(set1[i][j][k]) then continue; fi;
#           if a in sets[image[k]] then
#               if bool then 
#                   Print("IMGMachine is not unique");
#                   continue; 
#               fi;  
#               bool := true;
#               a := k;
#               Add(q,a,1);
#           fi; 
#
#       od;
#
#    od;
#    Add(set2,q);
#  od;
 

  for i in [1..Length(set2)] do
     q:=[];

     for j in [1..Length(set2[i])] do

       Add(q,
              Int( (Length(ExtRepOfObj(Image(puzzle.covering,gens1[set2[i][j]])))/2+1)/Length(gens0)/2 ),
            1);
           
     od;
     Add(set3,q);
  od;


  for i in [1..Length(set2)] do

     q:=0;
     a:=1;
     for j in [1..Length(set3[i])] do       
       q := a*set3[i][j]+q;  
       a := a*degree;  
     od;
     Add(set4,[]);
     for j in set1[i] do
        Add(set4[i],q/(a-1)+j/degree);
     od;
  od;



#Print(set2,"\n\n");
#Print(set3,"\n\n",degree);
  return [degree,set4,[]];
end;