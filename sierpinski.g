FindCantor := function(M,multicurve,limit)
    # search for a cantor multicurve, i.e. a curve lifting to a multiple of itself
    local len, w, x, mat, row, i, j, c, d, group, pi, gens, peripheral;
    
    len := Length(multicurve);
    gens := GeneratorsOfGroup(StateSet(M));
    group := FreeGroup(Length(gens)-1);
    c := IMGRelator(M);
    pi := GroupHomomorphismByImagesNC(StateSet(M),group,List([1..Length(gens)],i->Subword(c,i,i)),Concatenation(GeneratorsOfGroup(group),[Product(List(Reversed(GeneratorsOfGroup(group)),Inverse))]));

    w := PUSHRECURSION@fr(pi,M);

    peripheral := List(GeneratorsOfSemigroup(StateSet(M)),x->CyclicallyReducedWord(x^pi));
    
    multicurve := List(multicurve,x->CyclicallyReducedWord(x^pi));
    mat := [];
    for i in multicurve do
        d := w(i);
        row := List([1..len],i->0);
        for i in Cycles(PermList(d[2]),AlphabetOfFRObject(M)) do
            c := CyclicallyReducedWord(Product(d[1]{i}));
#            Error(c,"\n");
            if ForAny(peripheral,x->IsConjugate(group,x,c)) then
                continue; # peripheral curve
            fi;
            j := First([1..len],j->IsConjugate(group,c,multicurve[j])
                       or IsConjugate(group,c^-1,multicurve[j]));
            if j<>fail then
                row[j] := row[j] + 1;
            elif len<limit then # add one more curve
                for j in mat do Add(j,0); od;
                Add(row,1);
                len := len+1;
                Add(multicurve,c);
            fi;
        od;
        Add(mat,row);
    od;
    Info(InfoFR,1,"Thurston matrix is ",mat);

    x := List(EquivalenceClasses(StronglyConnectedComponents(BinaryRelationOnPoints(List([1..len],x->Filtered([1..len],y->IsPosRat(mat[x][y])))))),Elements);
    for i in x do
        if PERRONMATRIX@fr(mat{i}{i}) then # there's an eigenvalue >= 1
            d := rec(machine := M,
                     obstruction := [],
                     matrix := mat{i}{i});
            for j in i do
                c := [PreImagesRepresentative(pi,multicurve[j])];
                Append(d.obstruction,c);
            od;
            return d;
        fi;
    od;
    return fail;
end;
