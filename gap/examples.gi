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

InstallGlobalFunction(PoirierExamples, function(arg)
    if arg=[1] then
        return PolynomialIMGMachine(2,[1/7],[]);
    elif arg=[2] then
        return PolynomialIMGMachine(2,[],[1/2]);
    elif arg=[3,1] then
        return PolynomialIMGMachine(2,[],[5/12]);
    elif arg=[3,2] then
        return PolynomialIMGMachine(2,[],[7/12]);
    elif arg=[4,1] then
        return PolynomialIMGMachine(3,[[3/4,1/12],[1/4,7/12]],[]);
    elif arg=[4,2] then
        return PolynomialIMGMachine(3,[[7/8,5/24],[5/8,7/24]],[]);
    elif arg=[4,3] then
        return PolynomialIMGMachine(3,[[1/8,19/24],[3/8,17/24]],[]);
    elif arg=[5] then
        return PolynomialIMGMachine(3,[[3/4,1/12],[3/8,17/24]],[]);
    elif arg=[6,1] then
        return PolynomialIMGMachine(4,[],[[1/4,3/4],[1/16,13/16],[5/16,9/16]]);
    elif arg=[6,2] then
        return PolynomialIMGMachine(4,[],[[1/4,3/4],[3/16,15/16],[7/16,11/16]]);
    elif arg=[7] then
        return PolynomialIMGMachine(5,[[0,4/5],[1/5,2/5,3/5]],[[1/5,4/5]]);
    elif arg=[9,1] then
        return PolynomialIMGMachine(3,[[0,1/3],[5/9,8/9]],[]);
    elif arg=[9,2] then
        return PolynomialIMGMachine(3,[[0,1/3]],[[5/9,8/9]]);
    fi;
end);
#############################################################################

#E examples.gi. . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
