function [fhat,S,V] = fPrinciple(f)
%FPRINCIPLE finds the value of f in the principle direction, fhat.

%Setup.
SFunc=GLaMMFunctional('S');
%Find S.
S=SFunc(f);
%Find eigenmatrix V.
V=Seigen(S);
%Remove orentation from f to find fhat.
fhat=frotate(f,permute(V,[2,1,3]));
end