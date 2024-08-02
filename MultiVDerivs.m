function [dV] = MultiVDerivs(S,dS)
%MULTIVDERIVS finds the derivative of the eigenvector matrix, V, at 
%multiple times.



%Choose Dimensions.
Dim=size(S,1:3);
%Preallocate.
a=zeros(1,1,Dim(3));
b=zeros(1,1,Dim(3));
da=zeros(1,1,Dim(3));
db=zeros(1,1,Dim(3));
%Check Symmetry.
if numel(Dim)<3 || Dim(3)==1
    if issymmetric(S(:,:,1))==0
        error('S should be symmetric.')
    end
elseif issymmetric(S(:,:,2))==0 || issymmetric(S(:,:,end))==0
    error('S should be symmetric.')
end
%Choose correct formula for dV.
if isequal(Dim(1:2),[2,2])
    %Define each element of S and their derivatives.
    Sxx=S(1,1,:);
    Syy=S(2,2,:);
    Sxy=S(1,2,:);
    dSxx=dS(1,1,:);
    dSyy=dS(2,2,:);
    dSxy=dS(1,2,:);
    %Define the larger Eigenvalue, Lambda_1, and its derivative.
    Lam1=(Sxx+Syy+sqrt(4*(Sxy.^2)+(Sxx-Syy).^2))/2;
    dLam1=(dSxx+dSyy+((4*dSxy.*Sxy+dSxx.*Sxx+dSyy.*Syy-dSxx.*Syy-dSyy.*Sxx)...
        ./sqrt(4*(Sxy.^2)+(Sxx-Syy).^2)))/2;
    for i=1:Dim(3)
        if Sxx(i)>Syy(i) %Case 1. (Lam1 > or = Sxx)
            a(i)=Lam1(i)-Syy(i);
            da(i)=dLam1(i)-dSyy(i);
            b(i)=Sxy(i);
            db(i)=dSxy(i);
        elseif Syy(i)>Sxx(i) %Case 2. (Lam1 > or = Syy)
            a(i)=Sxy(i);
            da(i)=dSxy(i);
            b(i)=Lam1(i)-Sxx(i);
            db(i)=dLam1(i)-dSxx(i);
        elseif Sxx(i)==Syy(i) && Sxy(i)==0 %Degenerate Eigenvalues case.
            %Note Lambda_1 is greater than both Sxx and Syy if Sxy is non-zero,
            %so no issues of 1/0 or zero matrices can arise unless Sxy=0.
            warning('dV:equil','Missing dVdt formula for equilibrium case.')
            a(i)=nan;
            da(i)=nan;
            b(i)=nan;
            db(i)=nan;
        else
            error('Undefined error.')
        end
    end
%Calculate dV.
M=[b,a;-a,b];
dV=M.*((da.*b)-(db.*a))./((a.^2+b.^2).^(3/2));
elseif isequal(Dim(1:2),[3,3])
    error('Formula for 3D dVdt missing.')
else
    error('Incorrect format for S.')
end