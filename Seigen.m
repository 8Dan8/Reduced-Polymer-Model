function [V,Val] = Seigen(S)
%FPRINCIPLE Finds the matrix of sorted eigenvectors, V, from the principle
%rotation, S.

%Preallocate. (use S for size only)
V=S;
for t=1:size(S,3)
    %Calculate eigen values and vectors of S.
    %(eig ensures orthonormal for symmetric input.)
    [Vecu,Valu]=eig(S(:,:,t));
    %Find principle direction. (Principle eigenvalue)
    [Val,ind]=sort(diag(Valu),'descend');
    %Create transformation matrix (Reorder eigenvectors).
    V(:,:,t)=Vecu(:,ind);
    %Ensure element V_11 is +ve.
    if V(1,1,t)~=0
        V(:,:,t)=sign(V(1,1,t)).*V(:,:,t);
    end
end
end