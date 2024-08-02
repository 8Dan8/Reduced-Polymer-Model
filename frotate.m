function [f] = frotate(fhat,V)
%FROTATE Rotates f hat back to f using sorted eignvectors V.

%Choose Dimensions.
Dim=size(fhat);
%Preallocate f.
f=cell(Dim);
%Transform f. (f=V*fhat*V^T)
if isequal(Dim,[2,2])
    for i=1:2
        for k=1:2
            f{i,k}=...
                V(i,1,:).*V(k,1,:).*fhat{1,1}+...
                V(i,2,:).*V(k,2,:).*fhat{2,2}+...
                ((V(i,1,:).*V(k,2,:))+(V(i,2,:).*V(k,1,:))).*fhat{1,2};
        end
    end
elseif isequal(Dim,[3,3])
    for i=1:3
        for k=1:3
            f{i,k}=...
                V(i,1,:).*V(k,1,:).*fhat{1,1}+...
                V(i,2,:).*V(k,2,:).*fhat{2,2}+...
                V(i,3,:).*V(k,3,:).*fhat{3,3}+...
                ((V(i,1,:).*V(k,2,:))+(V(i,2,:).*V(k,1,:))).*fhat{1,2}+...
                ((V(i,1,:).*V(k,3,:))+(V(i,3,:).*V(k,1,:))).*fhat{1,3}+...
                ((V(i,2,:).*V(k,3,:))+(V(i,3,:).*V(k,2,:))).*fhat{2,3};
        end
    end
else
    error('Incorrect dimensions of f hat.')
end
end