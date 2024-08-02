function [dG,df]=dGdt(Functional,f,Time,N,Z,FlowType,StepET_R,cnu)

%Find number of functionals.
FuncNo=size(Functional,2);
%Preallocate.
DFunctional=cell(1,FuncNo);
pqpart=cell(1,FuncNo);
dG=cell(1,FuncNo);
pqmatrix=@(N) 0;
for i=1:FuncNo
    %Convert functional strings to function handles if nessecary.
    if isa(Functional{1,i},'function_handle')==0
        %Functional name input.
        [DFunctional{i},pqpart{i}]=GLaMMFunctionalddt(Functional{1,i},Z,N);
    elseif size(Functional,1)>1
        %FunctionalSet input (from Linear Combination).
        DFunctional{i}=Functional{(end-1),i};
        pqpart{i}=Functional{end,i};
    else
        %Significant speedup if pqmatrix is defined.
        warning('Using default pqmatrix, expect long calculation times.')
        error('Differental of Functional undefined.')
        %DFunctional{i}=Functional{i};
        %pqpart{i}=@(N) ones(N);
    end
    %Sum to find which pq elements need to be evaluated.
    pqmatrix=@(N) pqmatrix(N)+pqpart{i}(N);
end
%Find df.
df=Multifderivs(f,Time,N,Z,FlowType,StepET_R,cnu,pqmatrix(N));
for i=1:FuncNo
    %Apply functionals to df.
    dG{i}=DFunctional{i}(f,df);
end
end