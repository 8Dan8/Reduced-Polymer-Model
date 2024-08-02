function [Time,f,I_k]=FunctionalEval(Input,N,Z,Functionals,OverrideType,Override,OutputFile)
%FUNCTIONEVAL Evaluates all functionals within "Functionals" using the
%selected input.
%
%Inputs:
%Input - Accepts ET_R values, file name strings, output matrices (M) or 
%{Time,f} cell pairs.
%N - N value for the data.
%Z - Z value for the data.
%Functionals - cell array containing a list of functionals or strings
%containing the names of pre-defined functionals in the GLaMMFunctional
%function.
%OverrideType - Optional parameters to set functional you wish to limit.
%Override - Vector of min and max cutoff values based on OverrideType.
%OutputFile - If defined, writes the values of I_k into a file named the
%same as this string variable.
%
%Outputs:
%Time - Vector containing timesteps.
%f - Cell indexed by i,j containing 3D arrays (s,s',t) of f_{ij}(s,s').
%I_k - Matrix containing evaluations of each functional at each timestep.

%Preallocate I_k.
k=numel(Functionals);
I_k=cell(k,1);
%Default to empty Override if undefined.
if exist('OverrideType','var')==0
    OverrideType=[];
    Override=[];
else
    if exist('Override','var')==0
        Override=[];
    end
end
%Find f for input data.
if isfloat(Input) && numel(Input)==1
    %Type 1: ET_R.
    %Change an ET_R 'Input' to a string.
    Input=GLaMMD('Fpq','Extension',1,Input,N,Z,10000);
end
if iscell(Input)==0
    %Type 2: 'string' | Type 3: M.
    %Evaluate the external (Type 2) or internal (Type 3) data matrix.
    [Time,f]=FijpqExtract(Input,N,Z,OverrideType,Override);
else
    %Type 4: {Time,f}.
    %Direct inputs, no evaluation nessecary.
    Time=Input{1};
    f=Input{2};
end
%Peform each functional on output f.
for i=1:k
    %Change any strings within 'Functionals' cell array to respective functionals.
    if isa(Functionals{i},'function_handle')==0
        Functionals{i}=GLaMMFunctional(Functionals{i},Z,N);
    end
    I_k{i}=Functionals{i}(f);
end
if exist('OutputFile','var')
    %Output to file if nessecary.
    dlmwrite(OutputFile,Functionals);
    dlmwrite(OutputFile,I_k);
end
end