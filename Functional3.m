function [G_1,G_2,G_3,f,Time]=Functional3(Input,N,Z,Functionals,Plot,OverrideType,Override)
%FUNCTIONAL3 Evaluates all the suitability of 3 chosen "Functionals" 
%to charecterise f.
%
%Inputs:
%Input - Cell array containing multiple input files. Accepts ET_R values, 
%file name strings, output matrices (M) or {Time,f} cell pairs.
%N - N value for the data, vector or scalar input.
%Z - Z value for the data, vector or scalar input.
%Functionals - Cell array containing a set of 3 functionals or strings
%containing the names of pre-defined functionals in the GLaMMFunctional
%function.
%Plot - Optional logical vairable that supresses the plot if 0.
%OverrideType - Optional parameters to set functional you wish to limit.
%Override - Vector of min and max cutoff values based on OverrideType.
%
%Outputs:
%G_1,G_2,G_3 - Cell array (for each input) containing evaluations of one 
%              functional at each timestep.
%f - Cell array (for each input) of Cell arrays indexed by i,j containing 
%    3D arrays (s,s',t) of f_{ij}(s,s').
%Time - Cell array (for each input) containing vectors detailing the
%       timesteps for each data.

%%
%Default to empty Override if undefined.
if exist('OverrideType','var')==0
    OverrideType=[];
    Override=[];
elseif exist('Override','var')==0
    Override=[];
end
%Define number of samples.
SampNo=numel(Input);
%Preallocation.
f=cell(1,SampNo);
G_1=cell(1,SampNo);
G_2=cell(1,SampNo);
G_3=cell(1,SampNo);
Time=cell(1,SampNo);
%Allow single value N and Z.
if numel(N)==1
    N=N*ones(1,SampNo);
end
if numel(Z)==1
    Z=Z*ones(1,SampNo);
end
%Evaluate each input.
for i=1:SampNo
    [Time{i},f{i},G_k]=FunctionalEval(Input{i},N(i),Z(i),...
        Functionals,OverrideType,Override);
    G_1{i}=G_k{1};
    G_2{i}=G_k{2};
    G_3{i}=G_k{3};
end
%%
%Plot.
if exist('Plot','var')==0 || Plot~=0
figure
hold on
for i=1:SampNo
    plot3(G_1{i},G_2{i},G_3{i});
end
hold off
%Label the plot.
if isa(Functionals{1},'function_handle')
    xlabel(func2str(Functionals{1}))
else
    xlabel(Functionals{1})
end
if isa(Functionals{2},'function_handle')
    ylabel(func2str(Functionals{2}))
else
    ylabel(Functionals{2})
end
if isa(Functionals{3},'function_handle')
    zlabel(func2str(Functionals{3}))
else
    zlabel(Functionals{3})
end
legend('location','northwest')
grid on
end
end

%Example useage:
%[SingleCell,ETR1]=TestData(1);
%[G1,G2,G3]=Functional3({SingleCell{3},SingleCell{10}},75,25,{'NSD','NSDCs','PrinOrienVar'});