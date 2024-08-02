function [Sigma,Sigdim] = Rolie2D(Time,Kappa,Z,Te,beta)
%ROLIE2D Solves the Rolie-Poly model using ODE45.

%Correct if inputs not in cell format.
if iscell(Kappa)==0
    Kappa={Kappa};
end
if iscell(Time)==0
    Time={Time};
end
%Preallocate.
Sigdim=cell(1,numel(Kappa));
Sigma=cell(1,numel(Kappa));
%Setup inital time and sigma values [Sxx,Syy,Sxy].
Init=0;
Inis=[1;1;0];
%Set constants.
Td=3*Te*Z^3;
Tdf=Td*(1-(3.38/sqrt(Z))+(4.17/Z)+(-1.55/(Z^(3/2))))/2;
TR=Te*Z^2;
Gf=1-(1.69/sqrt(Z))+(2/Z)+(-1.24/(Z^(3/2)));
for i=1:numel(Kappa)
    %Adjust time vector. (for multi-flow)
    TimeShift=Time{i}-Init;
    %Solve Rolie-Poly ODE.
    soln=ode45(@(t,s) RolieODE2D(s,Kappa{i},Tdf,TR,beta),[0,TimeShift(end)],Inis);
    %Record outputs.
    Sigdim{i}=deval(soln,TimeShift);
    Sigma{i}=0.8*Gf.*Sigdim{i};
    %Initialise next flow. (for multi-flow)
    Init=Time{i}(end);
    Inis=Sigdim{i}(:,numel(TimeShift));
end
end