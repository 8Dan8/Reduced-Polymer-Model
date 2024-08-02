function f=RolieODE2D(s,kappa,Td,TR,beta)
%ROLIEODE2D is a function of the Rolie-Poly model for use with ODE45.

delta=-0.5;
%Read input [Sigxx,Sigyy,Sigxy].
sigma=[s(1),s(3);s(3),s(2)];
%Calculate ODE.
ks=kappa*sigma;
sk=sigma*(kappa');
Trsigma=s(1)+2*s(2);
dsigma=ks+sk-((sigma-eye(2))/Td)-(2*(1-sqrt(3/Trsigma))/TR)*...
    (sigma+(beta*((Trsigma/3)^delta)*(sigma-eye(2))));
%Ouput results.
f=zeros(3,1);
f(1)=dsigma(1,1);
f(2)=dsigma(2,2);
f(3)=dsigma(1,2);