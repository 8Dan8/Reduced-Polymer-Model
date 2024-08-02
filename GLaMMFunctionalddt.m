function [DFunctional,pqmatrix] = GLaMMFunctionalddt(Name,Z,Nin)
%GLAMMFUNCTIONALDDT attempts to find dG/dt(f,df) for a functional called 'name'.
%
%List of defined scalar functionals:
%'Stretch Ratio' - Integral of sqrt(Trf(s,s,t)) over s.
%'Stretch Ratio Counter-s' - Integral of sqrt(Trf(s,Z-s,t)) over s.
%'Sqrt Trace Diag' - Sqrt of integral of abs(Trf(s,s,t)) over s.
%'Sqrt Trace Diag Counter-s' - Sqrt of integral of abs(Trf(s,Z-s,t)) over s.
%'Log Trace C_11' - Log of Trace C_11.
%'Trace C_11' - Integral of Trf(s,s',t)*sin(ks)*sin(ks') over s and s'.
%'Total Trace Integral' - Intergral of Trf(s,s',t) over s and s'.
%'Sigma_xx' - Proportional to integral of f_xx(s,s,t) over s.
%'Sigma_yy' - Proportional to integral of f_yy(s,s,t) over s.
%'Sigma_xy' - Proportional to integral of f_xy(s,s,t) over s.
%'NSD' - Proportional to integral of f_xx(s,s,t)-f_yy(s,s,t) over s.
%'NSD Counter-s' - Prop to integral of f_xx(s,Z-s,t)-f_yy(s,Z-s,t) over s.
%'fZ2xx' - f_xx(s,s,t) at s=Z/2.
%'fZ2yy' - f_yy(s,s,t) at s=Z/2.
%
%List of defined matrix functionals:
%'S' - Normalised end to end vector, 'rr'.
%
%List of defined transformations:
%'TO' - Translates so that time=0 is at origin.
%'log()'
%'sqrt()'
%'().^2'
%
%Outputs:
%Functional - Function handle with inputs 'f' and 'df' that calculates the 
%             value of the derivative of the functional given 'f' and 'df'.
%pqmatrix - Function handle with input 'N' that returns a logical matrix
%           detailing pq pairs required to evaluate the functional.


if regexpi(Name,'ABC')
    %Product of three funcs.
    Names=regexpi(Name(4:end),'[^*]*(?=\^)','match');
    %Define exponents.
    P=cellfun(@(x) str2double(x),regexpi(Name,'(?<=\^)\d*','match'));
    %Error check.
    if numel(Names)~=3 || numel(P)~=3
        error('ABC products only accept three functionals.');
    end
    %Create product.
    [FuncA,pqA]=GLaMMFunctional(Names{1},Z);
    [FuncB,pqB]=GLaMMFunctional(Names{2},Z);
    [FuncC,pqC]=GLaMMFunctional(Names{3},Z);
    dFuncA=GLaMMFunctionalddt(Names{1},Z);
    dFuncB=GLaMMFunctionalddt(Names{2},Z);
    dFuncC=GLaMMFunctionalddt(Names{3},Z);
    DFunctional=@(f,df) ABCfuncddt(f,df,P,FuncA,FuncB,FuncC,dFuncA,dFuncB,dFuncC);
    pqmatrix=@(N) pqA(N)+pqB(N)+pqC(N);
    FName=[Names{1},'^',P(1),'*',Names{2},'^',P(2),'*',Names{3},'^',P(3)];
elseif numel(regexpi(Name,'\*'))
    %Product of two funcs.
    Names=regexpi(Name,'[^*]*','match');
    %Error check.
    if numel(Names)~=2
        error('GLaMMFunctional only accepts a product between two functionals');
    end
    %Create product.
    [FuncA,pqA]=GLaMMFunctional(Names{1},Z);
    [FuncB,pqB]=GLaMMFunctional(Names{2},Z);
    dFuncA=GLaMMFunctionalddt(Names{1},Z);
    dFuncB=GLaMMFunctionalddt(Names{2},Z);
    DFunctional=@(f,df) (dFuncA(f,df).*FuncB(f))+(FuncA(f).*dFuncB(f,df));
    pqmatrix=@(N) pqA(N)+pqB(N);
    FName=[Names{1},'*',Names{2}];
elseif regexpi(Name,'Stretch\s*Ratio')
    %Stretch ratio.
    DFunctional=@(f,df) StretchRatioddt(f,df);
    pqmatrix=@(N) eye(N);
    FName='Stretch Ratio';
elseif regexpi(Name,'Sqrt\s*Tr\w*\s*Diag\w*\s*C.*s')
    %Counter diagonal "Sqrt Trace Diagonal".
    DFunctional=@(f,df) SqrtTraceDiagonalCounterddt(f,df);
    pqmatrix=@(N) Counterdiag(N);
    FName='Sqrt Trace Diag Counter-s';
elseif regexpi(Name,'Pri\w*\s*Ori\w*\s*Var\w*')
    %Principle orientation variance.
    warning('Principle orientation variance not defined in GLaMMFunctionalddt.')
    DFunctional=[];
    pqmatrix=@(N) eye(N);
    FName=[];
elseif regexpi(Name,'Pri\w*\s*Sin\w*\s*Sum\w*')
    %Principle sin^2 sum.
    DFunctional=@(fhat,dfhat) PrinSin2Sumddt(fhat,dfhat,Z);
    pqmatrix=@(N) eye(N);
    FName='Principle Sin2 Sum';
elseif regexpi(Name,'fxy\s*Var\w*')
    %fxy variance.
    warning('fxy variance not defined in GLaMMFunctionalddt.')
    DFunctional=[];
    pqmatrix=@(N) eye(N);
    FName=[];
elseif regexpi(Name,'Sig\w*\s*xx')
    %Sigma_xx.
    DFunctional=@(f,df) Sigmaxxddt(df);
    pqmatrix=@(N) eye(N);
    FName='Sigma_xx';
elseif regexpi(Name,'Sig\w*\s*yy')
    %Sigma_yy.
    DFunctional=@(f,df) Sigmayyddt(df);
    pqmatrix=@(N) eye(N);
    FName='Sigma_yy';
elseif regexpi(Name,'N\w*\s*S\w*\s*D\w*\s*C.*s.*Z')
    %Counter diagonal "NSD".
    DFunctional=@(f,df) CounterNSDddt(df)./Z;
    pqmatrix=@(N) Counterdiag(N);
    FName='NSD Counter-s';
elseif regexpi(Name,'N\w*\s*S\w*\s*D.*Z')
    %1st normal stress difference.
    DFunctional=@(f,df) (Sigmaxxddt(df)-Sigmayyddt(df))./Z;
    pqmatrix=@(N) eye(N);
    FName='NSD';
elseif regexpi(Name,'N\w*\s*S\w*\s*D\w*\s*C.*s')
    %Counter diagonal "NSD".
    DFunctional=@(f,df) CounterNSDddt(df);
    pqmatrix=@(N) Counterdiag(N);
    FName='NSD Counter-s';
elseif regexpi(Name,'N\w*\s*S\w*\s*D')
    %1st normal stress difference.
    DFunctional=@(f,df) Sigmaxxddt(df)-Sigmayyddt(df);
    pqmatrix=@(N) eye(N);
    FName='NSD';
elseif ~isempty(regexpi(Name,'S','once')) || ~isempty(regexpi(Name,'uu','once'))
    %
    DFunctional=@(f,df) UUddt(f,df);
    pqmatrix=@(N) eye(N);
    FName='S';
elseif regexpi(Name,'Zero')
    %Debug functional
    DFunctional=@(f,df) zeros(1,size(f{1},3));
    pqmatrix=@(N) zeros(N);
    FName='Zero';
else
    %Check if GLaMMFunctional recognises functional name.
    GLaMMFunctional(Name,Z);
    %Error catching for unset d/dt.
    error('Functional found, but d/dt of functional undefined.')
end
%Transformations.
if ~isempty(regexpi(Name,'\(.*\(','once')) || ~isempty(regexpi(Name,'\(\*\)','once'))
    %Error catching.
    error('GLaMMFunctionalddt is not yet able to perform more than one transformation.')
end
if isempty(regexpi(Name,'\*','once'))
    if regexpi(Name,'TO')
        %Translate to origin needs to be considered in Fname.
        FName=[FName,'TO'];
    end
    if regexpi(Name,'\(.+\)')
        %Define regular functional if required.
        if exist('Nin','var')
            Functional=GLaMMFunctional(FName,Z,Nin);
        else
            Functional=GLaMMFunctional(FName,Z);
        end
    end
    %Note that origin translation does not change dG/dt, only G.
    if regexpi(Name,'log\(.+\)')
        %Log transformation.
        DFunctional=@(f,df) DFunctional(f,df)./Functional(f);
    end
    if regexpi(Name,'sqrt\(.+\)')
        %Square root transformation.
        DFunctional=@(f,df) DFunctional(f,df)./(2*sqrt(Functional(f)));
    end
    if regexpi(Name,'\(.+\).*2')
        %Square transformation.
        DFunctional=@(f,df) DFunctional(f,df).*(2*Functional(f));
    end
end
end

%%

function sol=ABCfuncddt(f,df,P,FuncA,FuncB,FuncC,dFuncA,dFuncB,dFuncC)
sol=((P(1)*dFuncA(f,df).*FuncB(f).*FuncC(f))+...
    (P(2)*dFuncB(f,df).*FuncA(f).*FuncC(f))+...
    (P(3)*dFuncC(f,df).*FuncA(f).*FuncB(f))).*...
    (FuncA(f).^(P(1)-1)).*(FuncB(f).^(P(2)-1)).*(FuncC(f).^(P(3)-1));
end

%d/dt of functionals.
function lambda=StretchRatioddt(f,df)
%Derivative of Stretch Ratio functional.

%Calculate N.
N=size(f{1},1);
%Predefine integrand.
%(Tr(f(0,0))=1, Tr(df(0,0)/dt)=0, thus intergrand is 0 at s=0).
integrand=zeros(N+1,size(f{1},3));
%Integrand is tr(df(s,s))/sqrt(tr(f(s,s))).
for s=1:N
    integrand(s+1,:)=(df{1,1}(s,s,:)+2*df{2,2}(s,s,:))./...
        sqrt(f{1,1}(s,s,:)+2*f{2,2}(s,s,:));
end
%Functional is 0.5*integral from 0 to Z.
lambda=0.5*(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function STDC=SqrtTraceDiagonalCounterddt(f,df)
%Derivative of Counter Sqrt Trace Diagonal functional.

%Calculate N.
N=size(f{1},1);
%Predefine integrands.
%(Tr(f(0,Z))=0, Tr(df(0,Z)/dt)=0, thus both intergrands are 0 at s=0).
integrandT=zeros(N+1,size(f{1},3));
integrandB=zeros(N+1,size(f{1},3));
%IntegrandT is tr(df(s,Z-s)), IntegrandB is tr(f(s,Z-s)).
for s=1:(N-1)
    integrandT(s+1,:)=df{1,1}(s,(N-s),:)+2*df{2,2}(s,(N-s),:);
    integrandB(s+1,:)=f{1,1}(s,(N-s),:)+2*f{2,2}(s,(N-s),:);
end
%Functional is 0.5*integralT/sqrt(integralB) both from 0 to Z.
STDC=0.5*((sum(integrandT,1)-0.5*(integrandT(1,:)+integrandT(N+1,:)))/N)./...
    sqrt((sum(integrandB,1)-0.5*(integrandB(1,:)+integrandB(N+1,:)))/N);
end

function Sum=PrinSin2Sumddt(fhat,dfhat,Z)
%Derivative of sum sin^2(angle) to principle direction.
%(Assumes f_11 is principle direction.)

%Set Dimensions.
Dim=2;
%Calculate N.
N=size(fhat{1},1);
%Calculate Tlen.
Tlen=size(fhat{1},3);
%Calculate points to ignore given proximity to equilbrium 'edge' cases.
svec=(1+ceil(N/(2*Z))):(N+1-ceil(N/(2*Z)));
%Preallocate. (Form is (a,b,t,s) for s=s'.)
fss=zeros(Dim,Dim,Tlen,N+1);
dfss=zeros(Dim,Dim,Tlen,N+1);
Vs=zeros(Dim,Dim,Tlen,N+1);
dVs=zeros(Dim,Dim,Tlen,N+1);
%Predefine fss,dfss at feq, Vs=I, dVs=0. (for edges only).
for i=1:Dim
    fss(i,i,:,:)=1/3;
    Vs(i,i,:,:)=1;
    %dfss(i,i,:,:)=0;
    %dVs(i,i,:,:)=0;
end
%Create summation terms.
for i=1:Dim
    %Extract from f.
    for k=i:Dim
        for s=1:N
            fss(i,k,:,s+1)=fhat{i,k}(s,s,:);
            dfss(i,k,:,s+1)=dfhat{i,k}(s,s,:);
        end
    end
    %Symmetry speedup.
    for j=i:Dim
        fss(j,i,:,:)=fss(i,j,:,:);
        dfss(j,i,:,:)=dfss(i,j,:,:);
    end
end
%Calculate S,dS at s.
Trfss=fss(1,1,:,:)+fss(2,2,:,:)+fss(Dim,Dim,:,:);
dTrfss=dfss(1,1,:,:)+dfss(2,2,:,:)+dfss(Dim,Dim,:,:);
Ss=fss./Trfss;
dSs=(dfss./Trfss)-(fss.*(dTrfss./(Trfss.^2)));
%Find the value of the first element in the first normalised eigenvector.
%(V11=cosTh)
for s=svec
    Vs(:,:,:,s)=Seigen(Ss(:,:,:,s));
    dVs(:,:,:,s)=MultiVDerivs(Ss(:,:,:,s),dSs(:,:,:,s));
end
CosTh=Vs(1,1,:,:);
dCosTh=dVs(1,1,:,:);
%Calculate derivative of sum of SinTh^2.
%This may or may not be var(SinTh), not sure if sum(SinTh(s))=0 is always
%true.
Sum=permute(sum(-2*CosTh.*dCosTh,4)/(numel(svec)),[1,3,2]);
end

function xx=Sigmaxxddt(df)
%Derivative of Sigma_xx.

%Calculate N.
N=size(df{1},1);
%Predefine integrand.
%(df_xx(0,0)/dt=0, thus intergrand is 0 at s=0).
integrand=zeros(N+1,size(df{1},3));
for s=1:N
    integrand(s+1,:)=df{1,1}(s,s,:);
end
Ge=1;
xx=((12*Ge)/5)*(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function yy=Sigmayyddt(df)
%Derivative of Sigma_yy.

%Calculate N.
N=size(df{1},1);
%Predefine integrand.
%(df_yy(0,0)/dt=0, thus intergrand is 0 at s=0).
integrand=zeros(N+1,size(df{1},3));
for s=1:N
    integrand(s+1,:)=df{2,2}(s,s,:);
end
Ge=1;
yy=((12*Ge)/5)*(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function CNSD=CounterNSDddt(df)
%Derivative of Counter Diagonal NSD.

%Calculate N.
N=size(df{1},1);
%Predefine integrand.
%(df_xx(0,Z)/dt=0, df_yy(0,Z)/dt=0, thus intergrand is 0 at s=0).
integrand=zeros(N+1,size(df{1},3));
for s=1:(N-1)
    integrand(s+1,:)=df{1,1}(s,N-s,:)-df{2,2}(s,N-s,:);
end
Ge=1;
CNSD=((12*Ge)/5)*(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function uuddt=UUddt(f,df)
%Derivative of normalised end to end functional.

%Set Dimensions.
Dim=2;
%Calculate N.
N=size(f{1},1);
%Calculate Tlen.
Tlen=size(f{1},3);
%Preallocate.
integrand=zeros(Dim,Dim,Tlen,N+1);
dintegrand=zeros(Dim,Dim,Tlen,N+1);
Trintf=zeros(1,1,Tlen);
Trintdf=zeros(1,1,Tlen);
%Predefine integrand as feq (for edges only).
for i=1:Dim
    integrand(i,i,:,:)=1/3;
    %dintergrand(i,i,:,:)=0;
end
%Create integrands
for i=1:Dim
    %Extract from f.
    for k=i:Dim
        for s=1:N
            integrand(i,k,:,s+1)=f{i,k}(s,s,:);
            dintegrand(i,k,:,s+1)=df{i,k}(s,s,:);
        end
    end
    %Symmetry speedup.
    for j=i:Dim
        integrand(j,i,:,:)=integrand(i,j,:,:);
        dintegrand(j,i,:,:)=dintegrand(i,j,:,:);
    end
end
%Calculate integrals with trapezium rule.
intf=(sum(integrand,4)-0.5*(integrand(:,:,:,1)+integrand(:,:,:,N+1)))/N;
intdf=(sum(dintegrand,4)-0.5*(dintegrand(:,:,:,1)+dintegrand(:,:,:,N+1)))/N;
if Dim==2
    for t=1:Tlen
        %Calculate trace. (ASSUMES yy=zz!!)
        Trintf(1,1,t)=trace(intf(:,:,t))+intf(2,2,t);
        Trintdf(1,1,t)=trace(intdf(:,:,t))+intdf(2,2,t);
    end
elseif Dim==3
    for t=1:Tlen
        %Calculate trace.
        Trintf(1,1,t)=trace(intf(:,:,t));
        Trintdf(1,1,t)=trace(intdf(:,:,t));
    end
end
%Use duu/dt formula.
uuddt=(intdf./Trintf)-(intf.*(Trintdf./(Trintf.^2)));
end

%%
function pqmatrix=Counterdiag(N)
%Logical matrix for q=Z-p. (note the exclusion of (0,N)=0 and (N,0)=0)
pqmatrix=zeros(N);
pqmatrix(1:N-1,1:N-1)=fliplr(eye(N-1));
end