function [Functional,pqmatrix] = GLaMMFunctional(Name,Z,Wij,Nin)
%GLAMMFUNCTIONAL attempts to find a functional called 'name'.
%
%The variable 'Name' should be a string that is one of the following,
%possibly with a single tranformation.
%List of defined scalar functionals:
%'Stretch Ratio' - Integral of sqrt(Trf(s,s,t)) over s.
%'Stretch Ratio Counter-s' - Integral of sqrt(Trf(s,Z-s,t)) over s.
%'Sqrt Trace Diag' - Sqrt of integral of abs(Trf(s,s,t)) over s.
%'Sqrt Trace Diag Counter-s' - Sqrt of integral of abs(Trf(s,Z-s,t)) over s.
%'Principle Orientation Variance' - Variance in direction about f_11.
%'Principle Sine2 Sum" - Alternate measure of variation in direction over s.
%'fxy Variance' - fxy variance about mean 0.
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
%'NSDWij' - Sum of Wij(f_xxij(t)-f_yyij(t)), requires Wij to be defined.
%
%List of defined matrix functionals:
%'S' - Normalised end to end vector, 'rr'.
%
%List of defined transformations:
%'TO' - Translates so that time=0 is at origin. (Speedup if Nin defined.)
%'log()'
%'sqrt()'
%'().^2'
%
%Outputs:
%Functional - Function handle with input 'f' that calculates the value of
%             the functional at 'f'.
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
    %Define Wij if unset.
    if exist('Wij','var')
        Wij=[];
    end
    %Create product.
    [FuncA,pqA]=GLaMMFunctional(Names{1},Z,Wij);
    [FuncB,pqB]=GLaMMFunctional(Names{2},Z,Wij);
    [FuncC,pqC]=GLaMMFunctional(Names{3},Z,Wij);
    Functional=@(f) ABCfunc(f,P,FuncA,FuncB,FuncC);
    pqmatrix=@(N) pqA(N)+pqB(N)+pqC(N);
elseif numel(regexpi(Name,'\*'))
    %Product of two funcs.
    Names=regexpi(Name,'[^*]*','match');
    %Error check.
    if numel(Names)~=2
        error('GLaMMFunctional only accepts a product between two functionals');
    end
    %Define Wij if unset.
    if exist('Wij','var')
        Wij=[];
    end
    %Create product.
    [FuncA,pqA]=GLaMMFunctional(Names{1},Z,Wij);
    [FuncB,pqB]=GLaMMFunctional(Names{2},Z,Wij);
    Functional=@(f) FuncA(f).*FuncB(f);
    pqmatrix=@(N) pqA(N)+pqB(N);
elseif regexpi(Name,'NSDWij')
    %General functional for xx-yy.
    Functional=@(f) NSDWij(f,Wij);
    pqmatrix=@(N) ones(N);
elseif regexpi(Name,'Stretch\s*Ratio\s*C.*s')
    %Counter diagonal "Stretch Ratio".
    Functional=@(f) StretchRatioCounter(f);
    pqmatrix=@(N) Counterdiag(N);
elseif regexpi(Name,'Stretch\s*Ratio')
    %Stretch ratio.
    Functional=@(f) StretchRatio(f);
    pqmatrix=@(N) eye(N);
elseif regexpi(Name,'Sqrt\s*Tr\w*\s*Diag\w*\s*C.*s')
    %Counter diagonal "Sqrt Trace Diagonal".
    Functional=@(f) SqrtTraceDiagonalCounter(f);
    pqmatrix=@(N) Counterdiag(N);
elseif regexpi(Name,'Sqrt\s*Tr\w*\s*Diag\w*')
    %Sqrt Trace Diagonal.
    Functional=@(f) SqrtTraceDiagonal(f);
    pqmatrix=@(N) eye(N);
elseif regexpi(Name,'Pri\w*\s*Ori\w*\s*Var\w*')
    %Principle orientation variance.
    Functional=@(fhat) PrinOrienVar(fhat,Z);
    pqmatrix=@(N) eye(N);
elseif regexpi(Name,'Pri\w*\s*Sin\w*\s*Sum\w*')
    %Principle sin^2 sum.
    Functional=@(fhat) PrinSin2Sum(fhat,Z);
    pqmatrix=@(N) eye(N);
elseif regexpi(Name,'fxy\s*Var\w*')
    %fxy variance.
    Functional=@(f) fxyVar(f);
    pqmatrix=@(N) eye(N);
elseif regexpi(Name,'l\w*Tr\w*\s*C\w*11')
    %Log trace C11.
    Functional=@(f) log(TraceC11(f,Z));
    pqmatrix=@(N) ones(N);
elseif regexpi(Name,'Tr\w*\s*C\w*11')
    %Trace C11.
    Functional=@(f) TraceC11(f,Z);
    pqmatrix=@(N) ones(N);
elseif regexpi(Name,'T\w*\s*T\w*\s*I\w*')
    %
    Functional=@(f) TotalIntegral(f,Z);
    pqmatrix=@(N) ones(N);
elseif regexpi(Name,'Sig\w*\s*xx')
    %Sigma_xx.
    Functional=@(f) Sigmaxx(f);
    pqmatrix=@(N) eye(N);
elseif regexpi(Name,'Sig\w*\s*yy')
    %Sigma_yy.
    Functional=@(f) Sigmayy(f);
    pqmatrix=@(N) eye(N);
elseif regexpi(Name,'Sig\w*\s*xy')
    %Sigma_xy.
    Functional=@(f) Sigmaxy(f);
    pqmatrix=@(N) eye(N);
elseif regexpi(Name,'N\w*\s*S\w*\s*D\w*\s*C.*s.*Z')
    %Counter diagonal "NSD".
    Functional=@(f) CounterNSD(f)./Z;
    pqmatrix=@(N) Counterdiag(N);
elseif regexpi(Name,'N\w*\s*S\w*\s*D.*Z')
    %1st normal stress difference.
    Functional=@(f) (Sigmaxx(f)-Sigmayy(f))./Z;
    pqmatrix=@(N) eye(N);
elseif regexpi(Name,'N\w*\s*S\w*\s*D\w*\s*P\w*\s*C.*s')
    %Counter diagonal "NSD". (With phantom)
    Functional=@(f) CounterpNSD(f);
    pqmatrix=@(N) Counterdiag(N);
elseif regexpi(Name,'N\w*\s*S\w*\s*D\w*\s*C.*s')
    %Counter diagonal "NSD".
    Functional=@(f) CounterNSD(f);
    pqmatrix=@(N) Counterdiag(N);
elseif regexpi(Name,'N\w*\s*S\w*\s*D')
    %1st normal stress difference.
    Functional=@(f) Sigmaxx(f)-Sigmayy(f);
    pqmatrix=@(N) eye(N);
elseif regexpi(Name,'fZ2xx')
    %
    Functional=@(f) permute(f{1,1}(ceil(end/2),ceil(end/2),:)-1,[1,3,2]);
    pqmatrix=@(N) z2(N);
elseif regexpi(Name,'fZ2yy')
    %
    Functional=@(f) permute(f{2,2}(ceil(end/2),ceil(end/2),:)-1,[1,3,2]);
    pqmatrix=@(N) z2(N);
elseif regexpi(Name,'S')
    %
    Functional=@(f) UU(f);
    pqmatrix=@(N) eye(N);
elseif regexpi(Name,'Zero')
    %Debug functional
    Functional=@(f) zeros(1,size(f{1},3));
    pqmatrix=@(N) zeros(N);
else
    %Error catching.
    error('Unrecognised functional.')
end
%Transformations.
if ~isempty(regexpi(Name,'\(.*\(','once')) || ~isempty(regexpi(Name,'\(.*\*.*\)','once'))
    %Error catching.
    error('GLaMMFunctional is not yet able to perform more than one transformation.')
    %This would require a method to determine the order of the transfomations.
end
if isempty(regexpi(Name,'\*','once'))
    if regexp(Name,'TO')
        %Origin translaton.
        if exist('Nin','var')
            Functional=@(f) Functional(f)-Functional(feqfunc(Nin,Z));
        else
            warning('Calculating origin translation without N, may result in slow down');
            Functional=@(f) Functional(f)-Functional(feqfunc(size(f{1},1),Z));
        end
    end
    if regexpi(Name,'log\(.+\)')
        %Log transformation.
        Functional=@(f) log(Functional(f));
    end
    if regexpi(Name,'sqrt\(.+\)')
        %Square root transformation.
        Functional=@(f) sqrt(Functional(f));
    end
    if regexpi(Name,'\(.+\).*2')
        %Square transformation.
        Functional=@(f) (Functional(f)).^2;
    end
end
end

%%
%Functionals

function sol=ABCfunc(f,P,FuncA,FuncB,FuncC)
sol=(FuncA(f).^P(1)).*(FuncB(f).^P(2)).*(FuncC(f).^P(3));
end

function Gamma=NSDWij(f,Wij)
%General functional for xx-yy.

%Find N/Z.
NdZ=size(f{1},1)/size(Wij,1);
%Group points for speed if required.
if mod(NdZ,1)~=0
    error('Dimensions of Wij are incorrect.')
end
if NdZ~=1
    Wij=repelem(Wij,NdZ,NdZ);
end
%Calculate N.
Gamma=permute(sum(Wij*(f{1,1}-f{2,2}),1,2),[2,3,1]);
end

function lambda=StretchRatio(f)
%Stretch Ratio functional.

%Calculate N.
N=size(f{1},1);
%Predefine integrand as feq (for edges only).
integrand=ones(N+1,size(f{1},3));
for s=1:N
    integrand(s+1,:)=sqrt(f{1,1}(s,s,:)+2*f{2,2}(s,s,:));
end
lambda=(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function lambda=StretchRatioCounter(f)
%Counter Diagonal Stretch Ratio functional.

%Calculate N.
N=size(f{1},1);
%Predefine integrand as feq (Non-zero elements in the middle will be overwritten).
integrand=zeros(N+1,size(f{1},3));
for s=1:(N-1)
    integrand(s+1,:)=sqrt(f{1,1}(s,(N-s),:)+2*f{2,2}(s,(N-s),:));
end
lambda=(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function STD=SqrtTraceDiagonal(f)
%Sqrt Trace Diagonal functional.

%Calculate N.
N=size(f{1},1);
%Predefine integrand as feq (for edges only).
integrand=ones(N+1,size(f{1},3));
for s=1:N
    integrand(s+1,:)=f{1,1}(s,s,:)+2*f{2,2}(s,s,:);
end
STD=sqrt((sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N);
end

function STDC=SqrtTraceDiagonalCounter(f)
%Counter Sqrt Trace Diagonal functional.

%Calculate N.
N=size(f{1},1);
%Predefine integrand as feq (Non-zero elements in the middle will be overwritten).
integrand=zeros(N+1,size(f{1},3));
for s=1:(N-1)
    integrand(s+1,:)=f{1,1}(s,(N-s),:)+2*f{2,2}(s,(N-s),:);
end
STDC=sqrt((sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N);
end

function Var=PrinOrienVar(fhat,Z)
%Variance in angle to principle direction.
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
Vs=ones(Dim,Dim,Tlen,N+1);
%Predefine fss as feq (for edges only).
for i=1:Dim
    fss(i,i,:,:)=1/3;
end
%Create integrands
for i=1:Dim
    %Extract from f.
    for k=i:Dim
        for s=1:N
            fss(i,k,:,s+1)=fhat{i,k}(s,s,:);
        end
    end
    %Symmetry speedup.
    for j=i:Dim
        fss(j,i,:,:)=fss(i,j,:,:);
    end
end
%Calculate S at s.
Trfss=fss(1,1,:,:)+fss(2,2,:,:)+fss(Dim,Dim,:,:);
Ss=fss./Trfss;
%Find the value of the first element in the first eigenvector. (V11=cosTh)
for s=svec
    Vs(:,:,:,s)=Seigen(Ss(:,:,:,s));
end
CosTh=Vs(1,1,:,:);
%Manually calculate var. (As we have data for |X| not X so E[X] is wrong.)
%Note this is not sample var, as we already have N+1 points.
Var=permute(sum((acos(CosTh)).^2,4)/(numel(svec)),[1,3,2]);
end

function Sum=PrinSin2Sum(fhat,Z)
%Sum in sin^2(angle) to principle direction.
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
Vs=ones(Dim,Dim,Tlen,N+1);
%Predefine fss as feq (for edges only).
for i=1:Dim
    fss(i,i,:,:)=1/3;
end
%Create integrands
for i=1:Dim
    %Extract from f.
    for k=i:Dim
        for s=1:N
            fss(i,k,:,s+1)=fhat{i,k}(s,s,:);
        end
    end
    %Symmetry speedup.
    for j=i:Dim
        fss(j,i,:,:)=fss(i,j,:,:);
    end
end
%Calculate S at s.
Trfss=fss(1,1,:,:)+fss(2,2,:,:)+fss(Dim,Dim,:,:);
Ss=fss./Trfss;
%Find the value of the first element in the first normalised eigenvector.
for s=svec
    Vs(:,:,:,s)=Seigen(Ss(:,:,:,s));
end
%(V11=cosTh)
CosTh=Vs(1,1,:,:);
%Calculate sum of SinTh^2.
%This may or may not be var(SinTh), not sure if sum(SinTh(s))=0 is always
%true.
Sum=permute(sum(1-(CosTh.^2),4)/(numel(svec)),[1,3,2]);
end

function Var=fxyVar(f)
%fxy variance functional.

%Calculate N.
N=size(f{1},1);
%Predefine vector.
fxyss=zeros(N+1,size(f{1},3));
%Extract fxy.
fxyss(1,:)=f{1,2}(N,N,:);
for s=1:N
    fxyss(s+1,:)=f{1,2}(s,s,:);
end
%Manually calculate var about mean 0.
Var=sum(fxyss.^2,1)/(N+1);
end

function TrC11=TraceC11(f,Z)
%Trace C_11 functional.

%Calculate N and time length.
[N,~,Tlength]=size(f{1});
%Preallocation.
feq=zeros(N);
C11xx=zeros(1,Tlength);
C11yy=zeros(1,Tlength);
%Create feq.
for p=1:N
    for q=1:N
        feq(p,q)=((abs(p-q)*Z/N)<0.5)/3;
    end
end
%Calculate elements of C_11.
for i=1:Tlength
    C11xx(i)=((sin(pi*(1:N)./N))*(f{1}(:,:,i)-feq)*(sin(pi*(1:N)./N).'));
    C11yy(i)=((sin(pi*(1:N)./N))*(f{2}(:,:,i)-feq)*(sin(pi*(1:N)./N).'));
end
%Calculate Trace C_11.
TrC11=C11xx+(2.*C11yy);
end

function lambda=TotalIntegral(f,Z)
%Total trace integral functional.

%Calculate N and time length.
[N,~,Tlength]=size(f{1});
%Preallocation.
feq=zeros(N+1);
%Create feq.
for p=0:N
    for q=0:N
        feq(p+1,q+1)=((abs(p-q)*Z/N)<0.5)/3;
    end
end
%Predefine integrand as feq.
integrand=feq.*ones(N+1,N+1,Tlength);
for s=1:N
    for sp=1:N
        integrand(s+1,sp+1,:)=(f{1,1}(s,sp,:)+2*f{2,2}(s,sp,:));
    end
end
lambda=permute(sum(integrand,[1,2])...
    -0.5*sum(integrand([1,N+1],:,:),[1,2])...
    -0.5*sum(integrand(:,[1,N+1],:),[1,2])...
    +0.25*sum(integrand([1,N+1],[1,N+1],:),[1,2]),[1,3,2])*((Z/N).^2);
end

function xx=Sigmaxx(f)
%Sigma_xx.

%Calculate N.
N=size(f{1},1);
%Predefine integrand as feq.
integrand=ones(N+1,size(f{1},3))/3;
for s=1:N
    integrand(s+1,:)=f{1,1}(s,s,:);
end
Ge=1;
xx=((12*Ge)/5)*(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function yy=Sigmayy(f)
%Sigma_yy.

%Calculate N.
N=size(f{1},1);
%Predefine integrand as feq.
integrand=ones(N+1,size(f{1},3))/3;
for s=1:N
    integrand(s+1,:)=f{2,2}(s,s,:);
end
Ge=1;
yy=((12*Ge)/5)*(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function xy=Sigmaxy(f)
%Sigma_xy.

%Calculate N.
N=size(f{1},1);
%Predefine integrand as feq.
integrand=zeros(N+1,size(f{1},3));
for s=1:N
    integrand(s+1,:)=f{1,2}(s,s,:);
end
Ge=1;
xy=((12*Ge)/5)*(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function CNSD=CounterNSD(f)
%Counter Diagonal NSD.

%Calculate N.
N=size(f{1},1);
%Predefine integrand as zeros (feq-feq).
integrand=zeros(N+1,size(f{1},3));
for s=1:(N-1)
    integrand(s+1,:)=f{1,1}(s,N-s,:)-f{2,2}(s,N-s,:);
end
Ge=1;
CNSD=((12*Ge)/5)*(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function CNSD=CounterpNSD(f)
%Counter Diagonal NSD (with phantom).

warning('This functional is not currently is use, did you mean NSDCs?')
%Calculate N.
N=size(f{1},1);
%Predefine integrand as zeros (feq-feq).
integrand=zeros(N+1,size(f{1},3));
for s=1:(N-1)
    integrand(s+1,:)=f{1,1}(s,N-s,:)-f{2,2}(s,N-s,:);
end
%Calculate phantom point value.
N2L=floor(N/2);
N2U=ceil(N/2);
Ppoint=permute((f{1,1}(N2U,N2U,:)-f{2,2}(N2U,N2U,:)+...
    f{1,1}(N2L,N2L,:)-f{2,2}(N2L,N2L,:))/2,[1,3,2]);
Pintdiff=(0.5*Ppoint)-(0.25*(integrand(N2L+1,:)+integrand(N2U+1,:)));
%+0.5Ppoint-0.25middlepoints.
Ge=1;
CNSD=(12*Ge/5)*(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:))+...
    Pintdiff)/N;
end

function uu=UU(f)
%Normalised end to end functional.

%Set Dimensions.
Dim=2;
%Calculate N.
N=size(f{1},1);
%Calculate Tlen.
Tlen=size(f{1},3);
%Preallocate. (Form is Integrand(a,b,t,s) for s=s'.)
integrand=zeros(Dim,Dim,Tlen,N+1);
Trintf=zeros(1,1,Tlen);
%Predefine integrand as feq (for edges only).
for i=1:Dim
    integrand(i,i,:,:)=1/3;
end
%Create integrands
for i=1:Dim
    %Extract from f.
    for k=i:Dim
        for s=1:N
            integrand(i,k,:,s+1)=f{i,k}(s,s,:);
        end
    end
    %Symmetry speedup.
    for j=i:Dim
        integrand(j,i,:,:)=integrand(i,j,:,:);
    end
end
%Calculate uu with trapezium rule.
intf=(sum(integrand,4)-0.5*(integrand(:,:,:,1)+integrand(:,:,:,N+1)))/N;
if Dim==2
    for t=1:Tlen
        %Calculate trace. (ASSUMES yy=zz!!)
        Trintf(1,1,t)=trace(intf(:,:,t))+intf(2,2,t);
    end
elseif Dim==3
    for t=1:Tlen
        %Calculate trace.
        Trintf(1,1,t)=trace(intf(:,:,t));
    end
end
uu=intf./Trintf;
end

%%
function pqmatrix=Counterdiag(N)
%Logical matrix for q=Z-p. (note the exclusion of (0,N)=0 and (N,0)=0)
pqmatrix=zeros(N);
pqmatrix(1:N-1,1:N-1)=fliplr(eye(N-1));
end

function pqmatrix=z2(N)
%Logical matrix for p=q=z/2.
pqmatrix=zeros(N);
pqmatrix(round(N/2),round(N/2))=1;
end

%%
% function Var=PrinOrienVar(f)
% %Old principle orientation variance functional.
% 
% %Calculate N.
% N=size(f{1},1);
% %Predefine vector.
% Rad=zeros(N+1,size(f{1},3));
% %Calculate angle in radians from 11 direction.
% Rad(1,:)=acos(f{1,1}(N,N,:)./sqrt((f{1,1}(N,N,:).^2)+(f{2,2}(N,N,:).^2)));
% for s=1:N
%     Rad(s+1,:)=acos(f{1,1}(s,s,:)./sqrt((f{1,1}(s,s,:).^2)+(f{2,2}(s,s,:).^2)));
% end
% %Manually calculate var. (As we have data for |X| not X so E[X] is wrong.)
% %Note this is not sample var, as we already have N+1 points.
% Var=sum(Rad.^2,1)/(N+1);
% end