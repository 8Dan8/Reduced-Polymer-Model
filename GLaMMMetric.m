function [Metric,Plotfuncs,Plotlabels] = GLaMMMetric(Name,Z_Tilde,Z)
%GLAMMMETRIC attempts to find a metric called 'name'.
%
%List of defined Metrics:
%'StressRMSE' - Root mean square error for each element of stress.
%(Root of sqerr(sigxx)+sqerr(sigyy)+2*sqerr(sigxy)/4 over time.)
%'StressSE' - Sum of square error for each element of stress.
%(sqerr(sigxx)+sqerr(sigyy)+2*sqerr(sigxy)/4 over time.)
%'DiagStressRMSE' - Root mean square error for diagonal elements of stress.
%(Root of sqerr(sigxx)+sqerr(sigyy)/2 over time.)
%'RelativefVolumeRMSE' - Normalised root mean square error of the volume under f(s,s').
%(Root of sum((int(err(f))/int(f)).^2,[i,j])/4 over time.)
%'fRMSE' - Root mean square error for all elements of f.
%(Root of sum(sqerr(f),[i,j,s,s'])/(4*(N^2)) over time.)
%'NSD' - Error in Normal stress difference.

%Should we consider both symmetrical elements?

if regexpi(Name,'Stress\s*RMSE')
    %Root mean square error of stress.
    Metric=@(f_pred,f) StressRMSE(f_pred,f,Z_Tilde,Z);
    Plotfuncs={@(f_pred) Sigmaxx(f_pred,Z_Tilde),...
        @(f_pred) Sigmayy(f_pred,Z_Tilde),...
        @(f_pred) Sigmaxy(f_pred,Z_Tilde);...
        @(f) Sigmaxx(f,Z),@(f) Sigmayy(f,Z),@(f) Sigmaxy(f,Z)};
    Plotlabels={'$\sigma_{xx}$','$\sigma_{yy}$','$\sigma_{xy}$'};
elseif regexpi(Name,'Stress\s*SE')
    %Square error of stress.
    Metric=@(f_pred,f) StressSE(f_pred,f,Z_Tilde,Z);
    Plotfuncs={@(f_pred) Sigmaxx(f_pred,Z_Tilde),...
        @(f_pred) Sigmayy(f_pred,Z_Tilde),...
        @(f_pred) Sigmaxy(f_pred,Z_Tilde);...
        @(f) Sigmaxx(f,Z),@(f) Sigmayy(f,Z),@(f) Sigmaxy(f,Z)};
    Plotlabels={'$\sigma_{xx}$','$\sigma_{yy}$','$\sigma_{xy}$'};
elseif regexpi(Name,'Diag\w*\s*Stress\s*RMSE')
    %Root mean square error of diagonal elements of stress.
    Metric=@(f_pred,f) DiagStressRMSE(f_pred,f,Z_Tilde,Z);
    Plotfuncs={@(f_pred) Sigmaxx(f_pred,Z_Tilde),...
        @(f_pred) Sigmayy(f_pred,Z_Tilde);...
        @(f) Sigmaxx(f,Z),@(f) Sigmayy(f,Z)};
    Plotlabels={'$\sigma_{xx}$','$\sigma_{yy}$'};
elseif regexpi(Name,'Rel\w*\s*f\s*\w*\s*RMSE')
    %Normalised root mean square error each volume of f.
    Metric=@(f_pred,f) RelfVolRMSE(f_pred,f,Z_Tilde,Z);
    Plotfuncs={};
    Plotlabels={};
elseif regexpi(Name,'f\s*RMSE')
    %Root mean square error of all elements of f.
    Metric=@(f_pred,f) fRMSE(f_pred,f,Z_Tilde,Z);
    Plotfuncs={};
    Plotlabels={};
elseif regexpi(Name,'N\w*\s*S\w*\s*D')
    %1st normal stress difference error.
    Metric=@(f_pred,f) abs((Sigmaxx(f_pred,Z_Tilde)-Sigmayy(f_pred,Z_Tilde))...
        -(Sigmaxx(f,Z)-Sigmayy(f,Z)));
    Plotfuncs={@(f_pred) Sigmaxx(f_pred,Z_Tilde),...
        @(f_pred) Sigmayy(f_pred,Z_Tilde);...
        @(f) Sigmaxx(f,Z),@(f) Sigmayy(f,Z)};
else
    %Error catching.
    error('Metric not defined.')
end
end

%%
%Metrics
function RMSE=StressRMSE(f_pred,f,Z_Tilde,Z)
%Stress root mean square error.

%Calculate error for each sigma element.
Exx=abs(Sigmaxx(f_pred,Z_Tilde)-Sigmaxx(f,Z));
Exy=abs(Sigmaxy(f_pred,Z_Tilde)-Sigmaxy(f,Z));
Eyy=abs(Sigmayy(f_pred,Z_Tilde)-Sigmayy(f,Z));
%Calculate RMSE, ignoring NaN values.
RMSE=sqrt((Exx.^2+2*(Exy.^2)+Eyy.^2)./4);
end

function SE=StressSE(f_pred,f,Z_Tilde,Z)
%Stress square error.

%Calculate error for each sigma element.
Exx=abs(Sigmaxx(f_pred,Z_Tilde)-Sigmaxx(f,Z));
Exy=abs(Sigmaxy(f_pred,Z_Tilde)-Sigmaxy(f,Z));
Eyy=abs(Sigmayy(f_pred,Z_Tilde)-Sigmayy(f,Z));
%Calculate SE, ignoring NaN values.
SE=(Exx.^2+2*(Exy.^2)+Eyy.^2)./4;
end

function RMSE=DiagStressRMSE(f_pred,f,Z_Tilde,Z)
%Trace(Stress) root mean square error.

%Calculate error for each sigma element.
Exx=abs(Sigmaxx(f_pred,Z_Tilde)-Sigmaxx(f,Z));
Eyy=abs(Sigmayy(f_pred,Z_Tilde)-Sigmayy(f,Z));
%Calculate RMSE, ignoring NaN values.
RMSE=sqrt((Exx.^2+Eyy.^2)./2);
end

function RMSE=RelfVolRMSE(f_pred,f,Z_Tilde,Z)
%Normalised f(s,s') volume root mean square error.

%Error checking.
if Z_Tilde~=Z
    error('Z must be equal between the two data sets to use the RelfVolRMSE metric.')
end
%Calculate error for each sigma element.
Efxx=abs(f_pred{1,1}-f{1,1});
Efxy=abs(f_pred{1,2}-f{1,2});
Efyy=abs(f_pred{2,2}-f{2,2});
Ef={Efxx,Efxy;Efxy,Efyy};
%Calculate f volumes.
fVol=fVolIntegral(f,Z);
EfVol=fVolIntegral(Ef,Z);
%Calculate RMSE, normalising f volumes.
RMSE=sqrt(nansum((EfVol./fVol).^2,1)./4);
end

function RMSE=fRMSE(f_pred,f,Z_Tilde,Z)
%f(s,s') root mean square error.

%Error checking.
if Z_Tilde~=Z
    error('Z must be equal between the two data sets to use the fRMSE metric.')
end
%Calculate error for each sigma element.
Efxx=abs(f_pred{1,1}-f{1,1});
Efxy=abs(f_pred{1,2}-f{1,2});
Efyy=abs(f_pred{2,2}-f{2,2});
%Calculate RMSE, ignoring NaN values.
RMSE=permute(sqrt(sum(Efxx.^2+2*(Efxy.^2)+Efyy.^2,[1,2])./(4*(Z^2))),[1,3,2]);
end

%%
%Functionals
function xx=Sigmaxx(f,Z)
%Sigma_xx.

%Calculate N.
N=size(f{1},1);
%Predefine integrand as feq.
integrand=ones(N+1,size(f{1},3));
for s=1:N
    integrand(s+1,:)=f{1,1}(s,s,:);
end
Ge=1;
xx=((12*Ge)/(5*Z))*(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function yy=Sigmayy(f,Z)
%Sigma_yy.

%Calculate N.
N=size(f{1},1);
%Predefine integrand as feq.
integrand=ones(N+1,size(f{1},3));
for s=1:N
    integrand(s+1,:)=f{2,2}(s,s,:);
end
Ge=1;
yy=((12*Ge)/(5*Z))*(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function xy=Sigmaxy(f,Z)
%Sigma_xy.

%Calculate N.
N=size(f{1},1);
%Predefine integrand as feq.
integrand=ones(N+1,size(f{1},3));
for s=1:N
    integrand(s+1,:)=f{1,2}(s,s,:);
end
Ge=1;
xy=((12*Ge)/(5*Z))*(sum(integrand,1)-0.5*(integrand(1,:)+integrand(N+1,:)))/N;
end

function fvol=fVolIntegral(f,Z)
%Total trace integral functional.

%Calculate N and time length.
[N,~,Tlength]=size(f{1});
%Preallocation.
feq=zeros(N+1);
integrand=cell(2);
fvol=zeros(4,Tlength);
%Create feq.
for p=0:N
    for q=0:N
        feq(p+1,q+1)=((abs(p-q)*Z/N)<0.5)/3;
    end
end
%Predefine integrand as feq.
integrand{1}=feq.*ones(N+1,N+1,Tlength);
integrand{4}=feq.*ones(N+1,N+1,Tlength);
integrand{2}=zeros(N+1,N+1,Tlength);
%loop over i,j.
for i=[1,2,4]
    for s=1:N
        for sp=1:N
            integrand{i}(s+1,sp+1,:)=(f{i}(s,sp,:));
        end
    end
    fvol(i,:)=permute(sum(integrand{i},[1,2])...
        -0.5*sum(integrand{i}([1,N+1],:,:),[1,2])...
        -0.5*sum(integrand{i}(:,[1,N+1],:),[1,2])...
        +0.25*sum(integrand{i}([1,N+1],[1,N+1],:),[1,2]),[1,3,2])/(N.^2);
end
%Apply symmetry.
fvol(3,:)=fvol(2,:);
end