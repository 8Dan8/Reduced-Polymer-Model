function [fPred,G1hPred,G2hPred,G3hPred,SPred,RedTime] = Principle3Runge(Traindata,N_Tilde,Z_Tilde,Testdata,N,Z,Functionals,FlowType,StepET_R,cnu,Indexstep,FixFunc,CoFunc,InterSc,PlotOpt,Rolie,SmoothN,Metric)
%PRINCIPLE3RUNGE Uses three chosen functionals alongside the principle
%direction to predict the evolution of the GLaMM model.
%
%Inputs:
%Traindata - Cell array containing multiple input files. Accepts ET_R values, 
%            file name strings, output matrices (M) or {Time,f} cell pairs.
%N_Tilde - N value for the training data, vector or scalar input.
%Z_Tilde - Z value for the training data, vector or scalar input.
%Testdata - Cell array containing multiple input files. Accepts ET_R values, 
%           file name strings, output matrices (M) or {Time,f} cell pairs.
%N - N value for the test data, vector or scalar input.
%Z - Z value for the test data, vector or scalar input.
%Functionals - Cell containing names of functionals G_1, G_2 and G_3.
%FlowType - String detailinging type of polymer flow, 'xxUE','yyUE' or
%           'xyShear' are supported.
%StepET_R - Cell array {number of steps(int), ET_R(vec), Start times(vec)}.
%           Alternately, a string or integer containing the number of steps
%           maybe be used. This defaults step time to every 100 units and 
%           step size to ET_R/number of steps.
%cnu - Value of cnu for the polymer flow.
%Optional inputs:
%Indexstep - 2 by k matrix where each column comprises of a dt value and 
%            start time for the RK4 method. These columns are ran in 
%            sequence. Note that dt value is measured with respect to GLaMM
%            data, so the default of [1;0] will start the RK4 method at
%            time 0 and use a dt value equal to the ouput from the GLaMM
%            model.
%FixFunc - Vector of 3 logical values to fix gamma_1, gamma_2, gamma_3 or S
%          to their true values respectively. (Defaults to [0,0,0,0])
%CoFunc - String containing covariance function to use for interpolation.
%         (Defaults to 'Oct')
%InterSc - String for the scaling method to use during interpolation 
%          distance, 'Relative' or 'Absolute' are the currently available
%          options. (Defaults to 'Abs')
%PlotOpt - Vector of 5 logical values to optionally display S, dGamma/dt,
%          Interpolation, interpolation error or f (without rotation) on 
%          plots respectively. (Defaults to [0,0,0,0,0])
%Rolie - Vector of values for [Tau_e, beta]. If non-empty, plots a
%        Rolie-poly estimation for the current data.
%Metric - N by 1 cell vector containing metrics to be evaluated. (Empty by
%         default)

%Warning for yy.
if regexpi(FlowType,'yy\s*u\w*\s*e\w*')
    error('S,dS and df formulae assume yy=zz.')
end
%Check input format for Functionals.
if size(Functionals,1)~=1 && size(Functionals,1)~=3
    error(['Incorrect format for Functionals input, please use either a'...
        ,' cell row vector containing functional names or a cell matrix'...
        ,' containing function handles as in the Functional set output '...
        ,'from GLaMMFunctionalLinC.'])
end
%Setup defaults for multi step flows.
if iscell(StepET_R)==0
    %If non-cell, assume no steps.
    StepET_R={1,StepET_R,0};
end
%Sanity check.
if StepET_R{3}(1)~=0
    error('First step must start at t=0.')
end
%Setup other defualts.
if exist('Indexstep','var')==0
   Indexstep=[1;0];
elseif size(Indexstep,1)==1
    %Correct size of indexstep matrix.
    Indexstep=[Indexstep;zeros(1,numel(Indexstep))];
elseif size(Indexstep,1)~=2
    error('Indexstep input is of the wrong size.')
end
if exist('FixFunc','var')==0
    FixFunc=[];
elseif isempty(FixFunc)==0
    if size(FixFunc,1)==1
        FixFunc=[FixFunc;zeros(size(FixFunc))];
    end
    if size(FixFunc,2)==3
        FixFunc=[FixFunc,zeros(2,1)];
    end
end
if exist('CoFunc','var')==0 || isempty(CoFunc)
    CoFunc='Oct';
end
if exist('InterSc','var')==0 || isempty(InterSc)
    InterSc='Abs';
end
if exist('PlotOpt','var')==0 || isempty(PlotOpt)
    PlotS=0;
    PlotdG=0;
    PlotI=0;
    PlotIerr=0;
    Plotf=0;
else
    PlotS=PlotOpt(1);
    PlotdG=PlotOpt(2);
    PlotI=PlotOpt(3);
    PlotIerr=PlotOpt(4);
    Plotf=PlotOpt(5);
end
if exist('SmoothN','var')==0
    SmoothN=[];
elseif numel(SmoothN)==1
    SmoothN=SmoothN*ones(1,3);
end
%Create feq.
hat=zeros(N);
for p=1:N
    for q=1:N
        hat(p,q)=((abs(p-q)*Z/N)<0.5)/3;
    end
end
feq={hat,zeros(N);zeros(N),hat};
%Extract gamma_1, gamma_2, gamma_3 and f for both datasets.
[G1_Tilde,G2_Tilde,G3_Tilde,f_Tilde]=Functional3(Traindata,N_Tilde,...
    Z_Tilde,Functionals(1,:),0);
[G1,G2,G3,f,Time]=Functional3(Testdata,N,Z,Functionals(1,:),0);
%Define S function to compute orentation.
SFunc=GLaMMFunctional('S');
%Check for errors then extract 1st cell element.
if size(Time)>1
    error('Multiple test files not supported.')
end
Time=Time{1};
Tlen=numel(Time);
%Calculate dt.
dt=Time-[0,Time(1:(end-1))];
%Restrcture f and G_i for interpolation. (significant speedup)
%Restructure Gi_Tilde.
%Create vectors for G1_Tilde and G2_Tilde.
G1TVec=cell2mat(G1_Tilde);
G2TVec=cell2mat(G2_Tilde);
G3TVec=cell2mat(G3_Tilde);
%Detemine SampNo.
SampNo=numel(G1_Tilde);
%Preallocate for fTVec.
fTVec=cell(2,2);
%Restructure.
Temp=cellfun(@(x) x{1,1},f_Tilde, 'UniformOutput',false);
fTVec{1,1}=cat(3,Temp{:});
Temp=cellfun(@(x) x{2,2},f_Tilde, 'UniformOutput',false);
fTVec{2,2}=cat(3,Temp{:});
Temp=cellfun(@(x) x{1,2},f_Tilde, 'UniformOutput',false);
fTVec{1,2}=cat(3,Temp{:});
clear('Temp')
clear('f_Tilde')
fTVec{2,1}=fTVec{1,2};
%End Restructure section.

%Setup functional derivatives.
%Find number of functionals.
FuncNo=size(Functionals,2);
%Preallocate. (Need at least indentity for dS calc.)
dGkFunc=cell(1,FuncNo);
pqpart=cell(1,FuncNo);
pqmatrix=@(N) eye(N);
for i=1:FuncNo
    %Convert functional strings to function handles if nessecary.
    if isa(Functionals{1,i},'function_handle')==0
        %Functional name input.
        [dGkFunc{i},pqpart{i}]=GLaMMFunctionalddt(Functionals{1,i},Z,N);
    else
        error('Differental of Functional undefined.')
    end
    %Sum to find which pq elements need to be evaluated.
    pqmatrix=@(N) pqmatrix(N)+pqpart{i}(N);
end
%Define functional derivative formulas.
dG1Func=dGkFunc{1};
dG2Func=dGkFunc{2};
dG3Func=dGkFunc{3};
dSFunc=GLaMMFunctionalddt('S',Z,N);
%%
%True Gamma plot.
%Calculate fhat.
[fhat,S]=fPrinciple(f{1});
%Calculate functionals in principle direction.
[~,~,GkhatTemp]=FunctionalEval({Time,fhat},N,Z,Functionals(1,:));
G1hat=GkhatTemp{1};
G2hat=GkhatTemp{2};
G3hat=GkhatTemp{3};
%Setup colours.
Green=[0,1,0]; %GLaMM in lab frame.
Blue=[0,0.447,0.741]; %GLaMM in principle frame.
Red=[0.85,0.325,0.098]; %New model in lab frame.
Orange=[0.929,0.694,0.125]; %New model in principle frame.
Purple=[0.63,0.07,0.48];
%Setup and plot true Gamma values.
clf
subplot(1+PlotdG,3,1)
hold on
if Plotf
    plot(Time,G1{1},'color',Green,'displayname',...
        '$\Gamma_i$ from GLaMM predictions (dt=1/90)')
end
plot(Time,G1hat,'color',Blue,'displayname',...
    '$\dot{\Gamma}_i$ from GLaMM predictions (dt=1/90)')
hold off
xlabel('Time')
ylabel('\Gamma_1')
legend('interpreter','latex','location','southeast')
subplot(1+PlotdG,3,2)
hold on
if Plotf
    plot(Time,G2{1},'color',Green)
end
plot(Time,G2hat,'color',Blue)
hold off
xlabel('Time')
ylabel('\Gamma_2')
subplot(1+PlotdG,3,3)
hold on
if Plotf
    plot(Time,G3{1},'color',Green)
end
plot(Time,G3hat,'color',Blue)
hold off
xlabel('Time')
ylabel('\Gamma_3')
%Plot derivatives if requested.
if PlotdG
    if isempty(FixFunc)==0 && sum(FixFunc(2,:))
        %2) Find eigenmatrix V.
        V=Seigen(S);
        %4) Find the derivative of f using the GLaMM model on f.
        [df,~]=Multifderivs(f{1},Time,N,Z,FlowType,StepET_R,cnu,pqmatrix(N));
        %5) Find the derivative of S using its formula on f and df.
        dS=dSFunc(f{1},df);
        %6) Find the derivative of V using its formula on S and dS.
        dV=MultiVDerivs(S,dS);
        %7) Find the derivative of fhat using its formula on f,df,S and dS.
        dfhat=MultifhatDerivs(f{1},df,V,dV);
        %8) Find the derivatives of the functionals using fhat and dfhat.
        dG1hat=dG1Func(fhat,dfhat);
        dG2hat=dG2Func(fhat,dfhat);
        dG3hat=dG3Func(fhat,dfhat);
    else
        %Note this dG is good enough for a plot, but is NOT the GLaMM dG/dt.
        dG1hat=[NaN,(G1hat(2:end)-G1hat(1:(end-1)))./dt(2:end)];
        dG2hat=[NaN,(G2hat(2:end)-G2hat(1:(end-1)))./dt(2:end)];
        dG3hat=[NaN,(G3hat(2:end)-G3hat(1:(end-1)))./dt(2:end)];
    end
    %Plot true Gamma derivatives.
    subplot(2,3,4)
    plot(Time(2:end),dG1hat(2:end),'color',Blue)
    xlabel('Time')
    ylabel('d\Gamma_1/dt')
    subplot(2,3,5)
    plot(Time(2:end),dG2hat(2:end),'color',Blue)
    xlabel('Time')
    ylabel('d\Gamma_2/dt')
    subplot(2,3,6)
    plot(Time(2:end),dG3hat(2:end),'color',Blue)
    xlabel('Time')
    ylabel('d\Gamma_3/dt')
end
%%
%Setup for Illustrative plot.
if PlotI
    figure
    plot3(G1hat,G2hat,G3hat)
    hold on
    axis manual
    for Samp=1:SampNo
        plot3(G1_Tilde{Samp},G2_Tilde{Samp},G3_Tilde{Samp},...
            'color',0.6*ones(1,3),'displayname','');
    end
    grid on
    xlabel('\Gamma_1')
    ylabel('\Gamma_2')
    zlabel('\Gamma_3')
%     %Optional Zoom
%     axis([0,10,0,0.5,0,0.0005])
    pause
end
%%
%Pred f plot. (For each timestep)
%Preallocate fPred, G1Pred, G2Pred, G3Pred, dG1Pred, dG2Pred and dG3Pred 
%based on number of timesteps to evaluate.
fPred=cell(1,size(Indexstep,2));
G1hPred=fPred;
G2hPred=fPred;
G3hPred=fPred;
SPred=fPred;
dG1hPred=fPred;
dG2hPred=fPred;
dG3hPred=fPred;
dSPred=fPred;
NoPts=fPred;
MatchIndexVec=fPred;
if Plotf
    %Preallocate
    G1Pred=fPred;
    G2Pred=fPred;
    G3Pred=fPred;
end
if PlotIerr
    %Preallocate.
    G1hTest=fPred;
    G2hTest=fPred;
    G3hTest=fPred;
end
if Plotf || PlotIerr
    %Setup Functionals.
    G1Func=GLaMMFunctional(Functionals{1,1},Z,N);
    G2Func=GLaMMFunctional(Functionals{1,2},Z,N);
    G3Func=GLaMMFunctional(Functionals{1,3},Z,N);
end
%Loop for different timesteps.
for k=1:size(Indexstep,2)
    StartDesc=[];
    
    %Reduce number of entries based on Indexstep. (No change if 1)
    %Find GLaMM timestep.
    GTStep=round(Time(1),2,'significant');
    %Indices to evaluate with GLaMM model.
    GIndexVec=0:1:floor(Indexstep(2,k)/GTStep);
    %Indices to evaluate with RK4 model.
    EIndexVec=(floor(Indexstep(2,k)/GTStep)+Indexstep(1,k)):Indexstep(1,k):Tlen;
    %Create vector of indices.
    IndexVec=[GIndexVec,EIndexVec];
    
    %Ensure that RK4 time vector contains switch times.
    for i=2:size(StepET_R{3})
        %Find index of actual switch time.
        SwIndTrue=find(Time==StepET_R{3}(i),1);
        %Find index of first time after the switch. (Note index of index)
        SwIndInd=find(IndexVec>=SwIndTrue,1);
        %Check if timestep to be corrected occurs during RK4 method.
        if Time(IndexVec(SwIndInd))>floor(Indexstep(2,k)/GTStep)
            %Shift all indicies after this first occurance.
            IndexVec(SwIndInd:end)=IndexVec(SwIndInd:end)+(SwIndTrue-IndexVec(SwIndInd));
        end
    end
    
    %Create vector of times.
    if mod(Indexstep(1,k),1)==0
        %Select indices for integer steps.
        RedTime=[0,Time(IndexVec(2:end))];
    else
        %Calculate indices for fractional indices.
        RedTime=zeros(1,length(IndexVec));
        for j=2:length(IndexVec)
            if IndexVec(j)<=1
                RedTime(j)=IndexVec(j)*Time(1);
            else
                RedTime(j)=(mod(IndexVec(j),1)*Time(ceil(IndexVec(j))))+...
                    ((1-mod(IndexVec(j),1))*Time(floor(IndexVec(j))));
            end
        end
        warning('Extractions from GLaMM will fail with fractional indexstep.')
    end
    RedTlen=numel(RedTime);
    Reddt=RedTime(2:end)-RedTime(1:(end-1));
    
    %Preallocate fPred, each GhPred, and each dGhPred for this 
    %timestep. (Note the index 1 is t=0)
    fPred{k}={zeros(N,N,RedTlen),zeros(N,N,RedTlen);...
        zeros(N,N,RedTlen),zeros(N,N,RedTlen)};
    G1hPred{k}=zeros(1,RedTlen);
    G2hPred{k}=zeros(1,RedTlen);
    G3hPred{k}=zeros(1,RedTlen);
    SPred{k}=zeros(2,2,RedTlen);
    dG1hPred{k}=NaN*ones(1,RedTlen-1);
    dG2hPred{k}=NaN*ones(1,RedTlen-1);
    dG3hPred{k}=NaN*ones(1,RedTlen-1);
    dSPred{k}=NaN*ones(2,2,RedTlen-1);
    NoPts{k}=zeros(1,RedTlen);
    if PlotIerr
        G1hTest{k}=NaN(1,RedTlen-1);
        G2hTest{k}=NaN(1,RedTlen-1);
        G3hTest{k}=NaN(1,RedTlen-1);
    end
    %Create predictions at time 0. (Functional applied to feq)
    [~,fTemp,GkZero]=FunctionalEval({0,feq},N,Z,Functionals(1,:));
    G1hPred{k}(1)=GkZero{1};
    G2hPred{k}(1)=GkZero{2};
    G3hPred{k}(1)=GkZero{3};
    SPred{k}(:,:,1)=SFunc(fTemp);
    fPred{k}{1,1}(:,:,1)=fTemp{1,1};
    fPred{k}{1,2}(:,:,1)=fTemp{1,2};
    fPred{k}{2,1}(:,:,1)=fTemp{2,1};
    fPred{k}{2,2}(:,:,1)=fTemp{2,2};
    %Use true values upto start time.
    %-Indexstep(2,k) is start time.
    %-ceil(Indexstep(2,k)/10) is the final index for the GLaMM model.
    %-ceil(Indexstep(2,k)/(10*Indexstep(1,k)))+1 is the final GLaMM index for this code.
    GFinIndex=length(GIndexVec);
    if Indexstep(2,k)~=0
        %Extract from GLaMM model up to start time.
        for j=2:GFinIndex %j is in terms of GLaMM index.
            G1hPred{k}(j)=G1hat(IndexVec(j));
            G2hPred{k}(j)=G2hat(IndexVec(j));
            G3hPred{k}(j)=G3hat(IndexVec(j));
            SPred{k}(:,:,j)=S(:,:,IndexVec(j));
            fPred{k}{1,1}(:,:,j)=f{1}{1,1}(:,:,IndexVec(j));
            fPred{k}{1,2}(:,:,j)=f{1}{1,2}(:,:,IndexVec(j));
            fPred{k}{2,1}(:,:,j)=f{1}{2,1}(:,:,IndexVec(j));
            fPred{k}{2,2}(:,:,j)=f{1}{2,2}(:,:,IndexVec(j));
        end
        %Record new fTemp if non-zero start time.
        fTemp{1,1}=fPred{k}{1,1}(:,:,GFinIndex);
        fTemp{1,2}=fPred{k}{1,2}(:,:,GFinIndex);
        fTemp{2,1}=fPred{k}{2,1}(:,:,GFinIndex);
        fTemp{2,2}=fPred{k}{2,2}(:,:,GFinIndex);
        %Amend labels.
        StartDesc=['(dt=',num2str(round(mean(Reddt(GFinIndex:end)),3,'sig')),...
                ', Start Time=',num2str(Indexstep(2,k))];
    end
    %Record Raw data and preallocate for linear regression smoothing.
    if isempty(SmoothN)==0
        G1hRaw=G1hPred{k};
        G2hRaw=G2hPred{k};
        G3hRaw=G3hPred{k};
        LinPred=cell(1,3);
    end
    
%     SmoothG3{k}=G3hPred{k}
    %Iteratively caluclate G1, G2 and G3 predictions using RK4 method 
    %with orientation.
    if GFinIndex<RedTlen
        %----Inital step----(f,Gamma_i and S known)
        %Inital index is GFinIndex.
        i=GFinIndex;
        %2) Find eigenmatrix V.
        VTemp=Seigen(SPred{k}(:,:,i));
        %Inv3) Remove orentation from f to find fhat.
        fhatTemp=frotate(fTemp,permute(VTemp,[2,1,3]));
        %4) Find the derivative of f using the GLaMM model on f.
        dfTemp=Multifderivs(fTemp,RedTime(i),N,Z,FlowType,StepET_R,cnu,pqmatrix(N));
        %5) Find the derivative of S using its formula on f and df.
        k1S=dSFunc(fTemp,dfTemp);
        %6) Find the derivative of V using its formula on S and dS.
        dVTemp=MultiVDerivs(SPred{k}(:,:,i),k1S);
        %7) Find the derivative of fhat using its formula on f,df,S and dS.
        dfhatTemp=MultifhatDerivs(fTemp,dfTemp,VTemp,dVTemp);
        %8) Find the derivatives of the functionals using fhat and dfhat.
        k1G1=dG1Func(fhatTemp,dfhatTemp);
        k1G2=dG2Func(fhatTemp,dfhatTemp);
        k1G3=dG3Func(fhatTemp,dfhatTemp);
        %Find the values of k2, k3 and k4 for Runge Kutta.
        [k2G1,k2G2,k2G3,k2S] = Prin3dG(RedTime(i)+(Reddt(i))/2,...
            fTVec,G1TVec,G2TVec,G3TVec,...
            G1hPred{k}(i)+(Reddt(i).*k1G1/2),...
            G2hPred{k}(i)+(Reddt(i).*k1G2/2),...
            G3hPred{k}(i)+(Reddt(i).*k1G3/2),...
            SPred{k}(:,:,i)+(Reddt(i).*k1S/2),...
            N,Z,dG1Func,dG2Func,dG3Func,dSFunc,FlowType,StepET_R,cnu,...
            0,CoFunc,InterSc,pqmatrix);
        [k3G1,k3G2,k3G3,k3S] = Prin3dG(RedTime(i)+(Reddt(i))/2,...
            fTVec,G1TVec,G2TVec,G3TVec,...
            G1hPred{k}(i)+(Reddt(i).*k2G1/2),...
            G2hPred{k}(i)+(Reddt(i).*k2G2/2),...
            G3hPred{k}(i)+(Reddt(i).*k2G3/2),...
            SPred{k}(:,:,i)+(Reddt(i).*k2S/2),...
            N,Z,dG1Func,dG2Func,dG3Func,dSFunc,FlowType,StepET_R,cnu,...
            0,CoFunc,InterSc,pqmatrix);
        [k4G1,k4G2,k4G3,k4S] = Prin3dG(RedTime(i)+Reddt(i),...
            fTVec,G1TVec,G2TVec,G3TVec,...
            G1hPred{k}(i)+(Reddt(i).*k3G1),...
            G2hPred{k}(i)+(Reddt(i).*k3G2),...
            G3hPred{k}(i)+(Reddt(i).*k3G3),...
            SPred{k}(:,:,i)+(Reddt(i).*k3S),...
            N,Z,dG1Func,dG2Func,dG3Func,dSFunc,FlowType,StepET_R,cnu,...
            0,CoFunc,InterSc,pqmatrix);
        %Calculate Runge-Kutta average dGihPred and dSPred.
        dG1hPred{k}(i)=(k1G1+(2*k2G1)+(2*k3G1)+k4G1)/6;
        dG2hPred{k}(i)=(k1G2+(2*k2G2)+(2*k3G2)+k4G2)/6;
        dG3hPred{k}(i)=(k1G3+(2*k2G3)+(2*k3G3)+k4G3)/6;
        dSPred{k}(:,:,i)=(k1S+(2*k2S)+(2*k3S)+k4S)/6;
        %(Optional) Overwrite chosen functionals with true values.
        %(Use IndexVec to find corresponding GLaMM)
        if isempty(FixFunc)==0 && sum(FixFunc(2,:))
            dG1hPred{k}(i)=(1-FixFunc(2,1))*dG1hPred{k}(i)+...
                FixFunc(2,1)*dG1hat(IndexVec(i));
            dG2hPred{k}(i)=(1-FixFunc(2,2))*dG2hPred{k}(i)+...
                FixFunc(2,2)*dG2hat(IndexVec(i));
            dG3hPred{k}(i)=(1-FixFunc(2,3))*dG3hPred{k}(i)+...
                FixFunc(2,3)*dG3hat(IndexVec(i));
            dSPred{k}(:,:,i)=(1-FixFunc(2,4))*dSPred{k}(:,:,i)+...
                FixFunc(2,4)*dS(:,:,IndexVec(i));
        end
        %9) Perform RK4 method for all functionals. (At T=t -> T=t+dt)
        G1hPred{k}(i+1)=G1hPred{k}(i)+(Reddt(i)*dG1hPred{k}(i));
        G2hPred{k}(i+1)=G2hPred{k}(i)+(Reddt(i)*dG2hPred{k}(i));
        G3hPred{k}(i+1)=G3hPred{k}(i)+(Reddt(i)*dG3hPred{k}(i));
        SPred{k}(:,:,i+1)=SPred{k}(:,:,i)+(Reddt(i)*dSPred{k}(:,:,i));
        %(Optional) Smoothen functional predictions using N previous
        %values with linear regression.
        if isempty(SmoothN)==0
            %Add current Gi to Raw results.
            G1hRaw(i+1)=G1hPred{k}(i+1);
            G2hRaw(i+1)=G2hPred{k}(i+1);
            G3hRaw(i+1)=G3hPred{k}(i+1);
            %Use linear regression on previous N points.
            Y={G1hRaw(max(1,i+2-SmoothN(1)):i+1)',...
                G2hRaw(max(1,i+2-SmoothN(2)):i+1)',...
                G3hRaw(max(1,i+2-SmoothN(3)):i+1)'};
            for j=1:3
                if SmoothN(j)>1
                    X=RedTime(max(1,i+2-SmoothN(j)):i+1)';
                    %Quadratic
                    FX=[ones(numel(X),1),X,X.^2];
                    B=((FX'*FX)^-1)*FX'*Y{j};
                    LinPred{j}=B(1)+(B(2)*RedTime(i+1))+...
                        (B(3)*(RedTime(i+1)^2));
                    %line
%                     FX=[ones(numel(X),1),X];
%                     B=((FX'*FX)^-1)*FX'*Y{j};
%                     LinPred{j}=B(1)+(B(2)*RedTime(i+1));
                end
            end
            %Overwrite results.
            if SmoothN(1)>1, G1hPred{k}(i+1)=LinPred{1}; end
            if SmoothN(2)>1, G2hPred{k}(i+1)=LinPred{2}; end
            if SmoothN(3)>1, G3hPred{k}(i+1)=LinPred{3}; end
%             if SmoothN(3)>1, SmoothG3{k}(i+1)=LinPred{3}; end
        end
        %(Optional) Overwrite chosen functionals with true values.
        %(Use IndexVec to find corresponding GLaMM)
        if isempty(FixFunc)==0
            G1hPred{k}(i+1)=(1-FixFunc(1,1))*G1hPred{k}(i+1)+...
                FixFunc(1,1)*G1hat(IndexVec(i+1));
            G2hPred{k}(i+1)=(1-FixFunc(1,2))*G2hPred{k}(i+1)+...
                FixFunc(1,2)*G2hat(IndexVec(i+1));
            G3hPred{k}(i+1)=(1-FixFunc(1,3))*G3hPred{k}(i+1)+...
                FixFunc(1,3)*G3hat(IndexVec(i+1));
            SPred{k}(:,:,i+1)=(1-FixFunc(1,4))*SPred{k}(:,:,i+1)+...
                FixFunc(1,4)*S(:,:,IndexVec(i+1));
        end
        %----Iterative loop----(Gamma_i and S known)
        for i=GFinIndex+1:(RedTlen-1)
            %Find the values of k1-k4 for Runge Kutta for each functional.
            [k1G1,k1G2,k1G3,k1S,k1f] = Prin3dG(RedTime(i),...
                fTVec,G1TVec,G2TVec,G3TVec,...
                G1hPred{k}(i),...
                G2hPred{k}(i),...
                G3hPred{k}(i),...
                SPred{k}(:,:,i),...
                N,Z,dG1Func,dG2Func,dG3Func,dSFunc,FlowType,StepET_R,cnu,...
                (i<=PlotI),CoFunc,InterSc,pqmatrix);
            [k2G1,k2G2,k2G3,k2S,k2f] = Prin3dG(RedTime(i)+(Reddt(i))/2,...
                fTVec,G1TVec,G2TVec,G3TVec,...
                G1hPred{k}(i)+(Reddt(i).*k1G1/2),...
                G2hPred{k}(i)+(Reddt(i).*k1G2/2),...
                G3hPred{k}(i)+(Reddt(i).*k1G3/2),...
                SPred{k}(:,:,i)+(Reddt(i).*k1S/2),...
                N,Z,dG1Func,dG2Func,dG3Func,dSFunc,FlowType,StepET_R,cnu,...
                0,CoFunc,InterSc,pqmatrix);
            [k3G1,k3G2,k3G3,k3S,k3f] = Prin3dG(RedTime(i)+(Reddt(i))/2,...
                fTVec,G1TVec,G2TVec,G3TVec,...
                G1hPred{k}(i)+(Reddt(i).*k2G1/2),...
                G2hPred{k}(i)+(Reddt(i).*k2G2/2),...
                G3hPred{k}(i)+(Reddt(i).*k2G3/2),...
                SPred{k}(:,:,i)+(Reddt(i).*k2S/2),...
                N,Z,dG1Func,dG2Func,dG3Func,dSFunc,FlowType,StepET_R,cnu,...
                0,CoFunc,InterSc,pqmatrix);
            [k4G1,k4G2,k4G3,k4S,k4f] = Prin3dG(RedTime(i)+Reddt(i),...
                fTVec,G1TVec,G2TVec,G3TVec,...
                G1hPred{k}(i)+(Reddt(i).*k3G1),...
                G2hPred{k}(i)+(Reddt(i).*k3G2),...
                G3hPred{k}(i)+(Reddt(i).*k3G3),...
                SPred{k}(:,:,i)+(Reddt(i).*k3S),...
                N,Z,dG1Func,dG2Func,dG3Func,dSFunc,FlowType,StepET_R,cnu,...
                0,CoFunc,InterSc,pqmatrix);
            %Calculate Runge-Kutta average dGihPred and dSPred.
            dG1hPred{k}(i)=(k1G1+(2*k2G1)+(2*k3G1)+k4G1)/6;
            dG2hPred{k}(i)=(k1G2+(2*k2G2)+(2*k3G2)+k4G2)/6;
            dG3hPred{k}(i)=(k1G3+(2*k2G3)+(2*k3G3)+k4G3)/6;
            dSPred{k}(:,:,i)=(k1S+(2*k2S)+(2*k3S)+k4S)/6;
            %(Optional) Overwrite chosen functionals with true values.
            %(Use IndexVec to find corresponding GLaMM)
            if isempty(FixFunc)==0 && sum(FixFunc(2,:))
                dG1hPred{k}(i)=(1-FixFunc(2,1))*dG1hPred{k}(i)+...
                    FixFunc(2,1)*dG1hat(IndexVec(i));
                dG2hPred{k}(i)=(1-FixFunc(2,2))*dG2hPred{k}(i)+...
                    FixFunc(2,2)*dG2hat(IndexVec(i));
                dG3hPred{k}(i)=(1-FixFunc(2,3))*dG3hPred{k}(i)+...
                    FixFunc(2,3)*dG3hat(IndexVec(i));
                dSPred{k}(:,:,i)=(1-FixFunc(2,4))*dSPred{k}(:,:,i)+...
                    FixFunc(2,4)*dS(:,:,IndexVec(i));
            end
            %9) Perform RK4 method for all functionals. (At T=t -> T=t+dt)
            G1hPred{k}(i+1)=G1hPred{k}(i)+(Reddt(i)*dG1hPred{k}(i));
            G2hPred{k}(i+1)=G2hPred{k}(i)+(Reddt(i)*dG2hPred{k}(i));
            G3hPred{k}(i+1)=G3hPred{k}(i)+(Reddt(i)*dG3hPred{k}(i));
            SPred{k}(:,:,i+1)=SPred{k}(:,:,i)+(Reddt(i)*dSPred{k}(:,:,i));
            %(Optional) Reformat fTemp and record f predictions.
            fTemp=cellfun(@(k1,k2,k3,k4) (k1+(2*k2)+(2*k3)+k4)/6,...
                k1f,k2f,k3f,k4f,'uniformoutput',false);
            fPred{k}{1,1}(:,:,i)=fTemp{1,1};
            fPred{k}{1,2}(:,:,i)=fTemp{1,2};
            fPred{k}{2,1}(:,:,i)=fTemp{2,1};
            fPred{k}{2,2}(:,:,i)=fTemp{2,2};
            %(Optional) Smoothen functional predictions using N previous
            %values with linear regression.
            if isempty(SmoothN)==0
                %Add current Gi to Raw results.
                G1hRaw(i+1)=G1hPred{k}(i+1);
                G2hRaw(i+1)=G2hPred{k}(i+1);
                G3hRaw(i+1)=G3hPred{k}(i+1);
                %Use linear regression on previous N points.
                Y={G1hRaw(max(1,i+2-SmoothN(1)):i+1)',...
                    G2hRaw(max(1,i+2-SmoothN(2)):i+1)',...
                    G3hRaw(max(1,i+2-SmoothN(3)):i+1)'};
                for j=1:3
                    if SmoothN(j)>1
                        X=RedTime(max(1,i+2-SmoothN(j)):i+1)';
                        %Quadratic
                        FX=[ones(numel(X),1),X,X.^2];
                        B=((FX'*FX)^-1)*FX'*Y{j};
                        LinPred{j}=B(1)+(B(2)*RedTime(i+1))+...
                            (B(3)*(RedTime(i+1)^2));
                        %line
%                         FX=[ones(numel(X),1),X];
%                         B=((FX'*FX)^-1)*FX'*Y{j};
%                         LinPred{j}=B(1)+(B(2)*RedTime(i+1));
                    end
                end
                %Overwrite results.
                if SmoothN(1)>1, G1hPred{k}(i+1)=LinPred{1}; end
                if SmoothN(2)>1, G2hPred{k}(i+1)=LinPred{2}; end
                if SmoothN(3)>1, G3hPred{k}(i+1)=LinPred{3}; end
%                 if SmoothN(3)>1, SmoothG3{k}(i+1)=LinPred{3}; end
            end
            %(Optional) Overwrite chosen functionals with true values.
            if isempty(FixFunc)==0
                G1hPred{k}(i+1)=(1-FixFunc(1,1))*G1hPred{k}(i+1)+...
                    FixFunc(1,1)*G1hat(IndexVec(i+1));
                G2hPred{k}(i+1)=(1-FixFunc(1,2))*G2hPred{k}(i+1)+...
                    FixFunc(1,2)*G2hat(IndexVec(i+1));
                G3hPred{k}(i+1)=(1-FixFunc(1,3))*G3hPred{k}(i+1)+...
                    FixFunc(1,3)*G3hat(IndexVec(i+1));
                SPred{k}(:,:,i+1)=(1-FixFunc(1,4))*SPred{k}(:,:,i+1)+...
                    FixFunc(1,4)*S(:,:,IndexVec(i+1));
            end
            %(Optional) Code to check interpolation errors. (Gi->fhat->Gi)
            if PlotIerr
                G1hTest{k}(i)=G1Func(fhatTemp);
                G2hTest{k}(i)=G2Func(fhatTemp);
                G3hTest{k}(i)=G3Func(fhatTemp);
            end
%             %Code to extract fhat predictions----------------------------
%             fhatPred{k}{1,1}(:,:,i)=fhatTemp{1,1};
%             fhatPred{k}{1,2}(:,:,i)=fhatTemp{1,2};
%             fhatPred{k}{2,1}(:,:,i)=fhatTemp{2,1};
%             fhatPred{k}{2,2}(:,:,i)=fhatTemp{2,2};
        end
    end
    %(Optional) Reformat fTemp and record f predictions.
    fPred{k}{1,1}(:,:,RedTlen)=NaN;
    fPred{k}{1,2}(:,:,RedTlen)=NaN;
    fPred{k}{2,1}(:,:,RedTlen)=NaN;
    fPred{k}{2,2}(:,:,RedTlen)=NaN;
    
    %Setup vector containing indices for when our model conicides with GLAMM.
    MatchIndexVec{k}=[1:(GFinIndex-1),(GFinIndex-1+(1/Indexstep(1,k))):...
        (1/Indexstep(1,k)):(RedTlen-1)];
    
    %Cleanup for Illustrative plots.
    if PlotI
        hold off
        figure(1)
    end
    %Locate extrapolated data.
    if regexpi(CoFunc,'Oct')
        Polate=find(NoPts{k}==8,1,'last');
        if isempty(Polate)
            Polate=GFinIndex;
        end
    else
        Polate=RedTlen;
    end
    %Plot RK4 using predicted Gamma.
    subplot(1+PlotdG,3,1)
    hold on
    if Plotf
        G1Pred{k}=G1Func(fPred{k});
        plot(RedTime(1:Polate),G1Pred{k}(1:Polate),...
            'color',Red,...
            'displayname',['$\Gamma_i$ from model predictions ',...
            StartDesc,')'])
        if Polate~=RedTlen
            plot(RedTime(Polate:(end-1)),G1Pred{k}(Polate:(end-1)),...
                'linestyle','--','color',Red,...
                'displayname',['$\Gamma_i$ from model predictions ',...
                StartDesc,', extrapolated)'])
        end
    end
    plot(RedTime(1:Polate),G1hPred{k}(1:Polate),...
        'color',Orange,...
        'displayname',['$\dot{\Gamma}_i$ from model predictions ',...
        StartDesc,')'])
    if Polate~=RedTlen
        plot(RedTime(Polate:end),G1hPred{k}(Polate:end),...
            'linestyle','--','color',Orange,...
            'displayname',['$\dot{\Gamma}_i$ from model predictions ',...
            StartDesc,', extrapolated)'])
    end
    hold off
    subplot(1+PlotdG,3,2)
    hold on
    if Plotf
        G2Pred{k}=G2Func(fPred{k});
        plot(RedTime(1:Polate),G2Pred{k}(1:Polate),'color',Red)
        if Polate~=RedTlen
            plot(RedTime(Polate:(end-1)),G2Pred{k}(Polate:(end-1)),...
                'linestyle','--','color',Red)
        end
    end
    plot(RedTime(1:Polate),G2hPred{k}(1:Polate),'color',Orange)
    if Polate~=RedTlen
        plot(RedTime(Polate:end),G2hPred{k}(Polate:end),'linestyle','--',...
            'color',Orange)
    end
    hold off
    subplot(1+PlotdG,3,3)
    hold on
    if Plotf
        G3Pred{k}=G3Func(fPred{k});
        plot(RedTime(1:Polate),G3Pred{k}(1:Polate),'color',Red)
        if Polate~=RedTlen
            plot(RedTime(Polate:(end-1)),G3Pred{k}(Polate:(end-1)),...
                'linestyle','--','color',Red)
        end
    end
    plot(RedTime(1:Polate),G3hPred{k}(1:Polate),'color',Orange)
    if Polate~=RedTlen
        plot(RedTime(Polate:end),G3hPred{k}(Polate:end),'linestyle','--',...
            'color',Orange)
    end
    hold off
    if PlotdG
        subplot(2,3,4)
        hold on
        %This commented method may not be suitable.
%         if Plotf
%             plot(RedTime(1:(end-2)),(G1Pred{k}(2:(end-1))-G1Pred{k}(1:(end-2)))./Reddt(1:end-1),'color',Red)
%         end
        plot(RedTime(1:(end-1)),dG1hPred{k},'color',Orange)
        hold off
        subplot(2,3,5)
        hold on
%         if Plotf
%             plot(RedTime(1:(end-2)),(G2Pred{k}(2:(end-1))-G2Pred{k}(1:(end-2)))./Reddt(1:end-1),'color',Red)
%         end
        plot(RedTime(1:(end-1)),dG2hPred{k},'color',Orange)
        hold off
        subplot(2,3,6)
        hold on
%         if Plotf
%             plot(RedTime(1:(end-2)),(G3Pred{k}(2:(end-1))-G3Pred{k}(1:(end-2)))./Reddt(1:end-1),'color',Red)
%         end
        plot(RedTime(1:(end-1)),dG3hPred{k},'color',Orange)
        hold off
    end
end
%%
% Rolie-Poly plot.
% %Override for transformation of G3.
% fhG3trans=Functional3Interpolator(fTVec,G1TVec,G2TVec,G3TVec,G1hPred{k},G2hPred{k},[GkZero{3},G3hat(IndexVec(2:end))],N,0,CoFunc,InterSc);
% VG3trans=Seigen(SPred{k});
% fPred{k}=frotate(fhG3trans,VG3trans);
% %Override for Smoothen G3 in transformation only.
% fhG3trans=Functional3Interpolator(fTVec,G1TVec,G2TVec,G3TVec,G1hPred{k},G2hPred{k},SmoothG3{k},N,0,CoFunc,InterSc);
% VG3trans=Seigen(SPred{k});
% fPred{k}=frotate(fhG3trans,VG3trans);
if exist('Rolie','var') && isempty(Rolie)==0
    %Define parameters.
    Taue=Rolie(1);
    beta=Rolie(2);
    %Preallocate.
    RTime=cell(1,StepET_R{1});
    Kappa=cell(1,StepET_R{1});
    RolieNSD=cell(1,StepET_R{1});
    RolieSigxy=cell(1,StepET_R{1});
    %Split time if multiple sources.
    for i=1:StepET_R{1}-1
        RTime{i}=[StepET_R{3}(i),...
            Time((StepET_R{3}(i)+0.01)<Time & Time<StepET_R{3}(i+1))];
    end
    RTime{StepET_R{1}}=[StepET_R{3}(StepET_R{1}),...
        Time((StepET_R{3}(StepET_R{1})+0.01)<Time)];
    %Create Kappa.
    TauR=Taue*(Z^2);
    for i=1:StepET_R{1}
        if regexpi(FlowType,'xx\s*u\w*\s*e\w*')
            %xx Extension [e,0,0;0,-e/2,0;0,0,-e/2]
            Kappa{i}=[StepET_R{2}(i)/TauR,0;0,-StepET_R{2}(i)/(2*TauR)];
        elseif regexpi(FlowType,'yy\s*u\w*\s*e\w*')
            %yy Extension [-e/2,0,0;0,e,0;0,0,-e/2]
            error('Trace assumes yy=zz in Rolie.')
        elseif regexpi(FlowType,'xy\s*sh\w*')
            %xy Shear [0,g,0;0,0,0;0,0,0]
            Kappa{i}=[0,StepET_R{2}(i)/TauR;0,0];
        else
            error('Unrecognised FlowType.')
        end
    end
    %Calculate the Rolie-Poly model.
    RolieSig=Rolie2D(RTime,Kappa,Z,Taue,beta);
    for i=1:StepET_R{1}
        RolieNSD{i}=(RolieSig{i}(1,:)-RolieSig{i}(2,:));
        RolieSigxy{i}=RolieSig{i}(3,:);
    end
    %Define functional.
    NSDFunc=GLaMMFunctional('NSD',Z,N);
    SigxyFunc=GLaMMFunctional('Sigxy',Z,N);
    %NSD plots.
    figure
    subplot(1,2,1)
    hold on
    %Plot GLaMM.
    NSDTrue=NSDFunc(f{1});
    plot(Time,NSDTrue,'color',Green,'displayname',...
        'NSD from GLaMM predictions (dt=1/90)')
    %Plot model predcitions.
    %Note this only displays the final indexstep result, a k loop could be 
    %added to display all.
    NSDPred=NSDFunc(fPred{k});
    plot(RedTime(1:Polate),NSDPred(1:Polate),...
        'color',Red,...
        'displayname',['NSD from model predictions ',...
        StartDesc,')'])
    if Polate~=RedTlen
        plot(RedTime(Polate:(end-1)),NSDPred(Polate:(end-1)),...
            'linestyle','--','color',Red,...
            'displayname',['NSD from model predictions ',...
            StartDesc,', extrapolated)'])
    end
    %Plot Rolie-Poly model.
    plot(cell2mat(RTime),cell2mat(RolieNSD),'color',Purple,...
        'displayname','NSD from Rolie-Poly model')
    hold off
    xlabel('Time')
    ylabel('NSD')
    legend

    %Sxy plots.
    subplot(1,2,2)
    hold on
    %Plot GLaMM.
    SigxyTrue=SigxyFunc(f{1});
    plot(Time,SigxyTrue,'color',Green,'displayname',...
        '\sigma_{xy} from GLaMM predictions (dt=1/90)')
    %Plot model predcitions.
    SigxyPred=SigxyFunc(fPred{k});
    plot(RedTime(1:Polate),SigxyPred(1:Polate),...
        'color',Red,...
        'displayname',['\sigma_{xy} from model predictions ',...
        StartDesc,')'])
    if Polate~=RedTlen
        plot(RedTime(Polate:(end-1)),SigxyPred(Polate:(end-1)),...
            'linestyle','--','color',Red,...
            'displayname',['\sigma_{xy} from model predictions ',...
            StartDesc,', extrapolated)'])
    end
    %Plot Rolie-Poly model.
    plot(cell2mat(RTime),cell2mat(RolieSigxy),'color',Purple,...
        'displayname','\sigma_{xy} from Rolie-Poly model')
    hold off
    xlabel('Time')
    ylabel('\sigma_{xy}')
    legend
end
%%
% Sxx and Sxy plots.
if PlotS
    figure
    %Preallocate.
    SxxPred=cell(1,size(Indexstep,2));
    SxyPred=SxxPred;
    dSxxPred=SxxPred;
    dSxyPred=SxxPred;
    %Define variables.
    Sxx=permute(S(1,1,:),[3,1,2]);
    Sxy=permute(S(1,2,:),[3,1,2]);
    for k=1:size(Indexstep,2)
        SxxPred{k}=permute(SPred{k}(1,1,:),[3,1,2]);
        SxyPred{k}=permute(SPred{k}(1,2,:),[3,1,2]);
        dSxxPred{k}=permute(dSPred{k}(1,1,:),[3,1,2]);
        dSxyPred{k}=permute(dSPred{k}(1,2,:),[3,1,2]);
    end
    %Plot true S values.
    subplot(1+PlotdG,2,1)
    plot(Time,Sxx,'color',Green,...
        'displayname','$S$ from GLaMM predictions (dt=1/90)')
    if Plotf
        %Only for testing.
        hold on
        [~,Shat]=fPrinciple(fhat);
        plot(Time,permute(Shat(1,1,:),[3,1,2]),'color',Blue,...
            'displayname','$\dot{S}$ from GLaMM predictions (dt=1/90)')
        hold off
    end
    xlabel('Time')
    ylabel('S_{xx}')
    legend('interpreter','latex','location','southeast')
    subplot(1+PlotdG,2,2)
    plot(Time,Sxy,'color',Green)
    if Plotf
        hold on
        plot(Time,permute(Shat(1,2,:),[3,1,2]),'color',Blue)
        hold off
    end
    xlabel('Time')
    ylabel('S_{xy}')
    if PlotdG
        %Note this dS is good enough for a plot, but is NOT the GLaMM dS/dt.
        subplot(2,2,3)
        plot(Time(2:end),((Sxx(2:end)-Sxx(1:(end-1)))'./dt(2:end)),...
            'color',Green)
        xlabel('Time')
        ylabel('dS_{xx}/dt')
        subplot(2,2,4)
        plot(Time(2:end),((Sxy(2:end)-Sxy(1:(end-1)))'./dt(2:end)),...
            'color',Green)
        xlabel('Time')
        ylabel('dS_{xy}/dt')
    end
    %Plot predicted S values.
    for k=1:size(Indexstep,2)
        subplot(1+PlotdG,2,1)
        hold on
        plot(RedTime(1:Polate),SxxPred{k}(1:Polate),'color',Red,...
            'displayname',['$S$ from model predictions ',...
            StartDesc,')'])
        if Polate~=RedTlen
            plot(RedTime(Polate:end),SxxPred{k}(Polate:end),...
                'linestyle','--','color',Red,...
                'displayname',['$S$ from model predictions ',...
                StartDesc,', extrapolated)'])
        end
        hold off
        subplot(1+PlotdG,2,2)
        hold on
        plot(RedTime(1:Polate),SxyPred{k}(1:Polate),'color',Red)
        if Polate~=RedTlen
            plot(RedTime(Polate:end),SxyPred{k}(Polate:end),...
                'linestyle','--','color',Red)
        end
        hold off
        if PlotdG
            subplot(2,2,3)
            hold on
            plot(RedTime(1:(end-1)),dSxxPred{k},'color',Red)
            hold off
            subplot(2,2,4)
            hold on
            plot(RedTime(1:(end-1)),dSxyPred{k},'color',Red)
            hold off
        end
    end
end
%%
% Interpolation error plot. 
if PlotIerr
    figure
    title('Interpolation error in $\dot\Gamma_i$')
    subplot(1,3,1)
    hold on
    for k=1:size(Indexstep,2)
        plot(RedTime(1:end-1),G1hTest{k}-G1hPred{k}(1:end-1),'color',Orange)
    end
    hold off
    xlabel('Time')
    ylabel('Error in $\dot\Gamma_1$','interpreter','latex')
    subplot(1,3,2)
    hold on
    for k=1:size(Indexstep,2)
        plot(RedTime(1:end-1),G2hTest{k}-G2hPred{k}(1:end-1),'color',Orange)
    end
    hold off
    xlabel('Time')
    ylabel('Error in $\dot\Gamma_2$','interpreter','latex')
    subplot(1,3,3)
    hold on
    for k=1:size(Indexstep,2)
        plot(RedTime(1:end-1),G3hTest{k}-G3hPred{k}(1:end-1),'color',Orange)
    end
    hold off
    xlabel('Time')
    ylabel('Error in $\dot\Gamma_3$','interpreter','latex')
end
%%
% Evaluate metric.
if exist('Metric','var') && isempty(Metric)==0
    error('Using older Metric code without orentation or G3.')
    for k=1:size(Indexstep,2)
        %Error Checking
        if Indexstep(1,k)>1
            error('Metric evaluation not currently calculable for Indexstep>1');
        end
        figure
        hold on
        %Calculate number of metrics.
        NoMetric=size(Metric,2);
        %Preallocate TScore and Score.
        TScore=cell(1,NoMetric);
        Score=zeros(1,NoMetric);
        for m=1:NoMetric
            if regexpi(Metric{1,m},'dG.dt')
                %dG_i/dt metric input. {'dG.dt',FlowType,StepET_R}
                if str2double(Metric{1,m}(3))==1
                    dGPred=dG1hPred{k};
                elseif str2double(Metric{1,m}(3))==2
                    dGPred=dG2hPred{k};
                else
                    error('Unrecognised Metric format.');
                end
                %Extract dG predcitions at times when true data exists.
                dGPredRed=dGPred(MatchIndexVec{k});
                %Extract true dG from GLaMM data.
                dGZero=dGdt(Functionals(:,str2double(Metric{1,m}(3))),...
                    feq,0,N,Z,FlowType,StepET_R,cnu);
                dGTrue=dGdt(Functionals(:,str2double(Metric{1,m}(3))),...
                    f{1},Time,N,Z,FlowType,StepET_R,cnu);
                %Reformat true dG data and include t=0.
                dGTrueRed=[dGZero{1},dGTrue{1}(1:(end-1))];
                %Overwrite if delayed start time.
                dGPredRed(1:(GFinIndex-1))=dGTrueRed(1:(GFinIndex-1));
                %Define metric times.
                MetTime=[0,Time(1:end-1)];
                %Calculate error over time. (Inc. t=0, Exc. t=end)
                %Change to relative dGdt if requested.
                if regexpi(Metric{1,m},'Rel')
                    MetricName=['Relative error in d\Gamma_',...
                        Metric{1,m}(3),'/dt '];
                    TScore{m}=(dGPredRed-dGTrueRed)./dGTrueRed;
                else
                    MetricName=['Error in d\Gamma_',Metric{1,m}(3),'/dt '];
                    TScore{m}=(dGPredRed-dGTrueRed);
                end
            elseif regexpi(Metric{1,m},'SumdGdt')
                %Normalised dG/dt sum metric input.
                MetricName='Sum of |error(d\Gamma_i/dt)/\Gamma_i|';
                %Extract dG predcitions at times when true data exists.
                dG1PredRed=dG1hPred{k}(MatchIndexVec{k});
                dG2PredRed=dG2hPred{k}(MatchIndexVec{k});
                %Extract true dG from GLaMM data.
                dGkZero=dGdt(Functionals,feq,0,N,Z,FlowType,StepET_R,cnu);
                dGkTrue=dGdt(Functionals,f{1},Time,N,Z,FlowType,StepET_R,cnu);
                %Reformat true G data and include t=0.
                G1Red=[GkZero{1},G1{1}(1:(end-1))];
                G2Red=[GkZero{2},G2{1}(1:(end-1))];
                %Reformat true dG data and include t=0.
                dG1TrueRed=[dGkZero{1},dGkTrue{1}(1:(end-1))];
                dG2TrueRed=[dGkZero{2},dGkTrue{2}(1:(end-1))];
                %Overwrite if delayed start time.
                dG1PredRed(1:(GFinIndex-1))=dG1TrueRed(1:(GFinIndex-1));
                dG2PredRed(1:(GFinIndex-1))=dG2TrueRed(1:(GFinIndex-1));
                %Calculate error over time. (Inc. t=0, Exc. t=end)
                TScore{m}=abs((dG1PredRed-dG1TrueRed)./G1Red)+...
                    abs((dG2PredRed-dG2TrueRed)./G2Red);
                MetTime=[0,Time(1:end-1)];
            else
                if ischar(Metric{1,m})
                    %Calculate metric for string inputs.
                    MetricName=Metric{1,m};
                    MetricFunc=GLaMMMetric(Metric{1,m},Z_Tilde,Z);
                else
                    %Direct metric input.
                    MetricName=func2str(Metric{1,m});
                    MetricFunc=Metric{1,m};
                end
                %Evaluate metric at each time which a GLaMM point exists.
                fPredRed={fPred{k}{1,1}(:,:,MatchIndexVec{k}+1),...
                    fPred{k}{1,2}(:,:,MatchIndexVec{k}+1);...
                    fPred{k}{2,1}(:,:,MatchIndexVec{k}+1),...
                    fPred{k}{2,2}(:,:,MatchIndexVec{k}+1)};
                TScore{m}=MetricFunc(fPredRed,f{1});
                MetTime=Time;
            end
            %Calculate score by taking RMS with repect to time.
            Score(m)=sqrt(nanmean((TScore{m}(GFinIndex:end)).^2));
            %Plot the prediction data.
            plot(MetTime,TScore{m},'DisplayName',[MetricName,' (RMS = ',num2str(Score(m)),')'])
        end
        hold off
        title(['Prediction ',num2str(k)])
        xlabel('Time')
        ylabel('Value of Metric')
        legend('location','northwest')
    end
end
end