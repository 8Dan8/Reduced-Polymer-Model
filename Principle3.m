function [fPred,G1hPred,G2hPred,SPred,RedTime] = Principle3(Traindata,N_Tilde,Z_Tilde,Testdata,N,Z,Functionals,CoFunc,InterSc,PlotOpt,Rolie,FlowType,StepET_R,cnu,RescV,TransMethod)
%PRINCIPLE3 Uses three chosen functionals alongside the principle direction
%to predict the GLaMM model given the true values of the functionals.
%This is comparable to PRINCIPLEEULER with fixfunc [1,1,1] in 2D.
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
%Optional inputs:
%CoFunc - String containing covariance function to use for interpolation.
%         (Defaults to 'Oct')
%InterSc - String for the scaling method to use during interpolation 
%          distance, 'Relative' or 'Absolute' are the currently available
%          options. (Defaults to 'Abs')
%PlotOpt - Vector of 4 logical values to optionally display S, dGamma/dt,
%          Interpolation or f (without rotation) on plots respectively. 
%          (Defaults to [0,0,0,0])
%Rolie - Vector of values for [Tau_e, beta]. If non-empty, plots a
%        Rolie-poly estimation for the current data.
%FlowType - String detailing type of polymer flow, 'xxUE','yyUE' or
%           'xyShear' are supported. (Required for Rolie or dG plots)
%StepET_R - Cell array {number of steps(int), ET_R(vec), Start times(vec)}.
%           Alternately, a string or integer containing the number of steps
%           maybe be used. This defaults step time to every 100 units and 
%           step size to ET_R/number of steps. (Required for Rolie or dG 
%           plots)
%cnu - Value of cnu for the polymer flow. (Required for Rolie or dG plots)
%RescV - 1x3 Logical vector denoting which functionals you wish to rescale.
%TransMethod - 


%Check input format for Functionals.
if size(Functionals,1)~=1 && size(Functionals,1)~=3
    error(['Incorrect format for Functionals input, please use either a'...
        ,' cell row vector containing functional names or a cell matrix'...
        ,' containing function handles as in the Functional set output '...
        ,'from GLaMMFunctionalLinC.'])
end

%Define defaults.
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
    Plotf=0;
else
    PlotS=PlotOpt(1);
    PlotdG=PlotOpt(2);
    PlotI=PlotOpt(3);
    Plotf=PlotOpt(4);
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
%Restrcture f and G_i for interpolation. (significant speedup)
%Restructure Gi_Tilde.
%Create vectors for G1_Tilde, G2_Tilde and G3_Tilde.
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
%Old Format.
% for Samp=1:SampNo
%     fTVec={cat(3,fTVec{1,1},f_Tilde{Samp}{1,1}),...
%         cat(3,fTVec{1,2},f_Tilde{Samp}{1,2});...
%         cat(3,fTVec{2,1},f_Tilde{Samp}{2,1}),...
%         cat(3,fTVec{2,2},f_Tilde{Samp}{2,2})};
% end
%End Restructure section.
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
figure(1)
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
    dt=(Time(2:end)-Time(1:(end-1)));
    %Note this dG is good enough for a plot, but is NOT the GLaMM dG/dt.
    subplot(2,3,4)
    plot(Time(2:end),(G1hat(2:end)-G1hat(1:(end-1)))./dt,...
        'color',Blue)
    xlabel('Time')
    ylabel('d\Gamma_1/dt')
    subplot(2,3,5)
    plot(Time(2:end),(G2hat(2:end)-G2hat(1:(end-1)))./dt,...
        'color',Blue)
    xlabel('Time')
    ylabel('d\Gamma_2/dt')
    subplot(2,3,6)
    plot(Time(2:end),(G3hat(2:end)-G3hat(1:(end-1)))./dt,...
        'color',Blue)
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
    fullscreen;
    drawnow;
    pause
    %Write to .gif file. (Optional)
%     [imind,cm]=rgb2ind(frame2im(getframe(gcf)),256);
%     imwrite(imind,cm,'Newest.gif');
end
%%
%Pred f plot.
%Setup time vector.
RedTime=[0,Time];
RedTlen=numel(RedTime);
%Preallocate fPred, each GhPred, and each dGhPred. (Note the index 1 is t=0)
fPred={zeros(N,N,RedTlen),zeros(N,N,RedTlen);...
    zeros(N,N,RedTlen),zeros(N,N,RedTlen)};
NoPts=zeros(1,RedTlen);

%Caluclate G1, G2 and G3 predictions with orientation.
%1) Interpolate fhat (no orientation) from functionals G1, G2 and G3.
[fhatPred,~,NoPts(2:end)]=Functional3Interpolator(fTVec,G1TVec,...
    G2TVec,G3TVec,G1hat,G2hat,G3hat,N,PlotI,CoFunc,InterSc);
if exist('RescV','var') && isempty(RescV)==0
    %Preallocate.
    RescF=Functionals;
    %Replace 'fixed scale' functionals to empty elements.
    for i=1:numel(Functionals)
        if RescV(i)==0
            RescF{i}=[];
        elseif RescV(i)~=1
            error('Incorrect format for RescV.')
        end
    end
    %Find current values for each functional.
    [~,~,GiIn]=FunctionalEval({Time,fhatPred},N,Z,Functionals(1,:));
    %Rescale chosen functionals.
    fhatPred = FRescale(fhatPred,GiIn,{G1hat,G2hat,G3hat},N,Z,RescF,TransMethod);
end
%2) Find eigenmatrix V.
VTemp=Seigen(S);
%3) Apply orientation to calculate f from fhat.
fTemp=frotate(fhatPred,VTemp);

%Calculate the value of each functionla at time 0.
[~,~,GkZero]=FunctionalEval({0,feq},N,Z,Functionals(1,:));
%Predict functionals in principle frame.
[~,~,GkhPredTemp]=FunctionalEval({Time,fhatPred},N,Z,Functionals(1,:));
G1hPred=[GkZero{1},GkhPredTemp{1}];
G2hPred=[GkZero{2},GkhPredTemp{2}];
G3hPred=[GkZero{3},GkhPredTemp{3}];
%Record fPred.
fPred{1,1}(:,:,1)=feq{1,1};
fPred{1,2}(:,:,1)=feq{1,2};
fPred{2,1}(:,:,1)=feq{2,1};
fPred{2,2}(:,:,1)=feq{2,2};
fPred{1,1}(:,:,2:end)=fTemp{1,1};
fPred{1,2}(:,:,2:end)=fTemp{1,2};
fPred{2,1}(:,:,2:end)=fTemp{2,1};
fPred{2,2}(:,:,2:end)=fTemp{2,2};
%Calculate SPred.
SPred=SFunc(fPred);

%Cleanup for Illustrative plots.
if PlotI
    hold off
    figure(1)
end
%Locate extrapolated data.
if regexpi(CoFunc,'Oct')
    Polate=find(NoPts==8,1,'last');
    if isempty(Polate)
        Polate=1;
    end
else
    Polate=RedTlen;
end
%Plot using predicted Gamma.
subplot(1+PlotdG,3,1)
hold on
if Plotf
    %Setup functionals.
    G1Func=GLaMMFunctional(Functionals{1,1},Z,N);
    G2Func=GLaMMFunctional(Functionals{1,2},Z,N);
    G3Func=GLaMMFunctional(Functionals{1,3},Z,N);
    %Plot each functional.
    G1Pred=G1Func(fPred);
    plot(RedTime(1:Polate),G1Pred(1:Polate),...
        'color',Red,...
        'displayname','$\Gamma_i$ from model predictions.')
    if Polate~=RedTlen
        plot(RedTime(Polate:(end-1)),G1Pred(Polate:(end-1)),...
            'linestyle','--','color',Red,...
            'displayname','$\Gamma_i$ from model predictions. (extrapolated)')
    end
end
plot(RedTime(1:Polate),G1hPred(1:Polate),...
    'color',Orange,...
    'displayname','$\dot{\Gamma}_i$ from model predictions.')
if Polate~=RedTlen
    plot(RedTime(Polate:end),G1hPred(Polate:end),...
        'linestyle','--','color',Orange,...
        'displayname','$\dot{\Gamma}_i$ from model predictions. (extrapolated)')
end
hold off
subplot(1+PlotdG,3,2)
hold on
if Plotf
    G2Pred=G2Func(fPred);
    plot(RedTime(1:Polate),G2Pred(1:Polate),'color',Red)
    if Polate~=RedTlen
        plot(RedTime(Polate:(end-1)),G2Pred(Polate:(end-1)),...
            'linestyle','--','color',Red)
    end
end
plot(RedTime(1:Polate),G2hPred(1:Polate),'color',Orange)
if Polate~=RedTlen
    plot(RedTime(Polate:end),G2hPred(Polate:end),'linestyle','--',...
        'color',Orange)
end
hold off
subplot(1+PlotdG,3,3)
hold on
if Plotf
    G3Pred=G3Func(fPred);
    plot(RedTime(1:Polate),G3Pred(1:Polate),'color',Red)
    if Polate~=RedTlen
        plot(RedTime(Polate:(end-1)),G3Pred(Polate:(end-1)),...
            'linestyle','--','color',Red)
    end
end
plot(RedTime(1:Polate),G3hPred(1:Polate),'color',Orange)
if Polate~=RedTlen
    plot(RedTime(Polate:end),G3hPred(Polate:end),'linestyle','--',...
        'color',Orange)
end
%Plot derivatives if requested.
if PlotdG
    %Preallocate. (Need at least indentity for dS calc.)
    FuncNo=3;
    %----------------------------------------------------------------------
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
    %Speedup by reducing density.
    DenReduce=10;
    DenInd=DenReduce:DenReduce:(RedTlen-1);
    TimeD=RedTime(DenInd+1);
    STempD=SPred(:,:,DenInd+1);
    VTempD=VTemp(:,:,DenInd);
    fTempD=cellfun(@(x)x(:,:,DenInd),fTemp,'UniformOutput',false);
    fhTempD=cellfun(@(x)x(:,:,DenInd),fhatPred,'UniformOutput',false);
    %Calculate derivative at t=0.
    warning('off','dV:equil')
    dGZero=dGdt(Functionals,feq,0,N,Z,FlowType,StepET_R,cnu);
    warning('on','dV:equil')
    %4) Find the derivative of f using the GLaMM model on f.
    [dfTemp,~]=Multifderivs(fTempD,TimeD,N,Z,FlowType,StepET_R,cnu,pqmatrix(N));
    %5) Find the derivative of S using its formula on f and df.
    dS=dSFunc(fTempD,dfTemp);
    %Debug: dS=(S(:,:,DenReduce:DenReduce:RedTlen)-S(:,:,(DenReduce:DenReduce:RedTlen)-1))./permute(dt((DenReduce:DenReduce:RedTlen)-1),[1,3,2]);
    %6) Find the derivative of V using its formula on S and dS.
    dVTemp=MultiVDerivs(STempD,dS);
    %7) Find the derivative of fhat using its formula on f,df,S and dS.
    dfhatTemp=MultifhatDerivs(fTempD,dfTemp,VTempD,dVTemp);
    %8) Find the derivatives of the functionals using fhat and dfhat.
    dG1hPred=dG1Func(fhTempD,dfhatTemp);
    dG2hPred=dG2Func(fhTempD,dfhatTemp);
    dG3hPred=dG3Func(fhTempD,dfhatTemp);
    %----------------------------------------------------------------------
    %Plot dGdt (for G1, G2 and G3).
    hold off
    subplot(2,3,4)
    hold on
    plot([0,TimeD],[dGZero{1},dG1hPred],'color',Orange)
    hold off
    subplot(2,3,5)
    hold on
    plot([0,TimeD],[dGZero{2},dG2hPred],'color',Orange)
    hold off
    subplot(2,3,6)
    hold on
    plot([0,TimeD],[dGZero{3},dG3hPred],'color',Orange)
% DenReduce=100;
% TimeD=RedTime((1+DenReduce):DenReduce:RedTlen);
% STempD=SPred(:,:,(1+DenReduce):DenReduce:RedTlen);
% VTempD=VTemp(:,:,DenReduce:DenReduce:(RedTlen-1));
% fTempD=cellfun(@(x)x(:,:,DenReduce:DenReduce:(RedTlen-1)),fTemp,...
% 'UniformOutput',false);
% fhTempD=cellfun(@(x)x(:,:,DenReduce:DenReduce:(RedTlen-1)),fhatPred,...
% 'UniformOutput',false);
% 
% G1D=G1hat(DenReduce:DenReduce:(RedTlen-1));;
% G2D=G2hat(DenReduce:DenReduce:(RedTlen-1));;
% G3D=G3hat(DenReduce:DenReduce:(RedTlen-1));;
% [fhat_pred,df_pred,dGi_pred]=Functional3InterpolatordMulti(fTVec,G1TVec,G2TVec,G3TVec,TimeD,G1D,G2D,G3D,STempD,N,Z,Functionals,FlowType,StepET_R,cnu,0,CoFunc,InterSc);
% plot([0,TimeD],[dGZero{3},dGi_pred{3}],'color',Red)
end
hold off
%%
% Rolie-Poly plot.
if exist('Rolie','var') && isempty(Rolie)==0
    %Warning for yy.
    if regexpi(FlowType,'yy\s*u\w*\s*e\w*')
        error('S,dS and df formulae assume yy=zz.')
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
    
    %Define parameters.
    Taue=Rolie(1);
    beta=Rolie(2);
    %Preallocate.
    RTime=cell(1,StepET_R{1});
    Kappa=cell(1,StepET_R{1});
    RolieNSD=cell(1,StepET_R{1});
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
    end
    %Plot Rolie-Poly model.
    subplot(1+PlotdG,2,1)
    hold on
    plot(cell2mat(RTime),cell2mat(RolieNSD),'color',Purple,...
        'displayname','NSD from Rolie-Poly model')
    hold off
end
%%
% Sxx and Sxy plots.
if PlotS
    figure
    %Define variables.
    Sxx=permute(S(1,1,:),[3,1,2]);
    Sxy=permute(S(1,2,:),[3,1,2]);
    SxxPred=permute(SPred(1,1,:),[3,1,2]);
    SxyPred=permute(SPred(1,2,:),[3,1,2]);
    %Plot true and predicted S values.
    subplot(1+PlotdG,2,1)
    %True Sxx
    plot(Time,Sxx,'color',Green,...
        'displayname','$S$ from GLaMM predictions (dt=1/90)')
    hold on
    if Plotf
        %Only for testing.
        [~,Shat]=fPrinciple(fhat);
        plot(Time,permute(Shat(1,1,:),[3,1,2]),'color',Blue,...
            'displayname','$\dot{S}$ from GLaMM predictions (dt=1/90)')
    end
    %Pred Sxx
    plot(RedTime(1:Polate),SxxPred(1:Polate),'color',Red,...
        'displayname','$S$ from model predictions.')
    if Polate~=RedTlen
        plot(RedTime(Polate:end),SxxPred(Polate:end),...
            'linestyle','--','color',Red,...
            'displayname','$S$ from model predictions.')
    end
    hold off
    xlabel('Time')
    ylabel('S_{xx}')
    legend('interpreter','latex','location','southeast')
    subplot(1+PlotdG,2,2)
    %True Sxy
    plot(Time,Sxy,'color',Green)
    hold on
    if Plotf
        plot(Time,permute(Shat(1,2,:),[3,1,2]),'color',Blue)
    end
    %Pred Sxy
    plot(RedTime(1:Polate),SxyPred(1:Polate),'color',Red)
    if Polate~=RedTlen
        plot(RedTime(Polate:end),SxyPred(Polate:end),...
            'linestyle','--','color',Red)
    end
    hold off
    xlabel('Time')
    ylabel('S_{xy}')
    if PlotdG
        %dSxx
        subplot(2,2,3)
        hold on
        %True dSxx
        %Note this dS is good enough for a plot, but is NOT the GLaMM dS/dt.
        plot(Time(2:end),...
            permute(S(1,1,2:end)-S(1,1,1:(end-1)),[3,1,2])./dt,...
            'color',Green)
        %Pred dSxx
        plot(TimeD,permute(dS(1,1,:),[3,1,2]),'color',Red)
        hold off
        %Label axis
        xlabel('Time')
        ylabel('dS_{xx}/dt')
        %dSxy
        subplot(2,2,4)
        hold on
        %True dSxy
        plot(Time(2:end),...
            permute(S(1,2,2:end)-S(1,2,1:(end-1)),[3,1,2])./dt,...
            'color',Green)
        %Pred dSxy
        plot(TimeD,permute(dS(1,2,:),[3,1,2]),'color',Red)
        hold off
        %Label axis
        xlabel('Time')
        ylabel('dS_{xy}/dt')
    end
end
%%
% Evaluate metric.
if exist('Metric','var') && isempty(Metric)==0
    error('Using older Metric code without orentation or G3.')
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
            %Calculate predicted dG if nessecary.
            if exist('dGPred','var')==0
                dGZero=dGdt(Functionals(:,str2double(Metric{1,m}(3))),...
                    feq,0,N,Z,FlowType,StepET_R,cnu);
                dGPred=dGdt(Functionals(:,str2double(Metric{1,m}(3))),...
                    fhatPred,RedTime,N,Z,FlowType,StepET_R,cnu);
                dGPredRed=[dGZero{1},dGPred{1}(1:(end-1))];
            else
                dGPredRed=[dGZero{1},...
                    dGPred{str2double(Metric{1,m}(3))}(1:(end-1))];
            end
            %Extract true dG from GLaMM data.
            dGTrue=dGdt(Functionals(:,str2double(Metric{1,m}(3))),...
                fhat,Time,N,Z,FlowType,StepET_R,cnu);
            %Reformat true dG data and include t=0.
            dGTrueRed=[dGZero{1},dGTrue{1}(1:(end-1))];
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
            TScore{m}=MetricFunc(fhatPred,fhat); %fhatPred may need changing to match time.
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

%Examples

%Principle3({{Time8SB,fhat8SB},{Time10SB,fhat10SB}},75,25,{{Time9SB,f9SB}},75,25,{'NSD','NSDCs','PrinOrienVar'});

%Traindata=TestData('Train');
%[Shearhatdata,ETRS]=TestData('Shearhat');
%Principle3([Traindata,Shearhatdata(ETRS~=9)],75,25,Sheardata(ETRS==9),75,25,{'NSD','NSDCs','PrinOrienVar'},'Oct');
%Principle3([Traindata,Shearhatdata(ETRS~=9)],75,25,Sheardata(ETRS==9),75,25,{'NSD','NSDCs','PrinOrienVar'},'Near');
%Principle3([Traindata,Shearhatdata(ETRS~=9)],75,25,Sheardata(ETRS==9),75,25,{'NSD','NSDCs','PrinOrienVar'},'Rad');

%G3 Rescaling needs fixing first!
% [Traindata,ETRT]=TestData('Train');
% [Sheardata,ETRS]=TestData('Shear');
% Shearhatdata=TestData('Shearhat');
% Principle3(Traindata(ETRT~=1.9),75,25,Traindata(ETRT==1.9),75,25,{'NSD','NSDCs','PrinOrienVar'},'Oct','Rel',[0,1,0,0],[],'xxUE',{1,1.9,0},0);
% Principle3(Traindata(ETRT~=1.9),75,25,Traindata(ETRT==1.9),75,25,{'NSD','NSDCs','PrinOrienVar'},'Oct','Rel',[0,1,0,0],[],'xxUE',{1,1.9,0},0,[1,1,0],'Multi');
% Principle3([Traindata,Shearhatdata(or(ETRS<=8,ETRS>=10))],75,25,Sheardata(ETRS==9),75,25,{'NSD','NSDCs','PrinOrienVar'},'Oct','Rel',[1,1,0,0],[],'xySH',{1,9,0},0.1);
% Principle3([Traindata,Shearhatdata(or(ETRS<=8,ETRS>=10))],75,25,Sheardata(ETRS==9),75,25,{'NSD','NSDCs','PrinOrienVar'},'Oct','Rel',[1,1,0,0],[],'xySH',{1,9,0},0.1,[1,1,0],'Multi');
% Principle3([Traindata,Shearhatdata(or(ETRS<=8,ETRS>=10))],75,25,Sheardata(ETRS==9),75,25,{'NSD','NSDCs','PrinOrienVar'},'Rad','Rel',[1,1,0,0],[],'xySH',{1,9,0},0.1);
% Principle3([Traindata,Shearhatdata(or(ETRS<=8,ETRS>=10))],75,25,Sheardata(ETRS==9),75,25,{'NSD','NSDCs','PrinOrienVar'},'Rad','Rel',[1,1,0,0],[],'xySH',{1,9,0},0.1,[1,1,0],'Multi');