function [Time,f] = FijpqExtract(File,N,Z,OverrideType,Override,Functionals)
%FIJPQEXTRACT is an updtated version of FpqExtract.
%
%Inputs:
%File - String containing name of file to open.
%N - N value from conditions.
%Z - Z value from conditions. (Only required if an OverrideType is defined.)
%OverrideType - Optional parameters to set functional you wish to limit.
%Override - Optional Vector of min and max cutoff values based on OverrideType.
%Functionals - Optional cell array containing names of functionals for use
%              with 'GRepeat' OverrideType.
%
%Outputs:
%Time - Vector of times at which f was evaluated.
%f - Cell array contaning 4 f_ij(s,s',Timeindex) type arrays.

%DOES NOT WORK FOR NON y,z SYMMETRIC FLOWS!
if ischar(File)==1
    M=dlmread(['GLaMM Results\',File]);
else
    M=File;
end
RecordNo=size(M,1)/(1+N^2);
RecordVec=1:RecordNo;
%Predifine size.
Time=zeros(1,RecordNo);
Fxx=zeros(N,N,RecordNo);
Fxy=zeros(N,N,RecordNo);
Fyy=zeros(N,N,RecordNo);
%Define functions that calculate the positions of Timestamps and f(s,:).
Tind=@(i) RecordVec(i)+(RecordVec(i)-1)*N^2;
FSind=@(i,s) RecordVec(i)+(RecordVec(i)-1)*N^2+(s-1)*N+(1:N);
for i=1:RecordNo
    Time(i)=M(Tind(i),1);
    for s=1:N
        Fxx(s,:,i)=M(FSind(i,s),1);
        Fxy(s,:,i)=M(FSind(i,s),2);
        Fyy(s,:,i)=M(FSind(i,s),3);
        %Colon is S'
        %Note that S,S' are actually s*Z/N where s=index(S).
    end
end
f={Fxx,Fxy;Fxy,Fyy};
if exist('OverrideType','var') && isempty(OverrideType)==0
    if regexpi(OverrideType,'GRepeat')
        %Preallocate.
        NewRange=zeros(1,RecordNo);
        %Evaluate functionals.
        [~,~,Gi]=FunctionalEval({Time,f},N,Z,Functionals,[],[]);
        %Set starting value.
        NewRange(1)=1;
        PrevGi=Gi(1);
        %Option for removal based on relative distance.
        if regexpi(OverrideType,'GRepeatRel')
            MinDist=(Override(1))*sqrt(sum(PrevGi.^2));
        else
            MinDist=Override(1);
        end
        %Trim data by only recording when a new data point is sufficently 
        %far from the previously recorded data point.
        for i=1:RecordNo
            if sqrt(sum((PrevGi-Gi(i)).^2))>=MinDist
                NewRange(i)=1;
                PrevGi=Gi(i);
                if regexpi(OverrideType,'GRepeatRel')
                    MinDist=(Override(1))*sqrt(sum(PrevGi.^2));
                end
            end
        end
    else
        if regexpi(OverrideType,'Density')
            %Density reduction.
            TestVec=zeros(1,RecordNo);
            for j=Override:Override:RecordNo
                TestVec(j)=1;
            end
            Override=[0.5,1.5];
        elseif regexpi(OverrideType,'Time')
            %Time type cutoff.
            TestVec=Time;
        elseif regexpi(OverrideType,'Steady')
            %Steady state cutoff.
            %Preallocate
            TestVec=ones(1,RecordNo);
            %Create f shifted one timsstep.
            fShift={cat(3,zeros(N),Fxx(:,:,1:(RecordNo-1))),...
                cat(3,zeros(N),Fxy(:,:,1:(RecordNo-1)));...
                cat(3,zeros(N),Fxy(:,:,1:(RecordNo-1))),...
                cat(3,zeros(N),Fyy(:,:,1:(RecordNo-1)))};
            %Create default tolerance.
            if exist('Override','var')==0 || isempty(Override)
                Override=[1e-3,inf];
            end
            %Evaluate fRMSE between each step.
            MFunc=GLaMMMetric('fRMSE',Z,Z);
            fRMSE=MFunc(f,fShift);
            %Evaluate average error for previous 5 timesteps.
            for j=5:RecordNo
                TestVec(j)=mean(fRMSE((j-4):j));
            end
            %Alternate test. (unused)
            %TestVec=sum(Fxx-cat(3,zeros(N),Fxx(:,:,1:(RecordNo-1))),[1,2]);
        else
            %Functional type cutoff.
            if regexpi(OverrideType,'Stretch\s*Ratio')
                %Stretch ratio type cutoff.
                %Create default stretch ratio min max.
                if exist('Override','var')==0 || isempty(Override)
                    Override=[1,6];
                elseif strcmpi(Override,'NaN')
                    %Removes elements that have an undefined stretch ratio.
                    Override=[-inf,inf];
                end
            else
                if exist('Override','var')==0 || isempty(Override)
                    error('Undefined override type.');
                end
            end
            %Calculate reduced range.
            OverrideFunc=GLaMMFunctional(OverrideType,Z);
            TestVec=OverrideFunc(f);
        end
        NewRange=find(TestVec>=Override(1) & TestVec<=Override(2));
    end
    %Code to ensure that both limits are within the range (Unused).
%     if NewRange(1)>1
%         NewRange=[NewRange(1)-1,NewRange];
%     end
%     if NewRange(end)<RecordNo
%         NewRange=[NewRange,NewRange(end)+1];
%     end
    %Setup new vectors.
    Time=Time(NewRange);
    f={Fxx(:,:,NewRange),Fxy(:,:,NewRange);Fxy(:,:,NewRange),Fyy(:,:,NewRange)};
end
end