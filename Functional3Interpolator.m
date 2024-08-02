function [f_pred,MinDis,NoPts]=Functional3Interpolator(fTVec,G1TVec,G2TVec,G3TVec,G1,G2,G3,N,Illustrative,CoFunc,InterSc,FuncOpt)
%FUNCTIONAL3INTERPOLATOR Predicts f(s,s') for a gamma triplet (G1,G2,G3) 
%with restructured dataset fTVec,G1TVec,G2TVec,G3TVec.
%
%Inputs:
%fTVec,G1TVec,G2TVec,G3TVec - Restructured training data for speed.
%G1,G2,G3 - Vector of evaluations of each functional at each timestep for 
%the test data.
%N - N value of test data.
%Illustrative - Logical variable, if true, shows points used in each 
%               interpolation.
%CoFunc - String containing covariance function to use, defaults to quad.
%InterSc -
%FuncOpt -
%
%Output:
%f_pred - Cell indexed by i,j containing 3D arrays (s,s',t) of predictions 
%for f_{ij}(s,s').
%MinDis - Distance to the closest point used in interpolation.
%NoPts - Number of points used for the interpolation. (0 is exact.)

% %Mpts test line.
% FuncOpt=[1000,1,2,1]
%Detemine Output length.
OutLength=numel(G1);
%Preallocate ouput.
NNT=zeros(N,N,OutLength);
f_pred={NNT,NNT;NNT,NNT};
MinDis=zeros(1,OutLength);
%Determine scale.
%Find max to rescale each functional.
maxG1=max(abs(G1TVec(~isinf(G1TVec))),[],'all');
maxG2=max(abs(G2TVec(~isinf(G2TVec))),[],'all');
maxG3=max(abs(G3TVec(~isinf(G3TVec))),[],'all');
if regexpi(InterSc,'Abs')
    %Set scale based on maximum value 
    %Correct if any of these values are zero. (Else distance returns NaN.)
    G1Sc=(maxG1+(maxG1==0))*ones(1,OutLength);
    G2Sc=(maxG2+(maxG2==0))*ones(1,OutLength);
    G3Sc=(maxG3+(maxG3==0))*ones(1,OutLength);
elseif regexpi(InterSc,'Rel')
    %Create vectors for scaling based on current size.
    G1Sc=max(G1,(maxG1+(maxG1==0))*(10^-10));
    G2Sc=max(G2,(maxG2+(maxG2==0))*(10^-10));
    G3Sc=max(G3,(maxG3+(maxG3==0))*(10^-10));
else
    error('Unrecognised interpolation scaling method.')
end
%Predict data.
if regexpi(CoFunc,'Oct\w*')
    NoPts=8*ones(1,OutLength);
    %Default PenaltyValue to inf.
    if exist('PenaltyValue','var')==0
        PenaltyValue=inf;
    end
    %Preallocate vectors for the chosen points.
    DisQ=zeros(1,8);
    LocQ=zeros(1,8);
    ValQ=cell(1,8);
    for i=1:OutLength
        %Find the relative gamma distance for all test data to Gi_pred.
        Dist=sqrt(((G1TVec-G1(i))/G1Sc(i)).^2+...
            ((G2TVec-G2(i))/G2Sc(i)).^2+...
            ((G3TVec-G3(i))/G3Sc(i)).^2);
        %Apply any distance transformations.
        if isempty(regexpi(CoFunc,'2','once'))==0 || isempty(regexpi(CoFunc,'Sq','once'))==0
            Dist=Dist.^2;
        elseif regexpi(CoFunc,'4')
            Dist=Dist.^4;
        elseif regexpi(CoFunc,'Exp')
            Dist=exp(Dist);
        end
        %Find quadrant for all training data at Gi_pred.
        Quadrant=1+(G1TVec<G1(i))+2*(G2TVec<G2(i))+4*(G3TVec<G3(i));
        for j=1:8
            if sum(Quadrant==j)>0 || PenaltyValue~=inf
                %Create penalty term for incorrect quadrants.
                Penalty=zeros(size(Dist));
                Penalty(Quadrant~=j)=PenaltyValue;
                %Find location and distance to closest point with penalty.
                [~,LocQ(j)]=min(Dist+Penalty);
                DisQ(j)=Dist(LocQ(j));
                %Find the value at that point.
                ValQ{j}={fTVec{1,1}(:,:,LocQ(j)),fTVec{1,2}(:,:,LocQ(j));...
                    fTVec{2,1}(:,:,LocQ(j)),fTVec{2,2}(:,:,LocQ(j))};
            else
                %Outside of data range, calculate based on remaining points.
                LocQ(j)=nan;
                DisQ(j)=nan;
            end
        end
        %Calculate output using weightings.
        if DisQ(1)==0
            for k=1:4 %xx,xy,yx,yy.
                %Error catching for exact matches.
                f_pred{k}(:,:,i)=ValQ{1}{k};
            end
            NoPts(i)=0;
        else
            for k=1:4 %xx,xy,yx,yy.
                for j=1:8 %Quadrants.
                    if ~isnan(DisQ(j))
                        f_pred{k}(:,:,i)=f_pred{k}(:,:,i)+...
                            (ValQ{j}{k}./DisQ(j))/sum(1./DisQ,'omitnan');
                    end
                end
            end
            %Record minimum distance and number of points used.
            MinDis(i)=min(DisQ);
            NoPts(i)=sum(1-isnan(DisQ));
        end
        %If you wish to extract points used for illustrative purposes.
        if i<=Illustrative && NoPts(i)~=0
            title(['Time index=',num2str(i),', NoPts=',num2str(NoPts(i))])
            delay=0.1;
            if NoPts(i)==0
                %Plot single point.
                InTether=plot3(G1(i),G2(i),G3(i),'color',[1, 11/17, 0],...
                    'marker','o');
            else
                %Find location indicies by removing quadrants with no results.
                Lind=LocQ(~isnan(DisQ));
                InTether=plot3([G1(i)*ones(1,NoPts(i));G1TVec(Lind)],...
                    [G2(i)*ones(1,NoPts(i));G2TVec(Lind)],...
                    [G3(i)*ones(1,NoPts(i));G3TVec(Lind)],...
                    'color',[1, 11/17, 0],'linestyle','--','marker','o',...
                    'markerindices',2);
            end
            %Write to .gif file. (Optional)
%             drawnow;
%             [imind,cm]=rgb2ind(frame2im(getframe(gcf)),256);
%             imwrite(imind,cm,'Newest.gif','gif','writemode','append');
            if delay==inf
                pause
            else
                pause(delay)
            end
            if delay~=0
                delete(InTether)
            end
        end
    end
elseif regexpi(CoFunc,'Near\w*')
    %Setup defaults
    if exist('FuncOpt','var')==0 || isempty(FuncOpt)
        MPts=1000;
    else
        MPts=FuncOpt(1);
    end
    NoPts=MPts*ones(1,OutLength);
    for i=1:OutLength
        %Find the relative gamma distance for all test data to Gi_pred.
        Dist=sqrt(((G1TVec-G1(i))/G1Sc(i)).^2+...
            ((G2TVec-G2(i))/G2Sc(i)).^2+...
            ((G3TVec-G3(i))/G3Sc(i)).^2);
        %Apply any distance transformations.
        if isempty(regexpi(CoFunc,'2','once'))==0 || isempty(regexpi(CoFunc,'Sq','once'))==0
            Dist=Dist.^2;
        elseif regexpi(CoFunc,'Exp')
            Dist=exp(Dist);
        end
        %Create list sorted by distance to test point.
        [~,SortV]=sort(Dist);
        %Find location (index) of closest M points, and their distance.
        LocN=SortV(1:MPts);
        DisN=permute(Dist(LocN),[1,3,2]);
        %Find the value at each point.
        ValN={fTVec{1,1}(:,:,LocN),fTVec{1,2}(:,:,LocN);...
            fTVec{2,1}(:,:,LocN),fTVec{2,2}(:,:,LocN)};
        %Calculate output using weightings.
        if DisN(1)==0
            for k=1:4 %xx,xy,yx,yy.
                %Error catching for exact matches.
                f_pred{k}(:,:,i)=ValN{k}(:,:,1);
            end
            NoPts(i)=0;
        else
            for k=1:4 %xx,xy,yx,yy.
                f_pred{k}(:,:,i)=sum((ValN{k}./DisN)/sum(1./DisN),3);
            end
        end
        %Record minimum distance.
        MinDis(i)=DisN(1);
        %If you wish to extract points used for illustrative purposes.
        if i<=Illustrative && NoPts(i)~=0
            title(['Time index=',num2str(i),', NoPts=',num2str(NoPts(i))])
            delay=0.1;
            if NoPts(i)==0
                %Plot single point.
                InTether=plot3(G1(i),G2(i),G3(i),'color',[1, 11/17, 0],...
                    'marker','o');
            else
                InTether=plot3([G1(i)*ones(1,MPts);G1TVec(LocN)],...
                    [G2(i)*ones(1,MPts);G2TVec(LocN)],...
                    [G3(i)*ones(1,MPts);G3TVec(LocN)],...
                    'color',[1, 11/17, 0],'linestyle','--','marker','o',...
                    'markerindices',2);
            end
            %Write to .gif file. (Optional)
%             drawnow;
%             [imind,cm]=rgb2ind(frame2im(getframe(gcf)),256);
%             imwrite(imind,cm,'Newest.gif','gif','writemode','append');
            if delay==inf
                pause
            else
                pause(delay)
            end
            if delay~=0
                delete(InTether)
            end
        end
    end
elseif regexpi(CoFunc,'Rad\w*')
    %Setup defaults
    if exist('FuncOpt','var')==0 || isempty(FuncOpt)
        MPts=100;
        if regexpi(InterSc,'Abs')
            Rad=0.1;
        elseif regexpi(InterSc,'Rel')
            Rad=0.5;
        end
    else
        MPts=FuncOpt(1);
        Rad=FuncOpt(2);
    end
    NoPts=MPts*ones(1,OutLength);
    %Create sphere for plots.
    [X,Y,Z]=sphere;
    for i=1:OutLength
        %Find the relative gamma distance for all test data to Gi_pred.
        Dist=sqrt(((G1TVec-G1(i))/G1Sc(i)).^2+...
            ((G2TVec-G2(i))/G2Sc(i)).^2+...
            ((G3TVec-G3(i))/G3Sc(i)).^2);
        %Apply any distance transformations.
        if isempty(regexpi(CoFunc,'2','once'))==0 || isempty(regexpi(CoFunc,'Sq','once'))==0
            Dist=Dist.^2;
        elseif regexpi(CoFunc,'Exp')
            Dist=exp(Dist);
        end
        %Find location of all points within a radius around test point.
        LocN=find(Dist<(Rad^2));
        %Determine number of points.
        NoPts(i)=numel(LocN);
        %Expand if no points found.
        while NoPts(i)==0
            Rad=2*Rad;
            warning(['Radius of sphere too small to find any relevant ',...
                'points, size has been doubled to ',num2str(Rad),'.'])
            LocN=find(Dist<(Rad^2));
            NoPts(i)=numel(LocN);
        end
        %Remove points if greater than max number of points.
        if NoPts(i)>MPts
            [~,SortV]=sort(Dist);
            LocN=SortV(1:MPts);
            NoPts(i)=MPts;
        end
        %Record distance for each point.
        DisN=permute(Dist(LocN),[1,3,2]);
        %Find the value at each point.
        ValN={fTVec{1,1}(:,:,LocN),fTVec{1,2}(:,:,LocN);...
            fTVec{2,1}(:,:,LocN),fTVec{2,2}(:,:,LocN)};
        %Calculate output using weightings.
        if sum(DisN==0)>0
            for k=1:4 %xx,xy,yx,yy.
                %Error catching for exact matches.
                f_pred{k}(:,:,i)=ValN{k}(:,:,DisN==0);
            end
            NoPts(i)=0;
        else
            for k=1:4 %xx,xy,yx,yy.
                f_pred{k}(:,:,i)=sum((ValN{k}./DisN)/sum(1./DisN),3);
            end
        end
        %Record minimum distance.
        MinDis(i)=min(DisN);
        %If you wish to extract points used for illustrative purposes.
        if i<=Illustrative && NoPts(i)~=0
            title(['Time index=',num2str(i),', NoPts=',num2str(NoPts(i))])
            delay=0.1;
            if NoPts(i)==0
                %Plot single point.
                InTether=plot3(G1(i),G2(i),G3(i),'color',[1, 11/17, 0],...
                    'marker','o');
            else
                %Plot tethers.
                InTether=plot3([G1(i)*ones(1,NoPts(i));G1TVec(LocN)],...
                    [G2(i)*ones(1,NoPts(i));G2TVec(LocN)],...
                    [G3(i)*ones(1,NoPts(i));G3TVec(LocN)],...
                    'color',[1, 11/17, 0],'linestyle','--','marker','o',...
                    'markerindices',2);
                %Rescale and plot sphere.
                Xi=X*Rad*G1Sc(i);
                Yi=Y*Rad*G2Sc(i);
                Zi=Z*Rad*G3Sc(i);
                Sphere=surf(Xi+G1(i),Yi+G2(i),Zi+G3(i),'edgecolor','flat',...
                    'edgealpha',0.1,'facealpha',0.1);
            end
            %Write to .gif file. (Optional)
%             drawnow;
%             [imind,cm]=rgb2ind(frame2im(getframe(gcf)),256);
%             imwrite(imind,cm,'Newest.gif','gif','writemode','append');
            if delay==inf
                pause
            else
                pause(delay)
            end
            if delay~=0
                delete(InTether)
                delete(Sphere)
            end
        end
    end
elseif regexpi(CoFunc,'OWR\w*')
    %Setup defaults
    if exist('FuncOpt','var')==0 || isempty(FuncOpt)
        MPts=100;
        if regexpi(InterSc,'Abs')
            Rad=0.1;
        elseif regexpi(InterSc,'Rel')
            Rad=0.5;
        end
    else
        MPts=FuncOpt(1);
        Rad=FuncOpt(2);
    end
    NoPts=MPts*ones(1,OutLength);
    %Create sphere for plots.
    [X,Y,Z]=sphere;
    %Preallocate vectors for the chosen points.
    for i=1:OutLength
        %Find the relative gamma distance for all test data to Gi_pred.
        Dist=sqrt(((G1TVec-G1(i))/G1Sc(i)).^2+...
            ((G2TVec-G2(i))/G2Sc(i)).^2+...
            ((G3TVec-G3(i))/G3Sc(i)).^2);
        %Apply any distance transformations.
        if isempty(regexpi(CoFunc,'2','once'))==0 || isempty(regexpi(CoFunc,'Sq','once'))==0
            Dist=Dist.^2;
        elseif regexpi(CoFunc,'4')
            Dist=Dist.^4;
        elseif regexpi(CoFunc,'Exp')
            Dist=exp(Dist);
        end
        %Find quadrant for all training data at Gi_pred.
        Quadrant=1+(G1TVec<G1(i))+2*(G2TVec<G2(i))+4*(G3TVec<G3(i));
        %Find location of all points within a radius around test point.
        LocN=find(Dist<(Rad^2));
        %Determine number of points.
        NoPts(i)=numel(LocN);
        %Expand if no points found.
        while NoPts(i)==0
            Rad=2*Rad;
            warning(['Radius of sphere too small to find any relevant ',...
                'points, size has been doubled to ',num2str(Rad),'.'])
            LocN=find(Dist<(Rad^2));
            NoPts(i)=numel(LocN);
        end
        %Remove points if greater than max number of points.
        if NoPts(i)>MPts
            [~,SortV]=sort(Dist);
            LocN=SortV(1:MPts);
            NoPts(i)=MPts;
        end
        %Determine Octant for each point.
        OctN=Quadrant(LocN);
        %Determine number of points in each Octant within radius.
        [ON,OV]=groupcounts(OctN');
        NoOctant=zeros(1,8);
        NoOctant(OV)=ON;
        %Record distance for each point.
        DisN=permute(Dist(LocN),[1,3,2]);
        %Calculate distance weighted by number of points.
        WDisN=DisN.*permute(NoOctant(OctN),[1,3,2]);
        %Find the value at each point.
        ValN={fTVec{1,1}(:,:,LocN),fTVec{1,2}(:,:,LocN);...
            fTVec{2,1}(:,:,LocN),fTVec{2,2}(:,:,LocN)};
        %Calculate output using weightings.
        if sum(DisN==0)>0
            for k=1:4 %xx,xy,yx,yy.
                %Error catching for exact matches.
                f_pred{k}(:,:,i)=ValN{k}(:,:,DisN==0);
            end
            NoPts(i)=0;
        else
            for k=1:4 %xx,xy,yx,yy.
                f_pred{k}(:,:,i)=sum((ValN{k}./WDisN)/sum(1./WDisN),3);
            end
        end
        %Record minimum distance.
        MinDis(i)=min(DisN);
        %If you wish to extract points used for illustrative purposes.
        if i<=Illustrative && NoPts(i)~=0
            title(['Time index=',num2str(i),', NoPts=',num2str(NoPts(i))])
            delay=0.1;
            if NoPts(i)==0
                %Plot single point.
                InTether=plot3(G1(i),G2(i),G3(i),'color',[1, 11/17, 0],...
                    'marker','o');
            else
                %Plot tethers.
                InTether=plot3([G1(i)*ones(1,NoPts(i));G1TVec(LocN)],...
                    [G2(i)*ones(1,NoPts(i));G2TVec(LocN)],...
                    [G3(i)*ones(1,NoPts(i));G3TVec(LocN)],...
                    'color',[1, 11/17, 0],'linestyle','--','marker','o',...
                    'markerindices',2);
                %Rescale and plot sphere.
                Xi=X*Rad*G1Sc(i);
                Yi=Y*Rad*G2Sc(i);
                Zi=Z*Rad*G3Sc(i);
                Sphere=surf(Xi+G1(i),Yi+G2(i),Zi+G3(i),'edgecolor','flat',...
                    'edgealpha',0.1,'facealpha',0.1);
            end
            if delay==inf
                pause
            else
                pause(delay)
            end
            if delay~=0
                delete(InTether)
                delete(Sphere)
            end
        end
    end
elseif regexpi(CoFunc,'ORAO\w*')
    %Setup defaults
    if exist('FuncOpt','var')==0 || isempty(FuncOpt)
        MPts=100;
        if regexpi(InterSc,'Abs')
            Rad=0.1;
        elseif regexpi(InterSc,'Rel')
            Rad=0.5;
        end
        if isempty(regexpi(CoFunc,'2','once'))==0 || isempty(regexpi(CoFunc,'Sq','once'))==0
            Alph=2;
        elseif regexpi(CoFunc,'4')
            Alph=4;
        elseif regexpi(CoFunc,'Exp')
            Alph=1;
        else
            Alph=1;
        end
        if regexpi(CoFunc,'ORAOO\w*')
            OctOver=1;
        else
            OctOver=0;
        end
    else
        MPts=FuncOpt(1);
        Rad=FuncOpt(2);
        Alph=FuncOpt(3);
        OctOver=FuncOpt(4);
    end
    NoPts=MPts*ones(1,OutLength);
    %Create sphere for plots.
    [X,Y,Z]=sphere;
    %Preallocate vectors for the chosen points.
    for i=1:OutLength
        %Find the relative gamma distance for all test data to Gi_pred.
        Dist=sqrt(((G1TVec-G1(i))/G1Sc(i)).^2+...
            ((G2TVec-G2(i))/G2Sc(i)).^2+...
            ((G3TVec-G3(i))/G3Sc(i)).^2);
        %Apply any distance transformations.
        if isempty(regexpi(CoFunc,'2','once'))==0 || isempty(regexpi(CoFunc,'Sq','once'))==0
            Dist=Dist.^2;
        elseif regexpi(CoFunc,'4')
            Dist=Dist.^4;
        elseif regexpi(CoFunc,'Exp')
            Dist=exp(Dist);
        end
        %Find Octant for all training data at Gi_pred.
        Octant=1+(G1TVec<G1(i))+2*(G2TVec<G2(i))+4*(G3TVec<G3(i));
        %Find location of all points within a radius around test point.
        LocN=find(Dist<(Rad^2));
        %Determine number of points.
        NoPts(i)=numel(LocN);
        %Expand if no points found.
        while NoPts(i)==0
            Rad=2*Rad;
            warning(['Radius of sphere too small to find any relevant ',...
                'points, size has been doubled to ',num2str(Rad),'.'])
            LocN=find(Dist<(Rad^2));
            NoPts(i)=numel(LocN);
        end
        %Remove points if greater than max number of points.
        if NoPts(i)>MPts
            [~,SortV]=sort(Dist);
            LocN=SortV(1:MPts);
            NoPts(i)=MPts;
        end
        %Determine Octant for each point.
        OctN=Octant(LocN);
        %Determine number of points in each Octant within radius.
        [ON,OV]=groupcounts(OctN');
        NoOctant=zeros(1,8);
        NoOctant(OV)=ON;
        %Override to guarantee at least one point from each octant.
        if OctOver
            %Find which octants are missing in radius, but exist in set.
            [~,TOV]=groupcounts(Octant');
            MisO=TOV(~ismember(TOV,OV));
            if ~isempty(MisO)
                %Find index of closest training point in desired octant.
                MisLoc=zeros(1,numel(MisO));
                for m=1:(numel(MisO))
                    LocCurO=find(Octant==MisO(m));
                    [~,MinLocCurO]=min(Dist(LocCurO));
                    MisLoc(m)=LocCurO(MinLocCurO);
                end
                %Append these indices onto the list of evaluated points.
                LocN=[LocN,MisLoc];
                OctN=[OctN,MisO'];
                NoPts(i)=numel(LocN);
                NoOctant(MisO)=NoOctant(MisO)+1;
            end
        end
        %Record distance for each point.
        DisN=permute(Dist(LocN),[1,3,2]);
        %Find the value at each point.
        ValN={fTVec{1,1}(:,:,LocN),fTVec{1,2}(:,:,LocN);...
            fTVec{2,1}(:,:,LocN),fTVec{2,2}(:,:,LocN)};
        %Preallocate for Octants.
        DisBar=nan(1,8);
        ValBar=cell(1,8);
        for j=1:8
            ValBar{j}=cell(2);
        end
        %Error catching for exact matches.
        if sum(DisN==0)>0
            for k=1:4 %xx,xy,yx,yy.
                f_pred{k}(:,:,i)=ValN{k}(:,:,1);
            end
        else


        %Calculate average value and distance for each Octant.
        for j=1:8 %Octants
            if NoOctant(j)~=0
                WeiOctN=DisN(OctN==j).^-Alph;
                DisBar(j)=sum((DisN(OctN==j).*WeiOctN)/sum(WeiOctN),3);
                for k=1:4 %xx,xy,yx,yy.
                    ValBar{j}{k}=sum((ValN{k}(:,:,OctN==j).*WeiOctN)/...
                        sum(WeiOctN),3);
                end
            end
        end

        %Use 'Octant' interpolator on the resultant points.
        for k=1:4 %xx,xy,yx,yy.
            for j=1:8 %Octants.
                if ~isnan(DisBar(j))
                    f_pred{k}(:,:,i)=f_pred{k}(:,:,i)+...
                        (ValBar{j}{k}./DisBar(j))/sum(1./DisBar,'omitnan');
                end
            end
        end
        %Record minimum distance and number of points used.
        MinDis(i)=min(DisBar);
        end

        %If you wish to extract points used for illustrative purposes.
        if i<=Illustrative && NoPts(i)~=0
            title(['Time index=',num2str(i),', NoPts=',num2str(NoPts(i));...
                'Number of Octants used=',num2str(numel(TOV)),...
                'Number of points outside radius=',num2str(numel(MisO))])
            delay=0.1;
            if NoPts(i)==0
                %Plot single point.
                InTether=plot3(G1(i),G2(i),G3(i),'color',[1, 11/17, 0],...
                    'marker','o');
            else
                %Plot tethers.
                InTether=plot3([G1(i)*ones(1,NoPts(i));G1TVec(LocN)],...
                    [G2(i)*ones(1,NoPts(i));G2TVec(LocN)],...
                    [G3(i)*ones(1,NoPts(i));G3TVec(LocN)],...
                    'color',[1, 11/17, 0],'linestyle','--','marker','o',...
                    'markerindices',2);
                %Rescale and plot sphere.
                Xi=X*Rad*G1Sc(i);
                Yi=Y*Rad*G2Sc(i);
                Zi=Z*Rad*G3Sc(i);
                Sphere=surf(Xi+G1(i),Yi+G2(i),Zi+G3(i),'edgecolor','flat',...
                    'edgealpha',0.1,'facealpha',0.1);
            end
            if delay==inf
                pause
            else
                pause(delay)
            end
            if delay~=0
                delete(InTether)
                delete(Sphere)
            end
        end
    end
elseif regexpi(CoFunc,'Q\w*Oct\w*')
    NoPts=8*ones(1,OutLength);
    %Default PenaltyValue to inf.
    if exist('PenaltyValue','var')==0
        PenaltyValue=inf;
    end
    %Preallocate vectors for the chosen points.
    DisQ=nan*ones(1,8);
    LocQ=zeros(1,8);
    ValQ=cell(1,8);
    for i=1:OutLength
        %Find the relative gamma distance for all test data to Gi_pred.
        Dist=sqrt(((G1TVec-G1(i))/G1Sc(i)).^2+...
            ((G2TVec-G2(i))/G2Sc(i)).^2+...
            ((G3TVec-G3(i))/G3Sc(i)).^2);
        %Apply any distance transformations.
        if isempty(regexpi(CoFunc,'2','once'))==0 || isempty(regexpi(CoFunc,'Sq','once'))==0
            Dist=Dist.^2;
        elseif regexpi(CoFunc,'Exp')
            Dist=exp(Dist);
        end
        %Revert to Quad if G3 is 0 for test point.
        if G3(i)==0
            Dist(G3TVec~=0)=inf;
        end
        Quadrant=1+(G1TVec<G1(i))+2*(G2TVec<G2(i))+4*(G3TVec<G3(i));
        for j=1:(8-(4*(G3(i)==0)))
            if sum(Quadrant==j)>0 || PenaltyValue~=inf
                %Create penalty term for incorrect quadrants.
                Penalty=zeros(size(Dist));
                Penalty(Quadrant~=j)=PenaltyValue;
                %Find location and distance to closest point with penalty.
                [~,LocQ(j)]=min(Dist+Penalty);
                DisQ(j)=Dist(LocQ(j));
                %Find the value at that point.
                ValQ{j}={fTVec{1,1}(:,:,LocQ(j)),fTVec{1,2}(:,:,LocQ(j));...
                    fTVec{2,1}(:,:,LocQ(j)),fTVec{2,2}(:,:,LocQ(j))};
            else
                %Outside of data range, calculate based on remaining points.
            end
        end
        %Calculate output using weightings.
        for k=1:4 %xx,xy,yx,yy.
            if DisQ(1)==0
                %Error catching for exact matches.
                f_pred{k}(:,:,i)=ValQ{1}{k};
                NoPts(i)=0;
            else
                for j=1:(8-(4*(G3(i)==0))) %Quadrants.
                    if ~isnan(DisQ(j))
                        f_pred{k}(:,:,i)=f_pred{k}(:,:,i)+...
                            (ValQ{j}{k}./DisQ(j))/nansum(1./DisQ);
                    end
                end
            end
        end
        %Record minimum distance and number of points used.
        MinDis(i)=min(DisQ);
        NoPts(i)=sum(1-isnan(DisQ));
    end
else
    error('Unrecognised Cofunc');
end
end