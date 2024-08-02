function File = GLaMMD(Data,Deform,Stepno,ET_R,N,Z,Tmax,CCR)
%GLAMMD Returns a GLaMM output file from database that fits the specified 
%input parameters.
%
%Inputs:
%Data - Choose to return the 'Trans' or 'Fpq' file.
%Deform - Chooses the type of deformation 'e' or 'extension' chooses
%         uniaxial extension flow, 's' or 'shear' chooses a shear flow.
%Stepno - Chooses 'single' or 'step' type, 1 and 3 can also be used.
%ET_R - Epsilon Tau_R value.
%N - N value.
%Z - Z value.
%Tmax - Largest T required.

%Set defaults.
if exist('Tmax','var')==0
    Tmax=1000;
end
if exist('CCR','var')==0
    CCR=0;
end
Flag=0;
%Compare Data.
if strcmpi(Data,'f') || strcmpi(Data,'Fpq')
    Data='Fpq';
elseif strcmpi(Data,'t') || strcmpi(Data,'Trans')
    Data='Trans';
else
    Flag=1;
end
%Compare deform.
if strcmpi(Deform,'e') || strcmpi(Deform,'Extension')
    Deform='Extension';
    Var='ET_R';
elseif strcmpi(Deform,'s') || strcmpi(Deform,'Shear')
    Deform='Shear';
    Var='GT_R';
else
    Flag=1;
end
%Compare Type.
if (numel(Stepno)==1 && Stepno==3) || strcmpi(Stepno,'Step')
    Stepno='step';
elseif (numel(Stepno)==1 && Stepno==2) || strcmpi(Stepno,'Split')
    Stepno='split';
elseif (numel(Stepno)==1 && Stepno==1) || strcmpi(Stepno,'Single')
    Stepno='single';
else
    Flag=1;
end
%Calculate CCR
if CCR==0
    Beta='B';
elseif CCR==0.1
    Beta='S';
else
    Flag=1;
end
%Calculate N/Z
if N==Z
    NZ=[];
elseif mod(N/Z,1)==0
    NZ=num2str(N/Z);
else
    Flag=1;
end
%Check length
if Tmax(1)==1000
    Length=[];
elseif Tmax(1)==10000
    Length='L';
elseif Tmax(1)==100000
    Length='LL';
else
    Flag=1;
end
%Test Section
if numel(Tmax)>1 && strcmp(Tmax(2),'T')
   Length=[Length,'T'];
end
%Check Flag, produce output.
if Flag
    error('Error: unrecognised inputs.')
else
    File=[Deform,'\',Data,num2str(Z),Stepno,Var,num2str(ET_R),Beta,NZ,Length,'.dat'];
end