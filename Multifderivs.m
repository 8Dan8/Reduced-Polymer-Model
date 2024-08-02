function [df,lambda,dfflow]=Multifderivs(f,Time,N,Z,FlowType,StepET_R,cnu,pqmatrix)
%MULTIFDERIVS finds the derivative of f at the specified times and 
%parameters for all locations specified within pqmatrix.
%
%Inputs:
%f - Cell array (in a,b) of matrices (in p,q,t).
%Time - Specifies the times at which f was evaluated.
%N - N from input parameters.
%Z - Z from input parameters.
%FlowType - Type of flow as a character string. 'xxuniext','yyuniext' and 
%           'xyshear' are supported.
%StepET_R - Cell array {number of steps(int), ET_R(vec), Start times(vec)}.
%cnu - scalar value for cnu.
%pqmatrix - N by N logical array to denote which p,q pairs to calculate. If
%           empty, all elements shall be calculated.
%
%Outputs:
%df - Cell array of matrices that contain df/dt. (In the same format as f)
%lambda - Vector of lambda over time.

%Default pq matrix to all points.
if exist('pqmatrix','var')==0
    pqmatrix=ones(N);
    warning('pqmatrix not defined, defaulting to all elements. (Slow)')
end
%Default cnu to 0.0.
if exist('cnu','var')==0 || isempty(cnu)
    cnu=0;
    warning('cnu not defined, defaulted to cnu=0.0.')
end
%Extract info about F.
Tlen=size(f{1,1},3);
A=int8(size(f,1));
%Create vectors for extdot. (Note StepET_R defines {stepno,ET_R,starttime})
stepindex=ones(1,Tlen);
for i=2:StepET_R{1}
    stepindex=stepindex+(Time>=StepET_R{3}(i));
end
extdot=(StepET_R{2}(stepindex)/(Z^2));
%Fixed paramters.
taue=1;
Z=int8(Z);
N=int8(N);
lambdam=10000000;

%Preallocate.
finital=zeros(N+1,N+1,A,A);
lambda=zeros(1,Tlen);
%Create inital f.
for p=0:N
    for q=0:N
        for a=1:A
            finital(p+1,q+1,a,a)=((abs(p-q)*Z/N)<0.5)/3;
        end
    end
end
%Preallocate df and dfflow.
df=cell(A);
for k=1:(A^2)
    df{k}=zeros(N,N,Tlen);
end
dfflow=df;
%Wrap F and turn into a cell vector for t.
for t=1:Tlen
    fpad=zeros(501,501,3,3); 
    fpad(1:(N+1),1:(N+1),1:A,1:A)=finital;
    for i=1:A
        for j=1:A
            %Note that F indexed from 0 for p and q in fortran.
            fpad(2:N+1,2:N+1,i,j)=f{i,j}(:,:,t);
        end
    end
    %Find the value of lambda for the current time.
    lambda(t)=flambda(fpad,cnu,taue,extdot(t),Z,N,lambdam);
    %Loop over these parameters.
    for p=1:N-1
        for q=1:N-1
            %Evaluate only at specific points.
            if pqmatrix(p,q)
                %Apply symmetries.
                if p>q && pqmatrix(q,p) 
                    %p>q so df{:,:}(q,:,:) are already calculated.
                    for a=1:A
                        for b=1:A
                            df{a,b}(p,q,t)=df{a,b}(q,p,t);
                        end
                    end
                elseif p+q>N && pqmatrix(N-q,N-p)
                    %p>N-q so df{:,:}(N-q,:,:) are already calculated.
                    for a=1:A
                        for b=1:A
                            df{a,b}(p,q,t)=df{a,b}(N-q,N-p,t);
                        end
                    end
                else
                    %Symmetric element not calculated, calculate element
                    %directly.
                    if regexpi(FlowType,'xx\s*u\w*\s*e\w*')
                        %xx Extension [e,0,0;0,-e/2,0;0,0,-e/2]
                        dfflow{1,1}(p,q,t)=2*extdot(t).*f{1,1}(p,q,t);
                        dfflow{2,2}(p,q,t)=-extdot(t).*f{2,2}(p,q,t);
                    elseif regexpi(FlowType,'yy\s*u\w*\s*e\w*')
                        %yy Extension [-e/2,0,0;0,e,0;0,0,-e/2]
                        dfflow{1,1}(p,q,t)=-extdot(t).*f{1,1}(p,q,t);
                        dfflow{2,2}(p,q,t)=2*extdot(t).*f{2,2}(p,q,t);
                        warning('Trace assumes yy=zz in df=fderivs.')
                    elseif regexpi(FlowType,'xy\s*sh\w*')
                        %xy Shear [0,g,0;0,0,0;0,0,0]
                        dfflow{1,1}(p,q,t)=2*extdot(t).*f{1,2}(p,q,t);
                        dfflow{1,2}(p,q,t)=extdot(t).*f{2,2}(p,q,t);
                        dfflow{2,1}(p,q,t)=dfflow{1,2}(p,q,t);
                    else
                        error('Unrecognised flow type.')
                    end
                    for a=1:A
                        %Use mex-function with generalised.
                        for b=a:A
                            df{a,b}(p,q,t)=fderivs(fpad,p,q,a,b,lambda(t),...
                                cnu,taue,0,Z,N,lambdam)+dfflow{a,b}(p,q,t);
                        end
                        %Utilise symmetry to save time.
                        for b=1:(a-1)
                            df{a,b}(p,q,t)=df{b,a}(p,q,t);
                        end
                    end
                end
            end
        end
    end
end
end

%Old mex function:
% for b=a:A
%     df{a,b}(p,q,t)=fderivs(fpad,p,q,a,b,lambda(t),...
%         cnu,taue,extdot(t),Z,N,lambdam);
% end

%Test example
%df=Multifderivs(f,Time,N,Z,'xxUE',2,0,eye(N));