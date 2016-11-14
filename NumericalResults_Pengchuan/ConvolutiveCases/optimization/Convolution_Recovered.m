% convolution method
% using fminsearchcon to solve global optimization
% everything same, except for optimization

clear all; close all;
global A11 A12 A21 A22 N q; 
showMix_choice=1;        % showmix_choice =1, play the mixed sounds
show_s_choice=1;         % if =1, show the computed s as well as v. if =0, only show v
check_corr_choice=1;     % check the correlation between clean_sig, mix and separated sources
check_corr_length=20;    % with time lag from 0 to check_corr_length-1
rho_maxaverage_choice=0; % if choice=0, use maximum, if choice=1, use average

% checkxAxNorm1_choice;  % if =1, divide norm(x1,1) and norm(x2,1) in checkx1Anx2
% hKnown_choice;         % if >=1, use clean signal to generate the mixed sound  
%                        % choice =1, q=3; choice =2, q=25; choice =3, q=50
%                        % if <=0, load pre-mixed data
% norm_choice;           % if =1, use L^1, if =2, use L^2;
%                        % if =0, use a_1^{11}=1=a_1^{22}                            
% normcoef;
% nn;                    % break the data into nn frames
% sparse_choice;         % if >=1, make use of the sparsity by estimating spikes positions 
%                        % =1, select highest and zero neighbor entry 
%                        % =2, select 3 consecutive positions around each large spikes
% sparsity;              % sparsity/4 is the number of spikes in each a^{ij} 
% N;                     % nn*(2N+1) will be the total number of equations form An. There are 4q unknowns. 
%                        % We will set N=nn*(2N+1) later.                  
%
% if hKnown_choice==0; premix_choice=1;  end    % % if h is not known, load premix_choice
% if norm_choice>=1; objfun_choice=1;  end      % % if objfun_choice=1, use |x1Ax2|-sigma(|x|-1), if =0, use |x1Ax2|/|x1||x2|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting the parameters; go directly down to "otherwise" case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

casenum=2;

switch casenum
    case 1 
        checkxAxNorm1_choice=0; % divide norm(x1,1) and norm(x2,1) in checkx1Anx2
        hKnown_choice=0; % load pre-mixed data
        norm_choice=1; % use L^1
        normcoef=0.002;  
        sparse_choice=1; % make use of the sparsity by estimating spikes positions; select highest and zero neighbor entry 
        sparsity=20; % sparsity/4 is the number of spikes in each a^{ij}   
        nn=1; % break the data into nn frames   
        if hKnown_choice==0; premix_choice=1;  end % load pre-mixed data
        if norm_choice>=1; objfun_choice=1;  end % use |x1Ax2|-sigma(|x|-1)
        N=150;  if nn==1; N=150;  elseif nn==41; N=4;   end
    otherwise
        checkxAxNorm1_choice=0; % divide norm(x1,1) and norm(x2,1) in checkx1Anx2
        hKnown_choice=1; % use clean signal to generate the mixed sound
        norm_choice=0; % use L^1  
        normcoef=0.002;  
        sparse_choice=1; % make use of the sparsity by estimating spikes positions; select highest and zero neighbor entry   
        sparsity=20;  % sparsity/4 is the number of spikes in each a^{ij}  
        nn=1; % break the data into nn frames   
        if hKnown_choice==0; premix_choice=1;  end  % load pre-mixed data
        if norm_choice>=1; objfun_choice=1;  end % use |x1Ax2|-sigma(|x|-1)
        N=200;  if nn==1; N=200;  elseif nn==41; N=4;   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load voice/music data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hKnown_choice>=1 % use clean signal to generate the mixed sound  
    [clean_sig,mixed_x,Xex,q,nonzero_ind_ex,nonzero_Xex]=loadmixH(casenum, hKnown_choice,norm_choice); 
    A0 = [Xex(1) Xex(q+1); Xex(2*q+1) Xex(3*q+1)];
    condA = cond(A0);
    % if there are more than 3 input, load randn random variable data
else
    mixed_x=loadmix(premix_choice);      q=70;
end

if showMix_choice==1;  show_play(mixed_x,2,['mixed']);  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate 2nd order statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=q+N-1;          % M is the size of total shifting in one direction. 2*M+L<=size(mixed_x,2)/nn;
L=floor(size(mixed_x,2)/nn)-2*M-50;  % L is the size of data used to estimate 2nd order statistics
rescaleL=100;      % Aij=Aij*rescaleL so that Aij is of order 1
A11=zeros(q,q,nn*(2*N+1)); A12=A11; A21=A11; A22=A11; 
B=zeros(4,2*M+1);

% Calculating the expectation of the recovered signals
for frm=1:nn           
    initL=(frm-1)*L;
	for i=1:2 
        for j=1:2
            for m=-M:M
                B(i*2+j-2,M+1+m)=mixed_x(i,initL+(M+1:M+L))*mixed_x(j,initL+(M+1-m:M+L-m))'/(L/rescaleL); %Equation 5.103 in reduced dimensions
                % y_j shift to left by m            
            end
        end
	end
    initN=(frm-1)*(2*N+1); % 2N+1 comes from n=-N:N; N-(-N)=2N
	
    % simpler format for writing the loop (i.e. i=1,j=1,i=1,j=2..)
    for k=1:q
        for m=1:q
           for n=-N:N
               % M+1+1-N-q > = 1. That's why M=q+N-1
               A11(k,m, initN+ N+1+n)=B(1,M+1+m+n-k); % y_j shift to left by m+n-k
               A12(k,m, initN+ N+1+n)=B(2,M+1+m+n-k); 
               A21(k,m, initN+ N+1+n)=B(3,M+1+m+n-k); 
               A22(k,m, initN+ N+1+n)=B(4,M+1+m+n-k); 
           end
       end
	end
end

C11=zeros(2*q-1,2*q-1,nn*(2*N+1)); C12=C11; C21=C11; C22=C11; 
D=zeros(4,2*M+1);

for frm=1:nn
    initL=(frm-1)*L;
	for i=1:2 
        for j=1:2
            for m=-M:M
                D(i*2+j-2,M+1+m)=clean_sig(i,initL+(M+1:M+L))*clean_sig(j,initL+(M+1-m:M+L-m))'/(L/rescaleL); 
                % y_j shift to left by m            
            end
        end
	end
    initN=(frm-1)*(2*N+1);
	for k=1:2*q-1
        for m=1:2*q-1
           for n=-2*(2*q-1)^2:2*(2*q-1)^2
               % M+1+1-N-q > = 1. That's why M=q+N-1
               C11(k,m, initN+ N+1+n)=D(1,M+1+m+n-k); % y_j shift to left by m+n-k
               C12(k,m, initN+ N+1+n)=D(2,M+1+m+n-k); 
               C21(k,m, initN+ N+1+n)=D(3,M+1+m+n-k); 
               C22(k,m, initN+ N+1+n)=D(4,M+1+m+n-k); 
           end
       end
	end
end
cs = zeros(4*(2*q-1)^2+1,2*(2*q-1)^2);
for i = -2*(2*q-1)^2:2*(2*q-1)^2
    cs(i+2*(2*q-1)^2+1,:) = [reshape(C11(:,:,initN+ N+1+i),1,(2*q-1)^2), reshape(C22(:,:,initN+ N+1+i),1,(2*q-1)^2)];
end
Cs = sqrt(norm(inv(cs'*cs)));

N=nn*(2*N+1); % N>=4*q-2, and now N is the number of equations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fval0 = CorrelationConv(Xex);
ObjectiveFunction = @CorrelationConv;
X0=zeros(4*q-2,1);   % Starting point
lb = -0.5*ones(1, 4*q-2);
ub = 0.5*ones(1, 4*q-2);
% change the parameters accordingly
options = optimset('MaxFunEvals',10^6, 'MaxIter', 10^7, 'TolFun',1e-7);
[x,fval,exitflag,output] = fminsearchcon(ObjectiveFunction,X0,lb,ub,[],[],[],options);
xopt = [1;x(1:q-1);x(q:2*q-1);x(2*q:3*q-1);1;x(3*q:4*q-2)];

%%% recover the sources %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% v... is recovered v, s... is recovered s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postprocess B A P
[vX,sX]=RecoverSource2(mixed_x,xopt,q,'from optimization method',show_s_choice);
if hKnown_choice>=1
    [BB, detB, P, normP, sigmaP] = postprocess(q, xopt, Xex);
	[vXex,sXex]=RecoverSource2(mixed_x,Xex,q,'from exact method',show_s_choice);
%    [envelope_Xs_exactposition,vXs_exactposition,sXs_exactposition]=...
%        RecoverSource(mixed_x,Xs_exactposition,q,'from sparse LM with exact init position',show_s_choice);
end

return;