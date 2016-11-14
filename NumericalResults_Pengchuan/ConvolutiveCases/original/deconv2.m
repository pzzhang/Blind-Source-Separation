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
  
casenum=2;

switch casenum
    case 1 
        checkxAxNorm1_choice=0;
        hKnown_choice=0; 
        norm_choice=1;   normcoef=0.002;  sparse_choice=1;  sparsity=20;  nn=1;    
        if hKnown_choice==0; premix_choice=1;  end 
        if norm_choice>=1; objfun_choice=1;  end 
        N=150;  if nn==1; N=150;  elseif nn==41; N=4;   end
    otherwise
        checkxAxNorm1_choice=0;
        hKnown_choice=1; 
        norm_choice=1;   normcoef=0.002;  sparse_choice=1;  sparsity=20;  nn=1;    
        if hKnown_choice==0; premix_choice=1;  end 
        if norm_choice>=1; objfun_choice=1;  end 
        N=200;  if nn==1; N=200;  elseif nn==41; N=4;   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load voice/music data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hKnown_choice>=1
    [clean_sig,mixed_x,Xex,q,nonzero_ind_exact,nonzero_Xex]=loadmixH(casenum, hKnown_choice,norm_choice);
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
L=floor(size(mixed_x,2)/nn)-2*M-50;  % L is the size of data used to estiamte 2nd order statistics
rescaleL=100;      % Aij=Aij*rescaleL so that Aij is of order 1
A11=zeros(q,q,nn*(2*N+1)); A12=A11; A21=A11; A22=A11; 
B=zeros(4,2*M+1);

for frm=1:nn
    initL=(frm-1)*L;
	for i=1:2 
        for j=1:2
            for m=-M:M
                B(i*2+j-2,M+1+m)=mixed_x(i,initL+(M+1:M+L))*mixed_x(j,initL+(M+1-m:M+L-m))'/(L/rescaleL); 
                % y_j shift to left by m            
            end
        end
	end
    initN=(frm-1)*(2*N+1);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call LM with full x. Solution is called X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fval0 = CorrelationConv(Xex);
opts=[1e-3, 1e-7, 1e-12, 1000, 1e-15];

if norm_choice==0   % means a_1^{11}=a_1^{22}=1
   X0=zeros(4*q-2,1);
   [X, inform] = marquardt('fun_marquardt', X0, opts);     
   X=[1;X(1:q-1);X(q:2*q-1);X(2*q:3*q-1);1;X(3*q:4*q-2)];
   X0=[1;X0(1:q-1);X0(q:2*q-1);X0(2*q:3*q-1);1;X0(3*q:4*q-2)];
else   
    X0=[1;zeros(q-1,1); zeros(2*q,1); 1;zeros(q-1,1)];
    if objfun_choice==1
        [X, inform] = marquardt('fun_L12marquardt', X0, opts, normcoef, norm_choice);
	elseif objfun_choice==0
        [X, inform] = marquardt('fun_Quotientmarquardt', X0, opts, norm_choice);
	end
end
inform_str=['info(6)=1 : Stopped by small gradient;  =2 : Stopped by small x-step;  = 3 : No. of iteration exceeded'];
disp(inform);  disp(inform_str);  disp(inform(6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call LM again with estimated position of nonzero spikes. Solution is called Xs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sparse_choice>=1      
	[nonzero_ind, xs0]=select_spikes(X,sparsity,sparse_choice);  % little x means truncated X and X is full. 
    % use (1,0,0,...) as initial data
    nonzero_ind=union(nonzero_ind,[1,3*q+1]);
    xtemp=[1,zeros(1,3*q-1),1,zeros(1,q-1)];  xs0=xtemp(nonzero_ind); 
    Xs0=zeros(4*q,1);  Xs0(nonzero_ind)=xs0;  % zero padding for drawing the figure later on
    
	% call marquardt again with estimated position of nonzero spikes
	[xs, inform] = marquardt('fun_sparseL12marquardt', xs0, opts, normcoef,norm_choice,nonzero_ind);
	disp(inform);  disp(inform_str);  disp(inform(6));
	Xs=zeros(4*q,1);   Xs(nonzero_ind)=xs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call LM again with exact position of nonzero spikes. Solution is called Xs_exactposition. s for sparse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hKnown_choice>=1
    xs0_exactposition = zeros(size(nonzero_Xex));
    [bb,ind]=sort(abs(nonzero_Xex));  xs0_exactposition(ind(end-1:end))=1;  % xs0_exactposition = something like [1,0,0,0,...,0] 

    % call marquardt again with exact position of nonzero spikes
    [xs_exactposition, inform] = marquardt('fun_sparseL12marquardt', xs0_exactposition, ...
                               opts, normcoef,norm_choice,nonzero_ind_exact);
	disp(inform);  disp(inform_str);  disp(inform(6));
    % little xs_... will be padded by zero to big Xs_...
	Xs_exactposition=zeros(4*q,1);   Xs_exactposition(nonzero_ind_exact)=xs_exactposition;  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% check the independence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hKnown_choice>=1
	[fval,x1Ax2]=checkx1Anx2(checkxAxNorm1_choice, Xex); 
	[fobj]=objfunval(Xex,normcoef,norm_choice,objfun_choice);
	disp('Exact Xex: (sum_n  |x1Anx2|^p)^1/p with p=[1/2,1,2] and obj-f(x)'); vpa([x1Ax2,fobj],3)
end

[fval,x1Ax2]=checkx1Anx2(checkxAxNorm1_choice, X); 
[fobj]=objfunval(X,normcoef,norm_choice,objfun_choice);
disp('full LM: (sum_n  |x1Anx2|^p)^1/p with p=[1/2,1,2] and obj-f(x)'); vpa([x1Ax2,fobj],3)

if hKnown_choice>=1
    [fval,x1Ax2]=checkx1Anx2(checkxAxNorm1_choice, Xs_exactposition); 
    [fobj]=objfunval(Xs_exactposition,normcoef,norm_choice,objfun_choice);
    disp('sparse LM with exact position: (sum_n  |x1Anx2|^p)^1/p with p=[1/2,1,2] and obj-f(x)'); vpa([x1Ax2,fobj],3)
end

if sparse_choice>=1; 
    [fval,x1Ax2]=checkx1Anx2(checkxAxNorm1_choice, Xs); 
    [fobj]=objfunval(Xs,normcoef,norm_choice,objfun_choice);
    disp('sparse LM with estimated position: (sum_n  |x1Anx2|^p)^1/p with p=[1/2,1,2] and obj-f(x)'); vpa([x1Ax2,fobj],3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% show weights %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
myfigure; 
draw_weights(X,'b',['FLM ']);
if hKnown_choice>=1; 
   draw_weights(Xex,'k.',['FLM ';'Ex  ']);
   draw_weights(Xs_exactposition,'g-.',['FLM ';'Ex  ';'SLMx']);
end 
if sparse_choice>=1;
   if hKnown_choice>=1; draw_weights(Xs,'r:',['FLM ';'Ex  ';'SLMx';'SLMs']);
   else   draw_weights(Xs,'r:',['FLM ';'SLMs']);
   end
end
    
%%% recover the sources %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% v... is recovered v, s... is recovered s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vX,sX]=RecoverSource2(mixed_x,X,q,'from full LM',show_s_choice);
if hKnown_choice>=1
    [BB, detB, P, normP, sigmaP] = postprocess(q, X, Xex);
	[vXex,sXex]=RecoverSource2(mixed_x,Xex,q,'from exact',show_s_choice);
%    [envelope_Xs_exactposition,vXs_exactposition,sXs_exactposition]=...
%        RecoverSource(mixed_x,Xs_exactposition,q,'from sparse LM with exact init position',show_s_choice);
end

return;

