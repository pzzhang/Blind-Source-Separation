function [V,est]=demix_soft_constrained_infoMax_siri(X,m)

%showMixChoice=1;
%showYChoice=0;

NS=m; % number of signals 

Ls=38000;
%[sig,fs]=wavread('s8.wav');
%cs(1,1:Ls)=sig(1:Ls);
%[sig,fs]=wavread('s1.wav');
%cs(2,1:Ls)=sig(1:Ls);
%[sig,fs]=wavread('s2.wav');
%cs(3,1:Ls)=sig(1:Ls);
%[sig,fs]=wavread('s4.wav');
%cs(4,1:Ls)=sig(1:Ls);
%     H=[-1.2944     0.7143
%        -1.3362    1.6236];

%M=[2,-1.5,-2,-1.4;
%    -1,-2,-1.2,-1;
%    1.3,1.4,-0.8,1;
%    2,-0.8,1.2,0.5];
%

%M=randn(2,2);

%mixed_x=M*cs;

mixed_x=X;

s=zeros(size(mixed_x));
L=10;
N=length(mixed_x);

%cputime1=cputime;

W=zeros(NS,NS,floor(N/L));
sig1=zeros(floor(N/L),1);
sig2=zeros(floor(N/L),1);
F1=zeros(floor(N/L)-1,1);
F2=zeros(floor(N/L)-1,1);
W(:,:,1)=eye(NS);
sig1(1)=0.01;
sig2(1)=0.01;
A=zeros(NS,NS);

y=W(:,:,1)*mixed_x(:,1:L); % y is the estiamted source
d=phi(y); % phi is -rho'/rho, should have phi*y' =Id
for j=1:NS
    for k=1:NS
        A(j,k)=d(j,:)*y(k,:)'/L;
    end
end
H=A*W(:,:,1); % This is H(1)
%nW=zeros(floor(N/L),1);
%nW(1)=norm(W(:,:,1),1);
for i=1:floor(N/L)-1;
    s(:,(i-1)*L+1:i*L)=y;

    nu=0.01;%1/sqrt(i+1); % this is nu(i)
    nu_ip1=nu;
    W(:,:,i+1)=(1+nu*sig1(i))*W(:,:,i)-nu*sig2(i)*H;

    y=W(:,:,i+1)*mixed_x(:,i*L+1:(i+1)*L);
    d=phi(y); % phi is -r'/r 
    for j=1:NS
        for k=1:NS
            A(j,k)=d(j,:)*y(k,:)'/L;
        end
    end
    H=A*W(:,:,i+1); % This is H(i+1)
    
    lam=0.01;%0.009;
    a=1;
    F1(i)=(1+nu*sig1(i))*norm(W(:,:,i+1),1)+lam*(1+nu*sig1(i))-a;%+lam*sig1(i)-a;
    F2(i)=sig2(i)*norm(H,1)+lam*sig2(i)-a;
    sig1(i+1)=((1+nu*sig1(i))*exp(-F1(i))-1)/nu_ip1;
    sig2(i+1)=sig2(i)*exp(-F2(i));
%    nW(i+1)=norm(W(:,:,i+1),1);
end

%disp(['total cpu time for the iteration is ',num2str(cputime-cputime1)]);

V=W(:,:,end);
for i=1:NS
 %   [a,id1]=max(abs(Wex(i,:)));
    [b,id2]=max(abs(V(i,:)));
%    Wex(i,:)=Wex(i,:)./(sign(Wex(i,id1))*a);
    V(i,:)=V(i,:)./(sign(V(i,id2))*b);    
end

est=V*mixed_x;


%disp('press any key to hear estimated source');
%    for i=1:NS
%        pause(); soundsc(est(i,:),fs);
%    end
return;


function d=phi(p)
% d=p+tanh(p);
% d=p-tanh(p);
d=p./(abs(p)+1e-16);

