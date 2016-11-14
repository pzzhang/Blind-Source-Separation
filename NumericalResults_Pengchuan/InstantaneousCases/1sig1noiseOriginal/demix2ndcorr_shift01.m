%tested  07/06/2013 (JX)

%% problem setting
clear all;
close all;
L=36000; % set uniform signal length
[sig1,fs]=wavread('data/s3.wav'); % load wav files; returns the sample rate (Fs) in Hertz used to encode the data 
cs(1,1:L)=sig1(1:L)/norm(sig1(1:L));
% [sig2,fs]=wavread('data/s1.wav'); % load wav files; returns the sample rate (Fs) in Hertz used to encode the data 
% cs(2,1:L)=sig2(1:L)/norm(sig2(1:L));
cs(2,1:L)=randn(1,L);
cs(2,1:L) = cs(2,1:L)/norm(cs(2,1:L));
% Cstable conditions
Lcs=length(cs);
C0=zeros(2,2,2); Cs_shift = zeros(2,2);
for shif=0:1 % compute correlation with shifts
    C0(:,:,shif+1)=cs(:,1:end-shif)*(cs(:,shif+1:end)')./(L-shif);
    Cs_shift(shif+1,:) = diag(C0(:,:,shif+1))';
end
Cs = norm(inv(Cs_shift),2);

figure(1) % plot two clean signals
plot(cs(1,:))
title('first clean signal');
figure(2) % plot two clean signals
plot(cs(2,:))
title('second clean signal');

A0=[-1.2944     0.7143
    -1.3362    1.6236]; % mixing matrix
%A0=randn(2,2); % mixing matrix can be replaced by this (a random matrix)
% Asigma assumpution
condA = cond(A0);
ms=A0*cs;

figure(3) % plot two mixtures
plot(ms(1,:))
title('first mixture');
figure(4) % plot two mixtures
plot(ms(2,:))
title('second mixture');

%% optimization setting
Lms=length(ms);
C=zeros(2,2,2); d=zeros(2,1); f=zeros(2,1); A=zeros(2,2); B=zeros(2,2,2);
for shif=0:1 % compute correlation with shifts
    C(:,:,shif+1)=ms(:,1:end-shif)*(ms(:,shif+1:end)')./(Lms-shif);
end
a=C(2,1,1)*C(2,2,2)-C(2,2,1)*C(2,1,2);
b=C(2,2,1)*C(1,1,2)+C(1,2,1)*C(2,1,2)-C(2,1,1)*C(1,2,2)-C(1,1,1)*C(2,2,2);
c=C(1,1,1)*C(1,2,2)-C(1,2,1)*C(1,1,2);
if b^2 >= 4*a*c
th(1)=(-b+sqrt(b^2-4*a*c))/(2*a); % define two tan theta values 
th(2)=(-b-sqrt(b^2-4*a*c))/(2*a);
else
    th(1)=-b/(2*a); th(2)=th(1);
end
for i=1:2
    f(i)=(C(1,1,1)-th(i)*C(2,1,1))/(C(1,2,1)-th(i)*C(2,2,1)); % computing tan phi 
end

E1 = C(1,1,1) - th(1)*C(2,1,1) - f(1)*C(1,2,1) + th(1)*f(1)*C(2,2,1);
E2 = C(1,1,2) - th(1)*C(2,1,2) - f(1)*C(1,2,2) + th(1)*f(1)*C(2,2,2);


th=atan(th); f=atan(f);
for i=1:2
    B(:,:,i)=[sin(f(i)),sin(th(i)); cos(f(i)),cos(th(i))]; % a11=sin(phi),a12=sin(theta),...
end
i = 2;
BB=[cos(th(i)),-sin(th(i)); -cos(f(i)),sin(f(i))]; % BB is the inverse of B
detB = det(BB);
rs1=BB*ms; % estimate sources
figure(5) % plot two received signals
plot(rs1(1,:))
title('first recovered signal');
figure(6) % plot two received signals
plot(rs1(2,:))
title('second recovered signal');

P=BB*A0; % approximate inverse times mixing matrix
vecP = abs(reshape(P,1,4));
vecP = sort(vecP, 'descend');
sigmaP = vecP(2)/vecP(3);
% fix column permuatation to ensure diagonal dominance
if abs(P(1,1))+abs(P(2,2)) < abs(P(1,2))+abs(P(2,1)) 
    A1=[0 1; 1 0]*P;
else
    A1 = P;
end
soundsc(rs1(1,:),fs);
disp('press any key to hear the 2nd speech'), pause;
soundsc(rs1(2,:),fs);

r=norm(A0(1,1)*cs(1,:))/norm(A0(1,2)*cs(2,:))+...
  norm(A0(2,2)*cs(2,:))/norm(A0(2,1)*cs(1,:));
SIRi=10*log10(r); %input SIR (in dB)
r=norm(A1(1,1)*cs(1,:))/norm(A1(1,2)*cs(2,:))+...
  norm(A1(2,2)*cs(2,:))/norm(A1(2,1)*cs(1,:));
SIRo=10*log10(r); % output SIR (in dB)
SIRI=SIRo-SIRi % SIR improvement (in dB) following Eq.(5.20)
