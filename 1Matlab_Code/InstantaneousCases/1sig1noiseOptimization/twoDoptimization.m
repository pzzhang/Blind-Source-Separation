%tested  07/06/2013 (JX)

%% problem setting
clear all;
close all;
L=36000; % set uniform signal length
[sig1,fs]=wavread('data/s3.wav'); % load wav files; returns the sample rate (Fs) in Hertz used to encode the data 
cs(1,1:L)=sig1(1:L)/norm(sig1(1:L));
% [sig2,fs]=wavread('data/s1.wav'); % load wav files; returns the sample rate (Fs) in Hertz used to encode the data 
% cs(2,1:L)=sig2(1:L)/norm(sig2(1:L));
load 'cr.mat';
cs(2,1:L)=cr;
cs(2,1:L) = cs(2,1:L)/norm(cs(2,1:L));
% Cstable conditions
Lcs=length(cs);
C0=zeros(2,2,2); Cs_shift = zeros(2,2);
for shif=1:2 % compute correlation with shifts
    C0(:,:,shif)=cs(:,1:end-shif)*(cs(:,shif+1:end)')./(L-shif);
    Cs_shift(shif,:) = diag(C0(:,:,shif))';
end
Cs = norm(inv(Cs_shift),2);

% figure(1) % plot two clean signals
% plot(cs(1,:))
% title('first clean signal');
% figure(2) % plot two clean signals
% plot(cs(2,:))
% title('second clean signal');

phi0 = pi/2*rand();
theta0 = pi/2*rand();
A0=[sin(phi0)     sin(theta0)
    cos(phi0)    cos(theta0)]; % mixing matrix
%A0=randn(2,2); % mixing matrix can be replaced by this (a random matrix)
% Asigma assumpution
condA = cond(A0);
ms=A0*cs;

% figure(3) % plot two mixtures
% plot(ms(1,:))
% title('first mixture');
% figure(4) % plot two mixtures
% plot(ms(2,:))
% title('second mixture');
fval0 = simple_objective2(A0, ms);

%% optimization setting
ObjectiveFunction = @(x) simple_objective(x, ms);
% [x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,X0)
lb = [0 0];
ub = pi/2*[1 1];

% randomly generate Ns points and get the minumum
Ns = 1000;
xsample = pi/2*rand(Ns,2);
fvalsample = zeros(Ns,1);
for i = 1:Ns
    fvalsample(i) = ObjectiveFunction(xsample(i,:));
end
[fvalsmin, index] = min(fvalsample);
% initial guess
X0 = xsample(index(1),:);

% change the parameters accordingly
options = optimset('MaxFunEvals',10^6, 'MaxIter', 10^7, 'TolFun',1e-13);
[x,fval,exitflag,output] = fminsearchcon(ObjectiveFunction,X0,lb,ub,[],[],[],options);

theta = x(1);
phi = x(2);

b11 = cos(theta);
b12 = -sin(theta);

b21 = -cos(phi);
b22 = sin(phi);

BB = [b11 b12; b21 b22];
detB = det(BB);
fval1 = simple_objective2(BB, ms);
rs1 = BB*ms;

% figure(5) % plot two received signals
% plot(rs1(1,:))
% title('first recovered signal');
% figure(6) % plot two received signals
% plot(rs1(2,:))
% title('second recovered signal');

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


% code used to call infomax
[BB,S]=demix_soft_constrained_infoMax_siri(ms,2);
A1=BB*A0;
P = BB*A0;
rs1 = BB*ms;

% soundsc(rs1(1,:),fs);
% disp('press any key to hear the 2nd speech'), pause;
% soundsc(rs1(2,:),fs);

r=min(norm(A0(1,1)*cs(1,:))/norm(A0(1,2)*cs(2,:)), norm(A0(2,2)*cs(2,:))/norm(A0(2,1)*cs(1,:)));
SIRi=10*log10(r); %input SIR (in dB)
r=min(norm(A1(1,1)*cs(1,:))/norm(A1(1,2)*cs(2,:)), norm(A1(2,2)*cs(2,:))/norm(A1(2,1)*cs(1,:)));
SIRo=10*log10(r); % output SIR (in dB)
SIRI=SIRo-SIRi % SIR improvement (in dB) following Eq.(5.20)