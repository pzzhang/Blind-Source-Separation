% three source separation using global optimization
% run this code with "simple_objective" and "fminsearchcon"
% may need to change the parameters in order for successful results
clear all 
close all
%% problem setting
L=36000; % set uniform signal length
[sig1,fs]=wavread('data/s3.wav');
cs(1,1:L)=sig1(1:L)/norm(sig1(1:L));
[sig2,fs]=wavread('data/s1.wav');
cs(2,1:L)=sig2(1:L)/norm(sig2(1:L));
[sig3,fs]=wavread('data/s2.wav');
cs(3,1:L)=sig3(1:L)/norm(sig3(1:L));
% cs(3,1:L)=randn(1,L);
% cs(3,1:L)=cs(3,1:L)/norm(cs(3,1:L));

% shift will be used
shiftarray = [0, 1, 2];
N = length(shiftarray);
% Cstable conditions
Lcs=length(cs);
C0=zeros(3,3,N); Cs_shift = zeros(N,N);
for i = 1:N % compute correlation with shifts
    shif = shiftarray(i);
    C0(:,:,i)=cs(:,1:end-shif)*(cs(:,shif+1:end)')./(Lcs-shif);
    Cs_shift(i,:) = diag(C0(:,:,i))';
end
Cs = norm(inv(Cs_shift),2);


figure(1) % plot two clean signals
plot(cs(1,:))
title('first clean signal');
figure(2) % plot two clean signals
plot(cs(2,:))
title('second clean signal');
figure(3) % plot two clean signals
plot(cs(3,:))
title('third clean signal');

theta1 = pi/6;
phi1 = pi/5;
theta2 = pi/4;
phi2 = pi/3;
theta3 = pi/2;
phi3 = pi/12;
x0 = [1/6 1/5 1/4 1/3 1/2 1/12]*pi;

b11 = cos(theta1)*cos(phi1);
b12 = sin(theta1)*cos(phi1);
b13 = sin(phi1);

b21 = cos(theta2)*cos(phi2);
b22 = sin(theta2)*cos(phi2);
b23 = sin(phi2);

b31 = cos(theta3)*cos(phi3);
b32 = sin(theta3)*cos(phi3);
b33 = sin(phi3);

x0 = [theta1 phi1 theta2 phi2 theta3 phi3];
B0 = [b11 b12 b13; b21 b22 b23; b31 b32 b33];
A0 = inv(B0);
%A0=randn(2,2); % mixing matrix can be replaced by this (a random matrix)
% Asigma assumpution
condA = cond(A0);
ms=A0*cs;
fval0 = simple_objective(x0, ms);

figure(4) % plot two mixtures
plot(ms(1,:))
title('first mixture');
figure(5) % plot two mixtures
plot(ms(2,:))
title('second mixture');
figure(6) % plot two mixtures
plot(ms(3,:))
title('third mixture');

%% optimization setting
ObjectiveFunction = @(x) simple_objective(x, ms);
% [x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,X0)
lb = [0 0 0 0 0 0];
ub = pi/2 * ones(1, 6);

% randomly generate Ns points and get the minumum
Ns = 1000;
xsample = pi/2*rand(Ns,6);
fvalsample = zeros(Ns,1);
for i = 1:Ns
    fvalsample(i) = ObjectiveFunction(xsample(i,:));
end
[fvalsmin, index] = min(fvalsample);
% initial guess
X0 = xsample(index(1),:);

% change the parameters accordingly
options = optimset('MaxFunEvals',10^6, 'MaxIter', 10^7, 'TolFun',1e-12);
[x,fval,exitflag,output] = fminsearchcon(ObjectiveFunction,X0,lb,ub,[],[],[],options);

theta1 = x(1);
phi1 = x(2);
theta2 = x(3);
phi2 = x(4);
theta3 = x(5);
phi3 = x(6);

b11 = cos(theta1)*cos(phi1);
b12 = sin(theta1)*cos(phi1);
b13 = sin(phi1);

b21 = cos(theta2)*cos(phi2);
b22 = sin(theta2)*cos(phi2);
b23 = sin(phi2);

b31 = cos(theta3)*cos(phi3);
b32 = sin(theta3)*cos(phi3);
b33 = sin(phi3);

BB = [b11 b12 b13; b21 b22 b23; b31 b32 b33];
detB = det(BB);
fval1 = simple_objective2(BB, ms);
P = BB*A0;
vecP = abs(reshape(P,1,9));
vecP = sort(vecP, 'descend');
sigmaP = vecP(3)/vecP(4);
rs1 = BB*ms;

figure(7) % plot two received signals
plot(rs1(1,:))
title('first recovered signal');
figure(8) % plot two received signals
plot(rs1(2,:))
title('second recovered signal');
figure(9) % plot two received signals
plot(rs1(3,:))
title('third recovered signal');

% v1 = b11*ms(1) + b12*ms(2) + b13*ms(3);
% v2 = b21*ms(1) + b22*ms(2) + b23*ms(3);
% v3 = b31*ms(1) + b32*ms(2) + b33*ms(3);

soundsc(rs1(1,:),fs);
disp('press any key to hear the 2nd speech'), pause;
soundsc(rs1(2,:),fs);
disp('press any key to hear the 3rd speech'), pause;
soundsc(rs1(3,:),fs);

% r=norm(A0(1,1)*cs(1,:))/(norm(A0(1,2)*cs(2,:))+norm(A0(1,3)*cs(3,:))/2)+...
%   norm(A0(2,2)*cs(2,:))/(norm(A0(2,1)*cs(1,:))+norm(A0(2,3)*cs(3,:))/2)+...
%   norm(A0(3,3)*cs(3,:))/(norm(A0(3,1)*cs(1,:))+norm(A0(3,2)*cs(2,:))/2);
% SIRi=10*log10(r); %input SIR (in dB)
% r=norm(P(1,1)*cs(1,:))/(norm(P(1,2)*cs(2,:))+norm(P(1,3)*cs(3,:))/2)+...
%   norm(P(2,2)*cs(2,:))/(norm(P(2,1)*cs(1,:))+norm(P(2,3)*cs(3,:))/2)+...
%   norm(P(3,3)*cs(3,:))/(norm(P(3,1)*cs(1,:))+norm(P(3,2)*cs(2,:))/2);
% SIRo=10*log10(r); % output SIR (in dB)
% SIRI=SIRo-SIRi % SIR improvement (in dB) following Eq.(5.20)