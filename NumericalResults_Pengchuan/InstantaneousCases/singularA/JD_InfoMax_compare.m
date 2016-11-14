%% problem setting
clear all;
close all;
L=36000; % set uniform signal length
[sig1,fs]=wavread('data/s3.wav'); % load wav files; returns the sample rate (Fs) in Hertz used to encode the data 
cs(1,1:L)=sig1(1:L)/norm(sig1(1:L));
%[sig2,fs]=wavread('data/s1.wav'); % load wav files; returns the sample rate (Fs) in Hertz used to encode the data 
%cs(2,1:L)=sig2(1:L)/norm(sig2(1:L));
% load 'cr.mat';
cs(2,1:L)= 0.1*randn(1,L);
cs(2,1:L) = cs(2,1:L)/norm(cs(2,1:L));
% Cstable conditions
Lcs=length(cs);
C0=zeros(2,2,2); Cs_shift = zeros(2,2);
for shif=1:2 % compute correlation with shifts
    C0(:,:,shif)=cs(:,1:end-shif)*(cs(:,shif+1:end)')./(L-shif);
    Cs_shift(shif,:) = diag(C0(:,:,shif))';
end
Cs = norm(inv(Cs_shift),2);

T = 100;
sigmaP = zeros(2, T);
tic
for t = 1:T
%     t
    testJointDiag;
    sigmaP(1,t) = sigmaP_JD;
    sigmaP(2,t) = sigmaP_InfoMax;
end
toc

figure(1)
semilogy(1:T, sigmaP(1,:),'r.','MarkerSize',20)
hold on
semilogy(1:T, sigmaP(2,:),'g*')
legend('SimulDiag','InfoMax')
title('SigmaP')
