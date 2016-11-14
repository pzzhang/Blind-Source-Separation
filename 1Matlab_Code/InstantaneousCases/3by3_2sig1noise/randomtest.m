% three source separation using global optimization
% run this code with "simple_objective" and "fminsearchcon"
% may need to change the parameters in order for successful results
clear all 
close all
%% problem setting
L=36000; % set uniform signal length
[sig1,fs]=wavread('data/bach.wav');
cs(1,1:L)=sig1(1:L)/norm(sig1(1:L));
[sig2,fs]=wavread('data/s1.wav');
cs(2,1:L)=sig2(1:L)/norm(sig2(1:L));
% [sig3,fs]=wavread('data/s2.wav');
% cs(3,1:L)=sig3(1:L)/norm(sig3(1:L));
%cs(3,1:L)=randn(1,L);
%cs(3,1:L)=cs(3,1:L)/norm(cs(3,1:L));
load 'cr.mat';
cs(3,1:L)=cr;
cs(3,1:L) = cs(3,1:L)/norm(cs(3,1:L));

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


% figure(1) % plot two clean signals
% plot(cs(1,:))
% title('first clean signal');
% figure(2) % plot two clean signals
% plot(cs(2,:))
% title('second clean signal');
% figure(3) % plot two clean signals
% plot(cs(3,:))
% title('third clean signal');


% generate random mixing matrix, sample size Nsample
Nsample = 10;
x0_seed = 2*rand(Nsample,6)-1;
% x0_all = pi/2*repmat([2,1,2,1,2,1],Nsample,1).*x0_seed;
% this line is for degenerate case.
% x0_all = pi/2*repmat([2,1,0.002,0.001,0.002,0.001],Nsample,1).*x0_seed;
% these two lines are for degenerate case 2.
x0_all = pi/2*repmat([2,1,2,1,2,1],Nsample,1).*x0_seed;
x0_all(:,5:6) = x0_all(:,3:4) + pi/2*repmat([0.002,0.001],Nsample,1).*x0_seed(:,5:6);
vs_all = cell(Nsample,2);
Siri_all = zeros(Nsample,2);
sigmaP_all = zeros(Nsample,2);
condA_all = zeros(Nsample,1);
exitflag_all = zeros(Nsample,1);

% for each sample, do the test
for t = 1 : Nsample
    t
    % form mixing matrix
    x0 = x0_all(t,:);
    theta1 = x0(1);
    phi1 = x0(2);
    theta2 = x0(3);
    phi2 = x0(4);
    theta3 = x0(5);
    phi3 = x0(6);

    b11 = cos(theta1)*cos(phi1);
    b12 = sin(theta1)*cos(phi1);
    b13 = sin(phi1);

    b21 = cos(theta2)*cos(phi2);
    b22 = sin(theta2)*cos(phi2);
    b23 = sin(phi2);

    b31 = cos(theta3)*cos(phi3);
    b32 = sin(theta3)*cos(phi3);
    b33 = sin(phi3);

    B0 = [b11 b12 b13; b21 b22 b23; b31 b32 b33];
    A0 = inv(B0);
    % Asigma assumpution
    condA = cond(A0);
    condA_all(t) = condA;
    %form mixing signal
    ms=A0*cs;
    fval0 = simple_objective(x0, ms);
    
    % do QR factorization
    [msQ, msR] = qr(ms',0);
    msQ = msQ';
    A0qr = (msQ*cs')/(cs*cs');
    x0qr = initialGuess( msQ,cs );
    
    % infomax
    disp('Infomax')
    [BBinf,Sinf]=demix_soft_constrained_infoMax_siri(ms,3);
    % [BBinf,Sinf]=demix_soft_constrained_infoMax_siri(msQ,3);
    vs_all{t,2} = Sinf;
    Pinf = Sinf*cs'/(cs*cs');
    vecP = abs(reshape(Pinf,1,9));
    vecP = sort(vecP, 'descend');
    sigmaP_all(t,2) = vecP(3)/vecP(4);
    Siriinf = Siricompute( cs, A0, Pinf );
    % Siriinf = Siricompute( cs, A0qr, Pinf );
    Siri_all(t,2) = Siriinf;
    
    % optimization
    % optimization setting
    disp('Optimization')
    [ BBopt, vsopt, ~, ~,exitflag,~ ] = demix_min_cov_3d( msQ, x0qr );
    vs_all{t,1} = vsopt;
    Popt = vsopt*cs'/(cs*cs');
    vecP = abs(reshape(Popt,1,9));
    vecP = sort(vecP, 'descend');
    sigmaP_all(t,1) = vecP(3)/vecP(4);
    Siriopt = Siricompute( cs, A0qr, Popt );
    Siri_all(t,1) = Siriopt;
    exitflag_all(t,1) = exitflag;
end

% plot the results
figure(1)
scatter(1:Nsample,Siri_all(:,1),[],exitflag_all,'filled')
hold on
plot(1:Nsample,Siri_all(:,2),'r*')
legend('opt','infomax')
title('SIRI')
hold off
figure(2)
scatter(1:Nsample,sigmaP_all(:,1),[],exitflag_all,'filled')
hold on
plot(1:Nsample,sigmaP_all(:,2),'r*')
legend('opt','infomax')
title('sigmaP')
hold off