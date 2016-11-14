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


% generate random mixing matrix, sample size Nsample
Nsample = 100;
x0_seed = 2*rand(Nsample,2)-1;
% x0_all = pi/2*repmat([1,1],Nsample,1).*x0_seed;
% these two lines are for degenerate case.
x0_all(:,1) = pi/2*x0_seed(:,1);
x0_all(:,2) = x0_all(:,1) + 0.001*x0_seed(:,2);
vs_all = cell(Nsample,3);
Siri_all = zeros(Nsample,3);
sigmaP_all = zeros(Nsample,3);
condA_all = zeros(Nsample,1);

% for each sample, do the test
for t = 1 : Nsample
    t
    % form mixing matrix
    x0 = x0_all(t,:);
    phi0 = x0(1);
    theta0 = x0(2);

    A0 = [sin(phi0), sin(theta0); cos(phi0), cos(theta0)];
    
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
    x0qr = initialGuess2D(cs,msQ);
    
    % infomax
    % [BBinf,Sinf]=demix_soft_constrained_infoMax_siri(ms,2);
    [BBinf,Sinf]=demix_soft_constrained_infoMax_siri(msQ,2);
    vs_all{t,2} = Sinf;
    Pinf = Sinf*cs'/(cs*cs');
    Pinfabs = abs(Pinf);
    sigmaP_all(t,2) = max(min(Pinfabs(1,1)/Pinfabs(1,2),Pinfabs(2,2)/Pinfabs(2,1)),min(Pinfabs(1,2)/Pinfabs(1,1),Pinfabs(2,1)/Pinfabs(2,2)));
    % Siriinf = Siricompute2D( cs, A0, Pinf );
    Siriinf = Siricompute2D( cs, A0qr, Pinf );
    Siri_all(t,2) = Siriinf;
    
    % optimization
    % optimization setting
    [ BBopt, vsopt, ~, ~,~,~ ] = demix_min_cov_2d( msQ, x0qr );
    vs_all{t,1} = vsopt;
    Popt = vsopt*cs'/(cs*cs');
    vecP = abs(reshape(Popt,1,4));
    vecP = sort(vecP, 'descend');
    sigmaP_all(t,1) = vecP(2)/vecP(3);
    Siriopt = Siricompute2D( cs, A0qr, Popt );
    Siri_all(t,1) = Siriopt;
    
    % analytic method
    % if zeroflag = 0, we use shifts 1,2
    % if zeroflag = 1, we use shifts 0,1
    zeroflag = 1;
    [ BBana, vsana, ~ ] = demix_analytic( ms, zeroflag );
    vs_all{t,3} = vsana;
    Pana = vsana*cs'/(cs*cs');
    vecP = abs(reshape(Pana,1,4));
    vecP = sort(vecP, 'descend');
    sigmaP_all(t,3) = vecP(2)/vecP(3);
    Siriana = Siricompute2D( cs, A0, Pana );
    Siri_all(t,3) = Siriana;
end

% plot the results
figure(1)
% plot(1:Nsample,Siri_all)
plot(1:Nsample,Siri_all(:,1),'ro',1:Nsample,Siri_all(:,2),'g*',1:Nsample,Siri_all(:,3),'b^')
legend('opt','infomax','analytic')
title('SIRI')
figure(2)
% plot(1:Nsample,sigmaP_all)
plot(1:Nsample,sigmaP_all(:,1),'ro',1:Nsample,sigmaP_all(:,2),'g*',1:Nsample,sigmaP_all(:,3),'b^')
legend('opt','infomax','analytic')
title('sigmaP')