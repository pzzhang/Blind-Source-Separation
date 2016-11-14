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

% generate mixing matrix
A0 = eye(3);
ms = A0*cs;

% start from different initial guess, find the minimum
Nsample = 10;
x_all = zeros(Nsample, 6);
fval_all = zeros(Nsample,1);
vs_all = cell(Nsample,1);
Siri_all = zeros(Nsample,1);
sigmaP_all = zeros(Nsample,1);
exitflag_all = zeros(Nsample,1);


ObjectiveFunction = @(x) simple_objective(x, ms);
lb = -pi/2*[2,1,2,1,2,1];
ub = pi/2*[2,1,2,1,2,1];
x0 = [0,0,pi/2,0,0,pi/2];
fvalref = ObjectiveFunction(x0)

% for each sample, do the test
for t = 1 : Nsample
    t
    % optimization
    % optimization setting
    fval = 10;
    trial = 0;
    while (fval > fvalref) && (trial <= 20)
        trial = trial + 1;
        % initial guess
        xperturb = 1;
        X0 = x0+xperturb*pi/4*(2*rand(1,6)-1);
        for i = 1 : 3
            X0(2*i-1) = X0(2*i-1)-2*pi*floor((X0(2*i-1)+pi)/(2*pi));
            X0(2*i) = X0(2*i)-pi*floor((X0(2*i)+pi/2)/(pi));
        end
        % optimization solver
        options = optimset('MaxFunEvals',10^6, 'MaxIter', 10^7, 'TolFun',1e-13, 'TolX', 1e-13);
        [x,fval,exitflag,output] = fminsearchcon(ObjectiveFunction,X0,lb,ub,[],[],[],options);
    end
    x_all(t,:) = x;
    fval_all(t,1) = fval;
    if fval > fvalref
        exitflag = 0;
    else
        exitflag = 1;
    end
    
    % form BB and vs
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
    vs = BB*ms;
    vs_all{t,1} = vs;
    Popt = vs*cs'/(cs*cs');
    vecP = abs(reshape(Popt,1,9));
    vecP = sort(vecP, 'descend');
    sigmaP_all(t,1) = vecP(3)/vecP(4);
    Siriopt = Siricompute( cs, A0, Popt );
    Siri_all(t,1) = Siriopt;
    exitflag_all(t,1) = exitflag;
end

% plot the results
% figure(1)
% scatter(1:Nsample,Siri_all(:,1),[],exitflag_all,'filled')
% title('SIRI')
figure(1)
scatter(1:Nsample,sigmaP_all(:,1),[],exitflag_all,'filled')
title('sigmaP')