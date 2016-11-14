% three source separation using global optimization
% run this code with "simple_objective" and "fminsearchcon"
% may need to change the parameters in order for successful results
clear all 
close all
global Cov_v Nshift d;
%% problem setting
d = 3;
L=36000; % set uniform signal length
[sig1,fs]=wavread('data/s3.wav');
cs(1,1:L)=sig1(1:L)/norm(sig1(1:L));
[sig2,fs]=wavread('data/s1.wav');
cs(2,1:L)=sig2(1:L)/norm(sig2(1:L));
% [sig3,fs]=wavread('data/s2.wav');
% cs(3,1:L)=sig3(1:L)/norm(sig3(1:L));
cs(3,1:L)=randn(1,L);
cs(3,1:L)=cs(3,1:L)/norm(cs(3,1:L));

% shift will be used
shiftarray = [0, 1, 2];
Nshift = length(shiftarray);
% Cstable conditions
Lcs=length(cs);
Cov_s=zeros(d,d,Nshift); Cs_shift = zeros(Nshift,Nshift);
for i = 1:Nshift % compute correlation with shifts
    shif = shiftarray(i);
    Cov_s(:,:,i)=cs(:,1:end-shif)*(cs(:,shif+1:end)')./(Lcs-shif);
    Cs_shift(i,:) = diag(Cov_s(:,:,i))';
end
Cs = norm(inv(Cs_shift),2);
% reference value
fref = 0;
for i = 1 : Nshift
    Cov_temp = Cov_s(:,:,i) - diag(diag(Cov_s(:,:,i)));
    fref = fref + 0.5*sum(sum(Cov_temp.^2));
end

% mix the source signals
A0=randn(d,d); % mixing matrix can be replaced by this (a random matrix)
% Asigma assumpution
condA = cond(A0);
ms=A0*cs;
save('mixedsignal.mat','ms','fs');
% QR factorization
[Qms, Rms] = qr(ms',0);
ms = Qms'; 
A0 = (ms*cs')/(cs*cs');
 
% compute the correlation for the mixtures
Cov_v=zeros(d,d,Nshift);
for i = 1:Nshift % compute correlation with shifts
    shif = shiftarray(i);
    Cov_v(:,:,i)=ms(:,1:end-shif)*(ms(:,shif+1:end)')./(Lcs-shif);
end

%% optimization setting
opts=[1e-3, 1e-14, 1e-14 , 1000, 1e-15];
X0=[zeros(d-1,1) pi/2*triu(ones(d-1,d-1))];
x0 = reshape(X0, d*(d-1),1);
f0 = 0.5*sum(fun_GeorgeMarquardt2(x0).^2);
[x, inform] = marquardt('fun_GeorgeMarquardt2', x0, opts);

% form the recovering matrix B
X = reshape(x, d-1, d);
BB = zeros(d, d);
for i = 1 : d
    BB(:,i) = Ang2Cart(X(:,i), d);
end
BB = BB';
detB = det(BB);
P = BB*A0;
% vecP = abs(reshape(P,1,9));
% vecP = sort(vecP, 'descend');
% sigmaP = vecP(3)/vecP(4);
rs1 = BB*ms;
save('recoveredsignal.mat','rs1','fs');

disp('press any key to hear the 1nd speech'), pause;
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