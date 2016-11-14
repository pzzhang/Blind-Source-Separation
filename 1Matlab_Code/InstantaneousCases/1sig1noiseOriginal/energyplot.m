%% problem setting
clear all;
close all;
L=36000; % set uniform signal length
[sig1,fs]=wavread('data/s3.wav'); % load wav files; returns the sample rate (Fs) in Hertz used to encode the data 
cs(1,1:L)=sig1(1:L)/norm(sig1(1:L));
% [sig2,fs]=wavread('data/s1.wav'); % load wav files; returns the sample rate (Fs) in Hertz used to encode the data 
% cs(2,1:L)=sig2(1:L)/norm(sig2(1:L));
load 'cr.mat';
cs(2,1:L) = cr;
cs(2,1:L) = cs(2,1:L)/norm(cs(2,1:L));
% Cstable conditions
Lcs=length(cs);
C0=zeros(2,2,2); Cs_shift = zeros(2,2);
for shif=0:1 % compute correlation with shifts
    C0(:,:,shif+1)=cs(:,1:end-shif)*(cs(:,shif+1:end)')./(L-shif);
    Cs_shift(shif+1,:) = diag(C0(:,:,shif+1))';
end
Cs = norm(inv(Cs_shift),2);

% regular case
% theta0 = pi/4;
% phi0 = - pi/4;
% degenerate case
theta0 = pi/4;
phi0 = pi/4 + 1e-3;

x0 = [theta0, phi0];
A0 = [cos(x0(1)), cos(x0(2)); sin(x0(1)), sin(x0(2))];
% Asigma assumpution
condA = cond(A0);
ms=A0*cs;

% mesh the two variables
h = 2*pi/100;
theta = 0 : h : 2*pi; % period
% theta = 2*pi; % period
phi = 0 : h : 2*pi;
% phi = 2*pi + 0.01;
N = length(theta);
GG = zeros(N,N);

% optimization method
for i = 1 : N
    for j = 1 : N
        stheta = theta(i);  
        sphi = phi(j);
        u = [cos(stheta - phi0); cos(stheta - theta0)];
        w = [cos(sphi - phi0); cos(sphi - theta0)];
        %Fvalue = (u'*C(:,:,1)*w)^2 + (u'*C(:,:,2)*w)^2;
        %Fvalue = (u'*C(:,:,1)*w)^2 + (u'*C(:,:,2)*w)^2 + sigma*(norm(u,1)+norm(w,1));
        Gvalue = (u'*C0(:,:,1)*w)^2 + (u'*C0(:,:,2)*w)^2; % objective function w/o penalty terms
        GG(i,j) = Gvalue;
    end
end
% shreshold G to the accuracy level
% Gmax = max(max(GG));
% for i = 1 : N
%     for j = 1 : N
%         if GG(i,j) < eps*Gmax
%             GG(i,j) = 0;
%         end
%     end
% end


% find the optimal point
Gmin = min(min(GG));
Gmax = max(max(GG));
[row, col] = find(GG == Gmin); %find where min occurs
theta_star = theta(row);
phi_star = phi(col);


% plot the surface
figure
surf(theta, phi, GG); %3D shaded surface plot
hold on % retains the current graph and adds another graph to it
% scatter3([th(1), theta_star(1)], [f(1), phi_star(1)], [Fvaluea, Fmin], 50, [3, 6]); 
% displays colored circles at the locations specified by the vectors X, Y, and Z
scatter3([theta_star(1)], [phi_star(1)], [Gmin], 50, [1,0,0], 'filled');
% axis([0 2*pi 0 2*pi])
axis([0 2*pi 0 2*pi 0 0.1*Gmax])
title('objective function before QR')
hold off % clears the existing graph and resets axes properties to default

figure
contour(theta, phi, GG, 20); %3D shaded surface plot
hold on % retains the current graph and adds another graph to it
% scatter3([th(1), theta_star(1)], [f(1), phi_star(1)], [Fvaluea, Fmin], 50, [3, 6]); 
% displays colored circles at the locations specified by the vectors X, Y, and Z
scatter([theta_star(1)], [phi_star(1)], 50, [1,0,0],'filled');
axis([0 2*pi 0 2*pi])
title('Energy Contour after QR')
hold off % clears the existing graph and resets axes properties to default