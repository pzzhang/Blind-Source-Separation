% calculates the energy function of the instantaneous case

demix2ndcorr;

% for i = 1 : 2
%     for j = 1 : 2
%         C(j,j,i) = -C(j,j,i);
%     end
% end
% two variables theta, phi
sigma = 10;
h = 2*pi/100;
theta = 0 : h : 2*pi; % period
% theta = 2*pi; % period
phi = 0 : h : 2*pi;
% phi = 2*pi + 0.01;
N = length(theta);
FF = zeros(N,N);
GG = zeros(N,N);

% optimization method
for i = 1 : N
    for j = 1 : N
        stheta = theta(i);  
        sphi = phi(j);
        u = [cos(stheta); sin(stheta)];
        w = [cos(sphi); sin(sphi)];
        %Fvalue = (u'*C(:,:,1)*w)^2 + (u'*C(:,:,2)*w)^2;
        %Fvalue = (u'*C(:,:,1)*w)^2 + (u'*C(:,:,2)*w)^2 + sigma*(norm(u,1)+norm(w,1));
        Fvalue = (u'*C(:,:,1)*w)^2 + (u'*C(:,:,2)*w)^2 + sigma*((norm(u,1)) + (norm(w,1))); % objective function
        Gvalue = (u'*C(:,:,1)*w)^2 + (u'*C(:,:,2)*w)^2; % objective function w/o penalty terms
        FF(i,j) = Fvalue;
        GG(i,j) = Gvalue;
    end
end
% shreshold G to the accuracy level
Gmax = max(max(GG));
for i = 1 : N
    for j = 1 : N
        if GG(i,j) < eps*Gmax
            GG(i,j) = 0;
        end
    end
end

%load 'angle1'
%analytical method
u = [cos(th(1)); sin(th(1))];
w = [cos(f(1)); sin(f(1))];
%Fvalue = (u'*C(:,:,1)*w)^2 + (u'*C(:,:,2)*w)^2;
Fvaluea = (u'*C(:,:,1)*w)^2 + (u'*C(:,:,2)*w)^2 + sigma*((norm(u,1)) + (norm(w,1)));
Gvaluea = (u'*C(:,:,1)*w)^2 + (u'*C(:,:,2)*w)^2;


% find the optimal point
Fmin = min(min(FF)); % finding the min from both rows & cols
Gmin = min(min(GG));
[row, col] = find(GG == Gmin); %find where min occurs
theta_star = theta(row);
phi_star = phi(col);

% output for theta_star and phi_star
th=theta_star; f=phi_star;
% we have multiple th and f here, and we just try the first one
L = length(th);
for i = 1 : L
    B(:,:,i) = [sin(f(i)),sin(th(i)); cos(f(i)),cos(th(i))];
end

% plot the surface
figure
surf(theta, phi, GG); %3D shaded surface plot
hold on % retains the current graph and adds another graph to it
% scatter3([th(1), theta_star(1)], [f(1), phi_star(1)], [Fvaluea, Fmin], 50, [3, 6]); 
% displays colored circles at the locations specified by the vectors X, Y, and Z
scatter3([th(1), theta_star(1)], [f(1), phi_star(1)], [Gvaluea, Gmin], 50, [3, 6]);
title('Surface of the objective function')
hold off % clears the existing graph and resets axes properties to default

BB = [u(1) -u(2); -w(1) w(2)];
rs2=BB*ms; % estimate sources
figure(3) % plot two received signals
subplot(2,1,1); plot(rs2(1,:))
subplot(2,1,2); plot(rs2(2,:))
title('Two recovered signals');

A1=BB*A0; % approximate inverse times mixing matrix
% fix column permuatation to ensure diagonal dominance
if abs(A1(1,1))+abs(A1(2,2)) < abs(A1(1,2))+abs(A1(2,1))
    A1=[0 1; 1 0]*A1;
end
soundsc(rs2(1,:),fs);
disp('press any key to hear the 2nd speech'), pause;
soundsc(rs2(2,:),fs);

r=norm(A0(1,1)*cs(1,:))/norm(A0(1,2)*cs(2,:))+...
  norm(A0(2,2)*cs(2,:))/norm(A0(2,1)*cs(1,:));
SIRi=10*log10(r); %input SIR (in dB)
r=norm(A1(1,1)*cs(1,:))/norm(A1(1,2)*cs(2,:))+...
  norm(A1(2,2)*cs(2,:))/norm(A1(2,1)*cs(1,:));
SIRo=10*log10(r); % output SIR (in dB)
SIRI=SIRo-SIRi % SIR improvement (in dB) following Eq.(5.20)
