function [f,J]=fun_GeorgeMarquardt2(x)
% d is the number of measured signals
% x: (d-1)*d angles, representing the recovering matrix
% return the function value f and Jacobian J=Df/Dx

global Cov_v Nshift d; 
% Cov is the correlation matrices for mixed signals: d*d*Nshift 

% the recover matrix
X = reshape(x, d-1, d);
A = zeros(d, d);
for i = 1 : d
    A(:,i) = Ang2Cart(X(:,i), d);
end

% the correlation of the recovered signals
Cov_r = zeros(size(Cov_v));
for n = 1 : Nshift
    Cov_r(:,:,n) = A' * Cov_v(:,:,n) * A;
end

% objective function
inddiag = (eye(d,d) == 1);
indoff = (1 : d^2)';
indoff(inddiag) = [];
D = (d-1)*d;
f = zeros(Nshift*D,1);
for n = 1 : Nshift
    inds = (n-1)*D+1 : n*D;
    Cov_rs = Cov_r(:,:,n);
    f(inds) = Cov_rs(indoff);
end

% Jacobi
if nargout > 1
    
    Jf2a = zeros(Nshift*D, d*d);
    for n = 1 : Nshift
        inds = (n-1)*D+1 : n*D;
        Jf2atemp = zeros(d*d, d*d);
        for i = 1 : d
            indi = (i-1)*d+1 : i*d;
            Jf2atemp(i:d:d^2, indi) = (Cov_v(:,:,n) * A)';
            Jf2atemp(indi, indi) = A' * Cov_v(:,:,n);
        end
        Jf2a(inds,:) = Jf2atemp(indoff,:);
    end
    
    J=zeros(Nshift*D, D);
    for i=1:d
        Jang2cart = JacobiAng2Cart(X(:,i),d);
        indi = (i-1)*d+1 : i*d;
        indires = (i-1)*(d-1)+1 : i*(d-1);
        J(:,indires) = Jf2a(:,indi)*Jang2cart;
    end
end

% % [X, info] = marquardt('fun', [-1.2,1])
% function  [f, J] = fun(x)
% f=[10*(x(2)-x(1)^2);1-x(1)];
% if nargout > 1
%    J=[-20*x(1) 10; -1 0];
% end
