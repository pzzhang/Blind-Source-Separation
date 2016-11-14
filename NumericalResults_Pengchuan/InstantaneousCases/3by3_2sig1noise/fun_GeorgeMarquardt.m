function [f,J]=fun_GeorgeMarquardt(x)
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
Cov_r_rev = zeros(size(Cov_v));
for n = 1 : Nshift
    Cov_r(:,:,n) = A' * Cov_v(:,:,n) * A;
    Cov_r_rev(:,:,n) = Cov_r(:,:,n) - diag(diag(Cov_r(:,:,n)));
end

% objective function
f = 0;
for n = 1 : Nshift
    f = f + sum(sum(Cov_r_rev(:,:,n).^2));
end

% gradient
if nargout > 1
    J=zeros((d-1)*d, 1);
    for i=1:d
        indtemp = (i-1)*(d-1)+1 : i*(d-1);
        Jang2cart = JacobiAng2Cart(X(:,i),d);
        Gf2cart = zeros(1, d);
        for n = 1 : Nshift
            Jtemp = Cov_v(:,:,n) * A * diag(Cov_r_rev(i,:,n));
            Jtemp = Jtemp' + diag(Cov_r_rev(:,i,n))*(A')*Cov_v(:,:,n);
            Gf2cart = Gf2cart + sum(Jtemp);
        end
        J(indtemp) = 2*Gf2cart*Jang2cart;
    end
end

% % [X, info] = marquardt('fun', [-1.2,1])
% function  [f, J] = fun(x)
% f=[10*(x(2)-x(1)^2);1-x(1)];
% if nargout > 1
%    J=[-20*x(1) 10; -1 0];
% end
