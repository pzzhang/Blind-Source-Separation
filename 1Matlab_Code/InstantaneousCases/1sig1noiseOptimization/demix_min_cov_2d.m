function [ BB, vs, detB, fval,exitflag,output ] = demix_min_cov_2d( ms, x0qr )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    % optimization setting
    ObjectiveFunction = @(x) simple_objective(x, ms);
    lb = -pi/2*[1,1];
    ub = pi/2*[1,1];
    fvalref = ObjectiveFunction(x0qr);
    fval = 1;
    while fval > fvalref
    % initial guess
    xperturb = 1;
    X0 = x0qr + xperturb*pi/4*(2*rand(1,2)-1);
    for i = 1 : 2
        X0(i) = pi*(X0(i)/(pi)-floor(X0(i)/(pi)))-pi/2;
    end
    % optimization solver
    options = optimset('MaxFunEvals',10^6, 'MaxIter', 10^7, 'TolFun',1e-12);
    [x,fval,exitflag,output] = fminsearchcon(ObjectiveFunction,X0,lb,ub,[],[],[],options);
    end
    
% form BB and vs
theta = x(1);
phi = x(2);

b11 = cos(theta);
b12 = -sin(theta);

b21 = -cos(phi);
b22 = sin(phi);

BB = [b11 b12; b21 b22];
detB = det(BB);
vs = BB*ms;

end

