function [ BB, vs, detB, fval,exitflag,output ] = demix_min_cov_3d( ms, x0qr )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    % optimization setting
    ObjectiveFunction = @(x) simple_objective(x, ms);
    lb = -pi/2*[2,1,2,1,2,1];
    ub = pi/2*[2,1,2,1,2,1];
    fvalref = ObjectiveFunction(x0qr);
    fval = 10;
    trial = 0;
    while (fval > fvalref) && (trial <= 20)
    trial = trial + 1;
    % initial guess
    xperturb = 0.5;
    X0 = x0qr + xperturb*pi/4*(2*rand(1,6)-1);
    for i = 1 : 3
        X0(2*i-1) = X0(2*i-1)-2*pi*floor((X0(2*i-1)+pi)/(2*pi));
        X0(2*i) = X0(2*i)-pi*floor((X0(2*i)+pi/2)/(pi));
    end
    % optimization solver
    options = optimset('MaxFunEvals',10^6, 'MaxIter', 10^7, 'TolFun',1e-13, 'TolX', 1e-13);
    [x,fval,exitflag,output] = fminsearchcon(ObjectiveFunction,X0,lb,ub,[],[],[],options);
    end
    
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
detB = det(BB);
vs = BB*ms;

end

