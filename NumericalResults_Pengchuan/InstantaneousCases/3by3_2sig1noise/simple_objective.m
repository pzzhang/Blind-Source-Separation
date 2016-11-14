function y = simple_objective(x, ms)

% old code

Lms = length(ms);
   
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
v = BB*ms;

%v1 = b11*ms(1) + b12*ms(2) + b13*ms(3);
%v2 = b21*ms(1) + b22*ms(2) + b23*ms(3);
%v3 = b31*ms(1) + b32*ms(2) + b33*ms(3);

energy = 0;
for shif = 0:10
    energy = energy + (v(1,1:end-shif)*(v(2,shif+1:end)')/(Lms - shif))^2 ...
            + (v(1,1:end-shif)*(v(3,shif+1:end)')/(Lms - shif))^2 ...
            + (v(2,1:end-shif)*(v(3,shif+1:end)')/(Lms - shif))^2;
end

y = energy;
