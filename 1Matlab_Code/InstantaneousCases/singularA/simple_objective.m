function y = simple_objective(x, ms)

% old code

Lms = length(ms);
   
theta = x(1);
phi = x(2);

b11 = cos(theta);
b12 = -sin(theta);

b21 = -cos(phi);
b22 = sin(phi);

BB = [b11 b12; b21 b22];
v = BB*ms;

%v1 = b11*ms(1) + b12*ms(2) + b13*ms(3);
%v2 = b21*ms(1) + b22*ms(2) + b23*ms(3);
%v3 = b31*ms(1) + b32*ms(2) + b33*ms(3);

energy = 0;
for shif = 0:10
    energy = energy + (v(1,1:end-shif)*(v(2,shif+1:end)')/(Lms - shif))^2;
end

y = energy;
