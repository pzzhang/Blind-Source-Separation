function [ a ] = Ang2Cart(x, d)
% x is (d-1) Euler angles
%   A is the corresponding unit vectors.

if d == 2
    a = [cos(x); sin(x)];
else if d > 2
        am1 = Ang2Cart(x(2:end), d-1);
        a = [cos(x(1)); sin(x(1))*am1];
    else
        error('dimension should not be less than 2!');
    end
end


end

