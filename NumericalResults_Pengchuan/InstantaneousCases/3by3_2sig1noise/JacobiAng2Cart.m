function [ J ] = JacobiAng2Cart( x, d )
%Jacobi matrix at x for Ang2Cart
%   J is a d*(d-1) matrix

if d == 2
    J = [-sin(x); cos(x)];
else if d > 2
        am1 = Ang2Cart(x(2:end), d-1);
        Jm1 = JacobiAng2Cart( x(2:end), d-1 );
        J = zeros(d, d-1);
        J(1,1) = -sin(x(1));
        J(2:end,1) = cos(x(1))*am1;
        J(2:end,2:end) = sin(x(1))*Jm1;
    else
        error('dimension should not be less than 2!');
    end
end


end

