function [ xinit ] = initialGuess2D( x,y )
%We want to find the angles to change from x to y
%   y = A x
%   A = [cos(theta)cos(phi) sin(theta)cos(phi) sin(phi);...]

% estimate A
A = y * x';

% estimate angles
xinit = [atan(A(1,1)/A(2,1)), atan(A(1,2)/A(2,2))];


end