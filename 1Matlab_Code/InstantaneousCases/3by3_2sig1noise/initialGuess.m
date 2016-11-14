function [ xinit ] = initialGuess( x,y )
%We want to find the angles to change from x to y
%   y = A x
%   A = [cos(theta)cos(phi) sin(theta)cos(phi) sin(phi);...]

% estimate A
A = y * x';

% estimate angles
angles = zeros(2,3);
for i = 1 : 3
    thetatemp = atan2(A(i,2),A(i,1));
    phitemp = atan2(A(i,3), sqrt(A(i,1)^2+A(i,2)^2));
    angles(1,i) = thetatemp;
    angles(2,i) = phitemp;
end

xinit = reshape(angles, 1, 6);


end

