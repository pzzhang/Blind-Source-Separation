function [ BB, vs ] = joint_diag_signal( ms, epsilon )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[d, L] = size(ms);
% compute the covariance matrix with shifts
K = 10;
C = zeros(d,d,K);
for k = 0:K-1 % compute correlation with shifts
    shif = k;
    C(:,:,k+1)=ms(:,1:end-shif)*(ms(:,shif+1:end)')./(L-shif);
end

A = JointDiag( epsilon, C );

BB = transpose(A);
vs = BB*ms;

end

