function [ SIRI ] = Siricompute2D( cs, A0, A1 )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
if abs(A1(1,1))+abs(A1(2,2)) < abs(A1(1,2))+abs(A1(2,1)) 
    A1=[0 1; 1 0]*A1;
end
r=min(norm(A0(1,1)*cs(1,:))/norm(A0(1,2)*cs(2,:)), norm(A0(2,2)*cs(2,:))/norm(A0(2,1)*cs(1,:)));
SIRi=10*log10(r); %input SIR (in dB)
r=min(norm(A1(1,1)*cs(1,:))/norm(A1(1,2)*cs(2,:)), norm(A1(2,2)*cs(2,:))/norm(A1(2,1)*cs(1,:)));
SIRo=10*log10(r); % output SIR (in dB)
SIRI=SIRo-SIRi; % SIR improvement (in dB) following Eq.(5.20)

end

