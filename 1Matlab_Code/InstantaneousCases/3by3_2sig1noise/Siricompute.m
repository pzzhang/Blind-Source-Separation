function [ SIRI ] = Siricompute( cs, A0, P )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

P0 = zeros(3,3);
for i = 1 : 3
    P0(i,:) = P(i,:)/max(abs(P(i,:)));
    for j = 1 : 3
        if abs(P0(i,j)) > 0.5
            P0(i,j) = 1;
        else
            P0(i,j) = 0;
        end
    end
end

csP = P0*cs;
P = P*P0';

r=norm(A0(1,1)*cs(1,:))/(norm(A0(1,2)*cs(2,:))+norm(A0(1,3)*cs(3,:))/2)+...
  norm(A0(2,2)*cs(2,:))/(norm(A0(2,1)*cs(1,:))+norm(A0(2,3)*cs(3,:))/2)+...
  norm(A0(3,3)*cs(3,:))/(norm(A0(3,1)*cs(1,:))+norm(A0(3,2)*cs(2,:))/2);
SIRi=10*log10(r); %input SIR (in dB)
r=norm(P(1,1)*csP(1,:))/(norm(P(1,2)*csP(2,:))+norm(P(1,3)*csP(3,:))/2)+...
  norm(P(2,2)*csP(2,:))/(norm(P(2,1)*csP(1,:))+norm(P(2,3)*csP(3,:))/2)+...
  norm(P(3,3)*csP(3,:))/(norm(P(3,1)*csP(1,:))+norm(P(3,2)*csP(2,:))/2);
SIRo=10*log10(r); % output SIR (in dB)
SIRI=SIRo-SIRi; % SIR improvement (in dB) following Eq.(5.20)

end

