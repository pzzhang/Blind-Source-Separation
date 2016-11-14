function [ BB, rs1, detB ] = demix_analytic( ms, zeroflag )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

%% optimization setting
Lms=length(ms);
C=zeros(2,2,2); d=zeros(2,1); f=zeros(2,1); A=zeros(2,2); B=zeros(2,2,2);
if zeroflag == 0
    for shif=1:2 % compute correlation with shifts
        C(:,:,shif)=ms(:,1:end-shif)*(ms(:,shif+1:end)')./(Lms-shif);
    end
else
    for shif=0:1 % compute correlation with shifts
        C(:,:,shif+1)=ms(:,1:end-shif)*(ms(:,shif+1:end)')./(Lms-shif);
    end
end
a=C(2,1,1)*C(2,2,2)-C(2,2,1)*C(2,1,2);
b=C(2,2,1)*C(1,1,2)+C(1,2,1)*C(2,1,2)-C(2,1,1)*C(1,2,2)-C(1,1,1)*C(2,2,2);
c=C(1,1,1)*C(1,2,2)-C(1,2,1)*C(1,1,2);
if b^2 >= 4*a*c
th(1)=(-b+sqrt(b^2-4*a*c))/(2*a); % define two tan theta values 
th(2)=(-b-sqrt(b^2-4*a*c))/(2*a);
else
    th(1)=-b/(2*a); th(2)=th(1);
end
for i=1:2
    f(i)=(C(1,1,1)-th(i)*C(2,1,1))/(C(1,2,1)-th(i)*C(2,2,1)); % computing tan phi 
end

th=atan(th); f=atan(f);
for i=1:2
    B(:,:,i)=[sin(f(i)),sin(th(i)); cos(f(i)),cos(th(i))]; % a11=sin(phi),a12=sin(theta),...
end
BB=inv(B(:,:,2)); % BB is the inverse of B
detB = det(BB);
rs1=BB*ms; % estimate sources

end

