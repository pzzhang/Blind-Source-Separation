% gets called in Convolution_Recovered
% calculates the objective function
function f = CorrelationConv(X)

global A11 A12 A21 A22 q N; 
% B=[-A11,A12;A21,-A22];

% add in the normalized entries.
X=[1;X(1:q-1);X(q:2*q-1);X(2*q:3*q-1);1;X(3*q:4*q-2)];

% x1=[X(3*q+1:4*q);X(q+1:2*q)];       x2=[X(2*q+1:3*q);X(1:q)];
a11=X(1:q)'; a12=X(q+1:2*q)'; a21=X(2*q+1:3*q)'; a22=X(3*q+1:4*q)';
f1 = zeros(N, 1);
for i=1:N 
    f1(i)= ( a12*A21(:,:,i) - a22*A11(:,:,i) )*a21' + ( a22*A12(:,:,i) - a12*A22(:,:,i) )*a11';
end
f=sum(abs(f1.^2));  % value objective function, except the sig(|x|-1) or /|x1||x2| term


% x1=x1./norm(x1,1);        x2=x2./norm(x2,1);

% pow=[2]; %[1/4,1/2,1,2,4];
% xAx=zeros(size(pow));
% for i=1:N
%     xAx = xAx + (abs(x1'*B(:,:,i)*x2)).^pow;
% end

% NN=1;  % NN=N;
% xAx=power(xAx/NN,1./pow);


% it should be f=a(3)^2; 
% f-xAx(3)^2