function f=objfunval(x,C,norm_choice,objfun_choice)
global A11 A12 A21 A22 N q; 
% Aij(k,m,n)=E(y_i(t-k),y_j(t-m-n))   
% N+2 is the number of equations and q is number of weights
f=zeros(N,1);
a11=x(1:q)'; a12=x(q+1:2*q)'; a21=x(2*q+1:3*q)'; a22=x(3*q+1:4*q)';
for i=1:N
    f(i)= ( a12*A21(:,:,i) - a22*A11(:,:,i) )*a21' + ( a22*A12(:,:,i) - a12*A22(:,:,i) )*a11';
end

if norm_choice>0
    x1norm=norm([a12,a22],norm_choice)^norm_choice;
    x2norm=norm([a11,a21],norm_choice)^norm_choice;
    if objfun_choice==1
        f=[C*x1norm-C; C*x2norm-C; f];
        f=sum(f.^2);
    elseif objfun_choice==0  % obj function is quotient 
        f=sum(f.^2);
        f=f./(x1norm*x2norm)^2;
    end
else  % a_11(1)=a_22(1)=1
    f=sum(f.^2);
end