function [f,J]=fun_Quotientmarquardt(x,norm_choice)
% for F=(f_1,...,f_{N}), return the function value f_i and Jacobian J=Df_i/Dx_j
global A11 A12 A21 A22 N q; 
% Aij(k,m,n)=E(y_i(t-k),y_j(t-m-n))   
% N+2 is the number of equations and q is number of weights
f=zeros(N,1);
a11=x(1:q)'; a12=x(q+1:2*q)'; a21=x(2*q+1:3*q)'; a22=x(3*q+1:4*q)';
x1norm=norm([a12,a22],norm_choice)^norm_choice;
x2norm=norm([a11,a21],norm_choice)^norm_choice;

for i=1:N
    f(i)= ( a12*A21(:,:,i) - a22*A11(:,:,i) )*a21' + ( a22*A12(:,:,i) - a12*A22(:,:,i) )*a11';
end

f=f./(x1norm*x2norm);

if nargout > 1
    J=zeros(4*q,N);
    for j=1:q
        J(j,:)    =a22*squeeze(A12(:,j,:)) - a12*squeeze(A22(:,j,:));
        J(q+j,:)  =a21*squeeze(A21(j,:,:)) - a11*squeeze(A22(j,:,:));
        J(2*q+j,:)=a12*squeeze(A21(:,j,:)) - a22*squeeze(A11(:,j,:));
        J(3*q+j,:)=a11*squeeze(A12(j,:,:)) - a21*squeeze(A11(j,:,:));
    end
    J=J'./(x1norm*x2norm);
    switch norm_choice
        case 1
            eps=1e-16;
            temp=x./(abs(x)+eps);
        case 2
            temp=2*x;
	end
    temp=temp';
    temp([1:q,2*q+1:3*q])=temp([1:q,2*q+1:3*q])./(x1norm*x2norm^2);
    temp([q+1:2*q,3*q+1:4*q])=temp([q+1:2*q,3*q+1:4*q])./(x1norm^2*x2norm);
    J=J-f*temp;
end
