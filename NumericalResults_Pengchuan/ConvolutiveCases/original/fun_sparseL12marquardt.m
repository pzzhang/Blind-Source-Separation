function [f,J]=fun_sparseL12marquardt(nonzero_x,C,norm_choice,nonzero_ind)
% for F=(f_1,...,f_{N+1},f_{N+2}), return the function value f_i and Jacobian J=Df_i/Dx_j
% the input nonzerox is only the non-zero element

global A11 A12 A21 A22 N q; 
% Aij(k,m,n)=E(y_i(t-k),y_j(t-m-n))   
% N+2 is the number of equations and q is number of weights before considering the sparseness
f=zeros(N,1);
x=zeros(4*q,1);
x(nonzero_ind)=nonzero_x;
a11=x(1:q)'; a12=x(q+1:2*q)'; a21=x(2*q+1:3*q)'; a22=x(3*q+1:4*q)';
for i=1:N
    f(i)= ( a12*A21(:,:,i) - a22*A11(:,:,i) )*a21' + ( a22*A12(:,:,i) - a12*A22(:,:,i) )*a11';
end

f=[C*norm([a11,a21],norm_choice)^norm_choice-C; C*norm([a12,a22],norm_choice)^norm_choice-C; f];

if nargout > 1
    J=zeros(4*q-2,N);
    for j=1:q
        J(j,:)    =a22*squeeze(A12(:,j,:)) - a12*squeeze(A22(:,j,:));
        J(q+j,:)  =a21*squeeze(A21(j,:,:)) - a11*squeeze(A22(j,:,:));
        J(2*q+j,:)=a12*squeeze(A21(:,j,:)) - a22*squeeze(A11(:,j,:));
        J(3*q+j,:)=a11*squeeze(A12(j,:,:)) - a21*squeeze(A11(j,:,:));
    end
    J=J';
end

switch norm_choice
    case 1
        eps=1e-16;
        temp=C*x./(abs(x)+eps);
    case 2
        temp=(2*C)*x;
end
temp=temp';
J=[temp(1:q),zeros(1,q),temp(2*q+1:3*q),zeros(1,q); ...
   zeros(1,q),temp(q+1:2*q),zeros(1,q),temp(3*q+1:4*q); J];

J=J(:,nonzero_ind);
