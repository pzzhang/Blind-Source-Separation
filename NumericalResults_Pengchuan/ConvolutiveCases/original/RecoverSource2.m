function [v,s]=RecoverSource2(mixed_x,X,q,str,show_s_choice)
a11=X(1:q)'; a12=X(q+1:2*q)'; a21=X(2*q+1:3*q)'; a22=X(3*q+1:4*q)';

v=zeros(size(mixed_x));
v(1,:) = filter( a22,1,mixed_x(1,:)) + filter(-a12,1,mixed_x(2,:));  
v(2,:) = filter(-a21,1,mixed_x(1,:)) + filter( a11,1,mixed_x(2,:));  

if show_s_choice==1
	h=zeros(1,2*q-1);
	a11=[a11,zeros(1,q-1)];  a12=[a12,zeros(1,q-1)];  a21=[a21,zeros(1,q-1)];  a22=[a22,zeros(1,q-1)];  
	for tau=1:2*q-1
        h(tau)=a11(1:tau)*a22(tau:-1:1)'-a12(1:tau)*a21(tau:-1:1)';
	end
	h=real(ifft(1./fft(h)));
	s=zeros(size(v));
	for i=1:2; s(i,:) = filter( h,1,v(i,:));  end
end
show_play(v,2,['v,  ',str]); 
if show_s_choice==1; show_play(s,2,['s,  ',str]); end

