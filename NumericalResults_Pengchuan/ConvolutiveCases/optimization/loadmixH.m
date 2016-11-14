% called from Convolution_Recovered
% chooses the number q(# of shifts)

function [clean_sig,mixed_x,Xex,q,nonzero_ind_ex,nonzero_Xex]=loadmixH(casenum, hKnown_choice,norm_choice)

p=0.9;
clean_sig1=wavread('data/s3.wav');
clean_sig1=filter([1,-p],1,clean_sig1);

switch casenum
    case 2
         clean_sig2=wavread('data/s7.wav'); 
         clean_sig2=filter([1,-p],1,clean_sig2);
    case 3  
         clean_sig2=resample(wavread('data/s2.wav'),16000,20000);  clean_sig2=clean_sig2./(4*max(abs(clean_sig2)));
         % clean_sig2=resample(wavread('data/songsh'),16000,44100); clean_sig2=clean_sig2(1e5:end);
         % clean_sig2=resample(wavread('data/maleone22'),16000,44100);    
         % clean_sig2=resample(wavread('data/femaleone36'),16000,44100);  
         
%          temp=clean_sig1;
%          clean_sig1=clean_sig2;
%          clean_sig2=temp;
end
   
N=min(size(clean_sig1,1),size(clean_sig2,1));
        
clean_sig1=clean_sig1(1:N);  clean_sig2=clean_sig2(1:N); 
clean_sig1=clean_sig1-mean(clean_sig1);  clean_sig2=clean_sig2-mean(clean_sig2);
	
switch hKnown_choice
    case 1
%		H11 = [1 -0.4  0.1];	H12 = [0.4 0.2 0 ];		 H21 = [-0.3 0.2 0];		H22 = [1 -0.3 0.1];
%		q=3;     % q is the estimated size of weights
% 		H11 = [1 -0.4 0.2 0.1]; H12 = [0 0.4 0.2 -0.1];  H21 = [0 -0.3 0.2 0];		H22 = [1 0.6 -0.3 0.1];
% 		q=4;     % q is the estimated size of weights
		H11 = [1 -0.4 0 0.2 0.1];   H12 = [0 0.4 0.3 0 -0.1];   H21 = [0 -0.3 0 0.2 0];     H22 = [1 0 0.5 -0.2 0.1];
		q=5;     % q is the estimated size of weights
% 		H11 = [1 0.8 0, 0.6 -0.4 0, 0.2 0 -0.1];     H12 = [0 0 -0.6, 0.4 -0.3 0, 0.2 -0.1 0];	
%       H21 = [-0.8, -0.3, 0, 0.3 0.2 0, 0.1 0 0]; 	 H22 = [1 0 0.6, 0 0 0.4, 0 -0.3 0.1];
% 		q=9;     % q is the estimated size of weights
%        H11 = [1 zeros(1,10) -0.4  zeros(1,10) 0.2 0 0];
%		H12 = [zeros(1,6) 0.4 zeros(1,7) -0.2 zeros(1,6) 0.1 0 0 0];
%		H21 = [zeros(1,7) -0.5 zeros(1,6) 0.2 zeros(1,7) 0.1 0 0];
%		H22 = [1 zeros(1,11) -0.3 zeros(1,11) 0.1];
%		q=25;     % q is the estimated size of weights
        
    case 2
		H11 = [1 zeros(1,10) -0.4  zeros(1,10) 0.2 0 0];
		H12 = [zeros(1,6) 0.4 zeros(1,7) -0.2 zeros(1,6) 0.1 0 0 0];
		H21 = [zeros(1,7) -0.5 zeros(1,6) 0.2 zeros(1,7) 0.1 0 0];
		H22 = [1 zeros(1,11) -0.3 zeros(1,11) 0.1];
		q=25;     % q is the estimated size of weights
    case 3
      H11 = [1 zeros(1,24) -0.4 zeros(1,19) 0.2 zeros(1,4)];
		H12 = [zeros(1,20) 0.4 zeros(1,7) -0.2 zeros(1,7) 0.1 zeros(1,13)];
		H21 = [zeros(1,10) 0.5 zeros(1,11) 0.3 zeros(1,11) 0.1 zeros(1,15)];
		H22 = [1 zeros(1,19) -0.3 zeros(1,17) 0.2 zeros(1,11)];
		q=50;   % q is the estimated size of weights
end

if norm_choice>0
   Hnorm1=norm([H11,H21],norm_choice); Hnorm2=norm([H12,H22],norm_choice);
   H11=H11./Hnorm1; H12=H12./Hnorm2; H21=H21./Hnorm1; H22=H22./Hnorm2;
end
		
ss1 = filter(H11,1,clean_sig1) + filter(H12,1,clean_sig2);  %ss1=ss1-mean(ss1); 
% ss1=ss1./max(abs(ss1));
ss2 = filter(H21,1,clean_sig1) + filter(H22,1,clean_sig2);  %ss2=ss2-mean(ss2); 
% ss2=ss2./max(abs(ss2));
% % Remark: after convolution, do not minus the mean, otherwise, the exact solution might change.
mixed_x=[ss1,ss2]';
clean_sig=[clean_sig1,clean_sig2]';

Xex=[H11,H12,H21,H22]'; 

nonzero_ind_ex=[]; nonzero_Xex=[];
x=zeros(4,q);	x(1,:)=H11;   x(2,:)=H12;  x(3,:)=H21;  x(4,:)=H22;
for i=1:4
    [y,ind]=sort(abs(x(i,:))); 
    j=1; while abs(y(j))<1e-6; j=j+1; end 
    nonzero_ind_ex=[nonzero_ind_ex, (i-1)*q + (ind(j:end))];
    nonzero_Xex=[nonzero_Xex, x(i,ind(j:end))];
end    

if casenum==3;  
    NSR=10*log10(norm(clean_sig2,2)/norm(clean_sig1,2))
end

