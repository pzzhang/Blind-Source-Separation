function [BB, detB, P, normP, sigmaP] = postprocess(q, xopt, Xex)
a11=Xex(1:q)'; a12=Xex(q+1:2*q)'; a21=Xex(2*q+1:3*q)'; a22=Xex(3*q+1:4*q)';
b11= xopt(3*q+1:4*q)'; b12=-xopt(q+1:2*q)'; b21=-xopt(2*q+1:3*q)'; b22=xopt(1:q)';
BB = [b11; b12; b21; b22];
detB = det([b11(1) b12(1); b21(1) b22(1)]);

	P=zeros(4,2*q-1);
    normP = zeros(2,2);
	a11=[a11,zeros(1,q-1)];  a12=[a12,zeros(1,q-1)];  a21=[a21,zeros(1,q-1)];  a22=[a22,zeros(1,q-1)];  
    b11=[b11,zeros(1,q-1)];  b12=[b12,zeros(1,q-1)];  b21=[b21,zeros(1,q-1)];  b22=[b22,zeros(1,q-1)]; 
	for tau=1:2*q-1
        P(1,tau)=b11(1:tau)*a11(tau:-1:1)'+b12(1:tau)*a21(tau:-1:1)';
        P(2,tau)=b11(1:tau)*a12(tau:-1:1)'+b12(1:tau)*a22(tau:-1:1)';
        P(3,tau)=b21(1:tau)*a11(tau:-1:1)'+b22(1:tau)*a21(tau:-1:1)';
        P(4,tau)=b21(1:tau)*a12(tau:-1:1)'+b22(1:tau)*a22(tau:-1:1)';
	end
	normP(1,1) = norm(P(1,:));
    normP(1,2) = norm(P(2,:));
    normP(2,1) = norm(P(3,:));
    normP(2,2) = norm(P(4,:));

sigmaP = min(normP(1,1),normP(2,2))/max(normP(1,2),normP(2,1));

end