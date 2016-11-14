function [ V ] = JointDiag( epsilon, M )
%Joint Diagonalize M
%   M is n*n*K; every M(:,:,k) is real symmetric

% get dimensions
[n,~,K] = size(M);

% initialize V
V = eye(n);

% Jacobi iteration
L = zeros(K,2);
[yoff, ytotal] = Moff ( M, n, K );
nsweep = 0;
while (yoff > epsilon*ytotal && nsweep < 30)
    for p = 1 : n
        for q = p+1 : n
            L(:,1) = squeeze(M(p,q,:));
            L(:,2) = squeeze(M(q,q,:)-M(p,p,:))/2;
            [~,~,W] = svd(L,'econ');
            w = W(:,2);
            if w(1) < 0
                w = -w;
            end
            c = sqrt((1+w(1))/2);
            s = w(2)/(2*c);
            R = [c, -s; s, c];
            V(:,[p,q]) = V(:,[p,q])*R;
            for k = 1 : K
                M(:,[p,q],k) = M(:,[p,q],k)*R;
                M([p,q],:,k) = R'*M([p,q],:,k);
            end
        end
    end
    [yoff, ytotal] = Moff ( M, n, K );
    nsweep = nsweep + 1;
end

end

function [ yoff, ytotal ] = Moff ( M, n, K )
% compute off diagonal values and total values
yfro2 = zeros(1,K);
ydiag2 = zeros(1,K);

for k = 1 : K
    yfro2(k) = norm(M(:,:,k),'fro')^2;
    ydiag2(k) = norm(diag(M(:,:,k)))^2;
end

ytotal = sum(yfro2);
yoff = ytotal - sum(ydiag2);

end