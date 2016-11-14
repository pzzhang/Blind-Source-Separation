function  [h, mu] = geth(A,g,mu)
% Solve  (Ah + mu*I)h = -g  with possible adjustment of  mu

% Version 04.01.24.  hbn@imm.dtu.dk

% Factorize with check of pos. def.
n = size(A,1);  chp = 1;
while  chp
  [R chp] = chol(A + mu*eye(n));
  if  chp == 0  % check for near singularity
    chp = rcond(R) < 1e-15;
  end
  if  chp,  mu = 10*mu; end
end

% Solve  (R'*R)h = -g
h = R \ (R' \ (-g));   