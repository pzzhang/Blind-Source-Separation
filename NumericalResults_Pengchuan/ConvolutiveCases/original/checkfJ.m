function  [err, F,f,J] = checkfJ(fun,x,varargin)
%CHECKFJ  Check Matlab function which is called by a 
% nonlinear least squares solver.

% Version 04.01.25.  hbn@imm.dtu.dk

err = 0;   F = NaN;  n = length(x);
if  nargout > 3    % Check  f  and  J
  [f J] = feval(fun,x,varargin{:});
  sf = size(f);   sJ = size(J);
  if  sf(2) ~= 1 | ~isreal(f) | any(isnan(f(:))) | any(isinf(f(:)))
    err = -2;  return, end
  if  ~isreal(J) | any(isnan(J(:))) | any(isinf(J(:)))
    err = -3;  return, end
  if  sJ(1) ~= sf(1) | sJ(2) ~= n
    err = -4;  return, end
  
else  % only check  f
  f = feval(fun,x,varargin{:});
  sf = size(f);   
  if  sf(2) ~= 1 | ~isreal(f) | any(isnan(f(:))) | any(isinf(f(:)))
    err = -2;  return, end
end

% Objective function
F = (f'*f)/2;
if  isinf(F),  err = -5; end