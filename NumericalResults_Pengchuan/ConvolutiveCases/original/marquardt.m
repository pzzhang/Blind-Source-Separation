function  [X, info, perf] = marquardt(fun, x0, opts, varargin)
%MARQUARDT  Levenberg-Marquardt's method for least squares.
% Find  xm = argmin{F(x)} , where  x  is an n-vector and
% F(x) = .5 * sum(f_i(x)^2) .
% The functions  f_i(x) (i=1,...,m) and the Jacobian matrix  J(x)  
% (with elements  J(i,j) = Df_i/Dx_j ) must be given by a MATLAB
% function with declaration
%            function  [f, J] = fun(x,p1,p2,...)
% p1,p2,... are parameters of the function.  In connection with 
% nonlinear data fitting they may be arrays with coordinates of 
% the data points.
%  
% Call
%    [X, info] = marquardt(fun, x0)
%    [X, info] = marquardt(fun, x0, opts, p1,p2,...)
%    [X, info, perf] = marquardt(.....)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for  xm .
% opts :  Vector with five elements.  
%         opts(1)  used in starting value for Marquardt parameter: 
%             mu = opts(1) * max(A0(i,i))  with  A0 = J(x0)'*J(x0)
%         opts(2:4)  used in stopping criteria:
%             ||F'||inf <= opts(2)                     or 
%             ||dx||2 <= opts(3)*(opts(3) + ||x||2)    or
%             no. of iteration steps exceeds  opts(4) .
%         opts(5)  lower bound on mu: 
%             mu = opts(5) * max(A(i,i))  with  A = J(x)'*J(x)
%         Default  opts = [1e-3 1e-4 1e-8 100 1e-14]
%         If the input opts has less than 5 elements, it is
%         augmented by the default values.
% p1,p2,..  are passed directly to the function FUN .    
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates
%         columnwise.  Otherwise, computed solution vector.
% info :  Performance information, vector with 6 elements:
%         info(1:4) = final values of 
%             [F(x)  ||F'||inf  ||dx||2  mu/max(A(i,i))] ,
%           where  A = J(x)'*J(x) .
%         info(5) = no. of iteration steps
%         info(6) = 1 : Stopped by small gradient
%                   2 :  Stopped by small x-step
%                   3 :  No. of iteration steps exceeded 
%                  -1 :  x is not a real valued vector
%                  -2 :  f is not a real valued column vector
%                  -3 :  J is not a real valued matrix
%                  -4 :  Dimension mismatch in x, f, J
%                  -5 :  Overflow during computation 
% perf :  Array, holding 
%         perf(1,:) = values of  F(x)
%         perf(2,:) = values of  || F'(x) ||inf
%         perf(3,:) = mu-values.
%
% Method
% Gauss-Newton with Levenberg-Marquardt damping, see eg
% H.B. Nielsen, "Damping parameter in Marquardt's method",
% IMM-REP-1999-05, IMM, DTU, April 1999.

%

% Check parameters and function call
if  nargin < 3,  opts = []; end
opts = checkopts(opts, [1e-3 1e-4 1e-8 100 1e-14]); 

if  nargin < 2,  stop = -1;
else
  [stop x n] = checkx(x0);   
  if  ~stop,  [stop F f J] = checkfJ(fun,x0,varargin{:}); end
end
if  ~stop
  g = J'*f;   ng = norm(g,inf);  A = J'*J;
  if  isinf(ng) | isinf(norm(A(:),inf)),  stop = -5; end
else
  F = NaN;  ng = NaN;
end
if  stop
  X = x0;  perf = [];  info = [F  ng  0  opts(1)  0  stop];
  return
end

%  Finish initialization
mu = opts(1) * max(diag(A));    kmax = opts(4);
Trace = nargout > 2;
if  Trace
  X = repmat(x,1,kmax+1);
  perf = repmat([F; ng; mu],1,kmax+1);
end 
k = 1;   nu = 2;   nh = 0;   stop = 0;

% Iterate
while   ~stop  
  if  ng <= opts(2),  stop = 1;  
  else
    mu = max(mu, opts(5)*max(diag(A)));
    [h mu] = geth(A,g,mu); 
    nh = norm(h);   nx = opts(3) + norm(x);
    if  nh <= opts(3)*nx,  stop = 2; end 
  end 
  if  ~stop
    xnew = x + h;   h = xnew - x;   dL = (h'*(mu*h - g))/2; 
    [stop Fn fn Jn] = checkfJ(fun, xnew, varargin{:});  
    if  ~stop      
      k = k + 1;  dF = F - Fn;  
      if  Trace
        X(:,k) = xnew;   perf(:,k) = [Fn norm(Jn'*fn,inf) mu]'; end
      if  (dL > 0) & (dF > 0)               % Update x and modify mu
        x = xnew;   F = Fn;  J = Jn;  f = fn; 
        A = J'*J;   g = J'*f;   ng = norm(g,inf);
        mu = mu * max(1/3, 1 - (2*dF/dL - 1)^3);   nu = 2;
      else                                  % Same  x, increase  mu
        mu = mu*nu;  nu = 2*nu; 
      end
      if      k > kmax,                           stop = 3; 
      elseif  isinf(ng) | isinf(norm(A(:),inf)),  stop = -5; end
    end
  end   
end
%  Set return values
if  Trace
  X = X(:,1:k);   perf = perf(:,1:k);
else,  X = x;  end
if  stop < 0,  F = NaN;  ng = NaN; end
info = [F  ng  nh  mu/max(diag(A))  k-1  stop];