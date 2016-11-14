function  opts = checkopts(opts, default)
%CHECKOPTS  Replace illegal values by default values.

% Version 04.01.25.  hbn@imm.dtu.dk

a = default;  la = length(a);  lo = length(opts);
for  i = 1 : min(la,lo)
  oi = opts(i);
  if  isreal(oi) & ~isinf(oi) & ~isnan(oi) & oi > 0
    a(i) = opts(i); 
  end
end
if  lo > la,  a = [a 1]; end % for linesearch purpose
opts = a;