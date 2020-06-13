function [kerval,kerneldiag,errKNull] = MaternKernel(t,x,theta)
[nt,d] = size(t);
[nx,dx] = size(x);
dth = length(theta);
if ~(dx == d) || ~((dth == d) || (dth == 1))
   error(['dim of t = ' int2str(d) ...
      ', dim of x = ' int2str(dx) ...
      ', and dim of theta = ' int2str(dth) ...
      ' do not all match'])
end
if nt > 0
   tmx = (reshape(t,[nt,1,d]) - reshape(x,[1,nx,d])) ...
      .* reshape(theta,[1 1 dth]);
   normtmx = reshape(sqrt(sum(tmx.^2,3)),[nt,nx]);
   kerval = (1 + normtmx) .*  exp(-normtmx);
else
   kerval = [];
   kerneldiag = @(x) ones(size(x,1),1);
   errKNull = 1;
end
