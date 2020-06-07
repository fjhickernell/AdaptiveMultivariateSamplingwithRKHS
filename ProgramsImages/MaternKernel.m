function kerval = MaternKernel(t,x,theta)
[nt,d] = size(t);
[nx,dx] = size(x);
dth = length(theta);
if ~(dx == d) || ~((dth == d) || (dth == 1))
   error(['dim of t = ' int2str(d) ...
      ', dim of x = ' int2str(dx) ...
      ', and dim of theta = ' int2str(dth) ...
      ' do not all match'])
end
tmx = (reshape(t,[nt,1,d]) - reshape(x,[1,nx,d])) ...
   .* reshape(theta,[1 1 dth]);
normtmx = reshape(sqrt(sum(tmx.^2,3)),[nt,nx]);
kerval = (1 + normtmx) .*  exp(-normtmx);
