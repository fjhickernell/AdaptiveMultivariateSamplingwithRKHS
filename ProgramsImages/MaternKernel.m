function kerval = MaternKernel(t,x,theta)
[nt,d] = size(t);
[nx,~] = size(x); 
tmx = reshape(t,[nt,1,d]) - reshape(x,[1,nx,d]);
normtmx = reshape(sqrt(sum(tmx.^2,3)),[nt,nx]);
kerval = (1 + theta*normtmx) .*  exp(-theta*normtmx);
