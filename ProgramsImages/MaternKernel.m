function [kerval,kdiageval,errKNull] = MaternKernel(t,x,theta)
    [nt,d] = size(t);
    [nx,dx] = size(x);
    dth = length(theta);
    if ~(dx == d) || ~((mod(dth, d) == 0) || (dth == 1))
        error(['dim of t = ' int2str(d) ...
        ', dim of x = ' int2str(dx) ...
        ', and dim of theta = ' int2str(dth) ...
        ' do not all match'])
    end
    
    if (dth == d)
       tmx = (reshape(t,[nt,1,d]) - reshape(x,[1,nx,d])) ...
           .* reshape(theta,[1 1 dth]);
       normtmx = reshape(sqrt(sum(tmx.^2,3)),[nt,nx]);
       kerval = (1 + normtmx) .*  exp(-normtmx);
       if nargout > 1
          kdiageval = ones(size(x,1),1);
          errKNull = max(kdiageval);
       end
    else
       dd = dth/2;
       thetaa = theta(1:dd);
       thetab = theta(dd+1:dth);
       tmx = (reshape(t,[nt,1,d]) - reshape(x,[1,nx,d])) ...
           .* reshape(thetaa,[1 1 dd]);
       normtmx = reshape(sqrt(sum(tmx.^2,3)),[nt,nx]);
       tmpl = reshape(sum((reshape(t,[nt,1,d]) + reshape(x,[1,nx,d]))...
           .* reshape(thetab,[1 1 dd]),3),[nt,nx]);
       kerval = exp(tmpl).*(1 + normtmx) .*  exp(-normtmx); 
       kdiageval = exp(x*(2*thetab'));
       if nargout > 1
          kdiageval = ones(size(x,1),1);
          errKNull = max(kdiageval);
       end
    end
end

