function [kerval,kdiageval,errKNull] = ...
   GaussKernel(t,x,theta)
   [nt,d] = size(t);
   nx = size(x,1);
   tmx = theta*(reshape(t,[nt,1,d]) - reshape(x,[1,nx,d]));
   normtmx2 = reshape(sum(tmx.^2,3),[nt,nx]);
   kerval = exp(-normtmx2);
   if nargout > 1
      kdiageval = ones(size(x,1),1);
      errKNull = max(kdiageval);
   end
end

