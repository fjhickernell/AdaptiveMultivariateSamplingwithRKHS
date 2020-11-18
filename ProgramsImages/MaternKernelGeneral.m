function [kerval,kdiageval,errKNull,theta] = ...
   MaternKernelGeneral(t,x,theta,transYes)
   nugget = 1e-10;
   if nargin < 4, transYes = false; end
   [nt,d] = size(t);
   [nx,dx] = size(x);
   dth = length(theta);
   if ~(dx == d) || ~((mod(dth-1, d) == 0) || (dth <= 2))
      error(['dim of t = ' int2str(d) ...
      ', dim of x = ' int2str(dx) ...
      ', and dim of theta = ' int2str(dth) ...
      ' do not all match'])
   end
   if (dth <= d + 1)
      if dth == 1
         thetascale = 0;
         dthp = 1;
      else
         thetascale = theta(2:dth);
         dthp = dth-1;
      end
      if transYes
         theta = exp(theta);
      end
      nu = theta(1);
      sc = sqrt(2*nu);
      tmx = (reshape(t,[nt,1,d]) - reshape(x,[1,nx,d])) ...
         .* reshape(thetascale*sc,[1 1 dthp]);
      normtmx = reshape(sqrt(sum(tmx.^2,3)),[nt,nx]);
      if nu < 40
         kerval = (2^(1-nu)/gamma(nu)) .* (normtmx.^nu) .* besselk(nu,normtmx);
         kerval(normtmx==0) = 1;
      else
         kerval = exp(-normtmx.^2/2);
      end
      if ((nt == nx) && (rcond(kerval) < nugget)) %add nugget
         kerval(1:nt+1:nt^2) = kerval(1:nt+1:nt^2) + nugget;
      end
      if nargout > 1
         kdiageval = ones(nx,1) + nugget;
         errKNull = max(kdiageval);
      end
   else
      dd = (dth-1)/2;
      if transYes
         theta = [exp(theta(1:dd+1)) theta(dd+2:dth)];
      end
      nu = theta(1);
      sc = sqrt(2*nu);
      thetaa = theta(2:dd+1);
      thetab = theta(dd+2:dth);
      tmx = (reshape(t,[nt,1,d]) - reshape(x,[1,nx,d])) ...
         .* reshape(thetaa*sc,[1 1 dd]);
      normtmx = reshape(sqrt(sum(tmx.^2,3)),[nt,nx]);
      tmpl = reshape(sum((reshape(t,[nt,1,d]) + reshape(x,[1,nx,d]))...
         .* reshape(thetab,[1 1 dd]),3),[nt,nx]);
      if nu < 40
         kerval = (2^(1-nu)/gamma(nu)) .* (normtmx.^nu) .* besselk(nu,normtmx);
         kerval(normtmx==0) = 1;
      else
         kerval = exp(-normtmx.^2/2);
      end
      kerval = exp(tmpl) .* kerval; 
      if nargout > 1
         kdiageval = exp(sum(x.*(2*thetab),2));
         errKNull = max(kdiageval);
      end
   end
   if any(any(isnan(kerval))) || ((nt == nx) && (rcond(kerval) < 1e-12))
      %keyboard
   end
end


