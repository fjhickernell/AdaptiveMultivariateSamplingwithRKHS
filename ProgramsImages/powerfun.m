function [powkeval,powpeval] = powerfun(xeval,xdata,kernel,kerneldiag,trend)
Kmat = kernel(xdata,xdata);
Kdateval = kernel(xdata,xeval);
powkeval = sqrt(kerneldiag(xeval) ...
   - sum(Kdateval .* (Kmat\Kdateval),1)');
if nargin < 5
   powpeval = [];
else
   Pmat = trend(xdata);
   Peval = trend(xeval);
   PTKinv = Pmat'/Kmat;
   PTKinvKdateval = PTKinv*Kdateval;
   powpeval = sqrt(sum((Peval-PTKinvKdateval').^2,2));
end

   
