function [Kmat, Kdateval, Kdiageval, errKNull, Pdat, Peval, PTKinv, PTKinvP, Mmat] = ...
KMP(xdata, xeval, kernel, trend)
Kmat = kernel(xdata,xdata);
if nargout > 1
   [Kdateval, Kdiageval, errKNull] = kernel(xdata,xeval);
end

if nargin >= 4
   Pdat = trend(xdata);
   Peval = trend(xeval);
   PTKinv = Pdat'/Kmat;
   PTKinvP = PTKinv*Pdat;
   Mmat = inv(Kmat) - PTKinv'*(PTKinvP\PTKinv);
end
