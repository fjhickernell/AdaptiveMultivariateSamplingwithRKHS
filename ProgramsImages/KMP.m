function [Kmat, Kdateval, Kdiageval, Pdat, Peval, PTKinv, PTKinvP, Mmat] = ...
   KMP(xdata, xeval, kernel, kerneldiag, trend)
Kmat = kernel(xdata,xdata);
Kdateval = kernel(xdata,xeval);
Kdiageval = kerneldiag(xeval);
if nargin >= 5
   Pdat = trend(xdata);
   Peval = trend(xeval);
   PTKinv = Pdat'/Kmat;
   PTKinvP = PTKinv*Pdat;
   Mmat = inv(Kmat) - PTKinv'*(PTKinvP\PTKinv);
end
