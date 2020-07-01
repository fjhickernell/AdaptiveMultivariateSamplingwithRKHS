function [Kmat, Kdateval, Kdiageval, Pdat, Peval, PTKinv, PTKinvP, Mmat] = ...
   KMP(xdata, xeval, kernel, trend)
Kmat = kernel(xdata,xdata);
[Kdateval, Kdiageval] = kernel(xdata,xeval);
if nargin >= 4
   Pdat = trend(xdata);
   Peval = trend(xeval);
   PTKinv = Pdat'/Kmat;
   PTKinvP = PTKinv*Pdat;
   Mmat = inv(Kmat) - PTKinv'*(PTKinvP\PTKinv);
end
