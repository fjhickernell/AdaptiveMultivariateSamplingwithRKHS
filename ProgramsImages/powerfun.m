function [errKXx,errKX,whKX,errPXx,errPX] = ...
   powerfun(Kmat, Kdateval, Kdiageval, cutoff, Peval, PTKinv)
%errKXx = sqrt(abs(Kdiageval - sum(Kdateval .* (Kmat\Kdateval),1)'));
if nargin <=3
   cutoff = 0;
end
errKXx = sqrt(abs(Kdiageval - sum(Kdateval .* (KinvY(Kmat,Kdateval,cutoff)),1)'));
[errKX,whKX] = max(errKXx);
if nargin >= 6
   PTKinvKdateval = PTKinv*Kdateval;
   errPXx = sqrt(abs(sum((Peval-PTKinvKdateval').^2,2)));
   errPX = max(errPXx);
end

   
