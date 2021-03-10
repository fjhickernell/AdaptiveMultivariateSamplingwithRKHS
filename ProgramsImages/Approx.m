function [Appx, fluctNorm, ErrBdx, ErrBd] = ...
   Approx(y, Kmat, Kdateval, errKXx, errKX, AX, cutoff)
if nargin <= 6
   cutoff = 0;
end
c = KinvY(Kmat,y,cutoff);
Appx = Kdateval'*c;
fluctNorm = sqrt(abs(y'*c));
ErrBdx = AX*errKXx*fluctNorm;
ErrBdx(isnan(ErrBdx))=0;
ErrBd = AX*errKX*fluctNorm;
