function [Appx, fluctNorm, ErrBdx, ErrBd] = ...
   Approx(y, Kmat, Kdateval, errKXx, errKX, AX )
%c = Kmat\y;
c = KinvY(Kmat,y);
Appx = Kdateval'*c;
fluctNorm = sqrt(abs(y'*c));
ErrBdx = AX*errKXx*fluctNorm;
ErrBdx(isnan(ErrBdx))=0;
ErrBd = AX*errKX*fluctNorm;
