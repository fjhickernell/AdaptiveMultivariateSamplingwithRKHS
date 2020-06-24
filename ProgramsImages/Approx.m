function [Appx, fluctNorm, ErrBdx, ErrBd] = ...
   Approx(y, Kmat, Kdateval, errKXx, errKX, AX )
c = Kmat\y;
Appx = Kdateval'*c;
fluctNorm = sqrt(abs(y'*c));
ErrBdx = AX*errKXx*fluctNorm;
ErrBdx(find(isnan(ErrBdx)==1))=0;
ErrBd = AX*errKX*fluctNorm;

