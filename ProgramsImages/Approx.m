function [Appx, fluctNorm, ErrBdx, ErrBd] = ...
   Approx(y, Kmat, Kdateval, errKXx, errKX, AX )
c = Kmat\y;
Appx = Kdateval'*c;
fluctNorm = sqrt(y'*c);
ErrBdx = AX*errKXx*fluctNorm;
ErrBd = AX*errKX*fluctNorm;

