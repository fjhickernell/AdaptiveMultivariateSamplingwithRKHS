function [Appx, fluctNorm, ErrBdx, ErrBd] = ...
   Approx(y, Kmat, Kdateval, errKXx, errKX, AX )
c = Kmat\y;
Appx = Kdateval'*c;
fluctNorm = sqrt(abs(y'*c));
ErrBdx = AX*errKXx*fluctNorm;
ErrBd = AX*errKX*fluctNorm;

