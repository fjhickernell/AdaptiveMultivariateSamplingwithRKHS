function [AX, BX] = ABfun(errKX,errKNull,Ainf,B0)
BX = errKX/errKNull;
if BX < B0
   AX = Ainf*B0/(B0 - BX);
else
   AX = inf;
end