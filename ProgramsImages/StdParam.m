function [thetavec,nth,xplot,nxplot,Ainf,B0,errFudge] = ...
   StdParam(d)
%STDPARAM sets the standard parameters for our examples
thetavec = [1 4 16]';
nth = size(thetavec,1);
xplot = (0:0.0002:1)';
nxplot = size(xplot,1);
Ainf = 0.5;
B0 = 0.05;
errFudge = 1000*eps;
end

