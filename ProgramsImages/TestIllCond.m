%% Test ill conditoning

clearvars, close all

f = @(x) exp(-5*x).*cos(10*x);
th = 1;
n = 16;
cutoff = 1e-10;
ntest = 1e4;
xdata = (0:(n-1))'/(n-1);
xeval = (0:(ntest-1))'/(ntest-1);
fdata = f(xdata);
feval = f(xeval);
kernel = @(t,x) GaussKernel(t,x,th); %without the 4th arg, no transform
[Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata, xeval, kernel);
[errKXx, errKX] = powerfun(Kmat, Kdateval, Kdiageval, cutoff);
AX = ABfun(errKX,errKNull,1,0.2);
[Appx, AppxNorm, ErrBdx, ErrBd] = ...
   Approx(fdata, Kmat, Kdateval, errKXx, errKX, AX, cutoff);
trueErrX = feval - Appx;
trueErr = max(abs(trueErrX))
ErrBd

figure(1)
plot(xeval,feval, ...
   xeval,Appx, ...
   xdata,fdata,'.')

figure(2)
semilogy(xeval,abs(trueErrX),...
   xeval,ErrBdx)
   