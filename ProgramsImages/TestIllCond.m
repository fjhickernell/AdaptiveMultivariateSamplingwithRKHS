%% Test ill conditoning

clearvars, close all
InitializeDisplay


%% Rank of Gram matrix with as n increases
thetavec = [0.01 0.1 1 10 100];
nth = size(thetavec,2);
mmax = 6;
nmax = 2^mmax;
xdata_all = (0:nmax)'/nmax;
xeval = 0;
rankKmat(nth,mmax+1) = 0;
rankKmatB(nth,mmax+1) = 0;
legtxt = cell(nth,1);
for ii = 1:nth
   th = thetavec(ii);
   kernel = @(t,x) GaussKernel(t,x,th); %without the 4th arg, no transform
   legtxt(ii) = {['\(\theta = ' num2str(th) '\)']};
   for m = 0:mmax
      nskip =2^(mmax-m);
      xdata = xdata_all(1:nskip:nmax+1);
      Kmat = KMP(xdata, xeval, kernel);
      rankKmat(ii,m+1) = rank(Kmat);
      [~,rankKmatB(ii,m+1)] = KinvY(Kmat,zeros(size(Kmat,1),1));,
   end
end
figure;
nvec = 2.^(0:mmax)+1;
%h = loglog(nvec,rankKmatB');
legend(h,legtxt,'location','northwest','box','off')
xlabel('\(n\)')
ylabel('rank(\(\mathsf{K}\))')
axis([1 10^ceil(log10(nmax)) 1 10^ceil(log10(nmax))])
print('rankKernel.eps','-depsc')


%% Plot bases functions

thetavec = [0.01 0.1 1 10 100];
nth = size(thetavec,2);
nplot = 1e3;
xeval = (0:nplot)'/nplot;
n = 32;
y0 = zeros(n,1);
xdata = (0:n)'/nmax;
rankKmat(nth,mmax+1) = 0;
rankKmatB(nth,mmax+1) = 0;
legtxt = cell(nth,1);
for ii = 1:nth
   th = thetavec(ii);
   kernel = @(t,x) GaussKernel(t,x,th); %without the 4th arg, no transform
   legtxt(ii) = {['\(\theta = ' num2str(th) '\)']};
   Kmat = KMP(xdata, xeval, kernel);
   [~,rankKmatB(ii,m+1)] = KinvY(Kmat,y0);
end
figure;
nvec = 2.^(0:mmax)+1;
%h = loglog(nvec,rankKmat');
h = loglog(nvec,rankKmatB');
legend(h,legtxt,'location','northwest','box','off')
xlabel('\(n\)')
ylabel('rank(\(\mathsf{K}\))')
axis([1 10^ceil(log10(nmax)) 1 10^ceil(log10(nmax))])
print('rankKernel.eps','-depsc')

return
      

      
      
   

%% Test function
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
   