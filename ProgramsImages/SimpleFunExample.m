%% Example of easy function

gail.InitializeWorkspaceDisplay

abstol = 0.02;
theta = 1;
Ainf = 0.1;
B0 = 0.2;
xeval = (0:0.001:1)';
f = @simpleFun;
feval = f(xeval);
colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];
nmax = 500;
xdata(nmax,1) = 0;
fdata(nmax,1) = 0;
kerneldiag = @(x) ones(size(x,1),1);
errKNull = 1;


%% Algorithm 1 Sample size is adaptive
kernel = @(t,x) MaternKernel(t,x,theta);
plotn = [0 1 3 8 nmax];
ploti = 2;
figure
plot(xeval,feval,'color',colorScheme(1,:));
hold on
for n = 1:nmax
   xdata(n) = seqFixedDes(n);
   fdata(n) = f(xdata(n));
   [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:n,:), xeval, kernel, kerneldiag);
   [errKXx, errKX] = powerfun(Kmat, Kdateval, Kdiageval);
   AX = ABfun(errKX,errKNull,Ainf,B0);
   [Appx, fluctNorm, ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX );
   if n == plotn(ploti)
      plot(xeval,Appx,'color',colorScheme(mod(ploti-1,6)+1,:))
      nrange = plotn(ploti-1)+1:n;
      plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(ploti-1,6)+1,:))
      ploti = ploti+1;
   end
   if ErrBd < abstol
      plot(xeval,Appx,'color',colorScheme(mod(ploti-1,6)+1,:))
      nrange = plotn(ploti-1)+1:n;
      plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(ploti-1,6)+1,:))
      break
   end
end
trueErr = max(abs(feval - Appx));
disp('Algorithm 1')
disp(['True error = ' num2str(trueErr,3) ', Error bound = ' num2str(ErrBd,3) ...
   ', Using ' int2str(n) ' data'])
disp(' ')
xlabel('\(x\)')


%% Algorithm 2 Sample location is adaptive
plotn = [0 1 3 8 nmax];
ploti = 2;
figure
plot(xeval,feval,'color',colorScheme(1,:));
hold on
for n = 1:nmax
   if n == 1
      xdata(1) = seqFixedDes(1);
   else
      xdata(n) = xeval(whKX);
   end
   fdata(n) = f(xdata(n));
   [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:n,:), xeval, kernel, kerneldiag);
   [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   AX = ABfun(errKX,errKNull,Ainf,B0);
   [Appx, fluctNorm, ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX );
   if n == plotn(ploti)
      plot(xeval,Appx,'color',colorScheme(mod(ploti-1,6)+1,:))
      nrange = plotn(ploti-1)+1:n;
      plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(ploti-1,6)+1,:))
      ploti = ploti+1;
   end
   if ErrBd < abstol
      plot(xeval,Appx,'color',colorScheme(mod(ploti-1,6)+1,:))
      nrange = plotn(ploti-1)+1:n;
      plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(ploti-1,6)+1,:))
      break
   end
end
trueErr = max(abs(feval - Appx));
disp('Algorithm 2')
disp(['True error = ' num2str(trueErr,3) ', Error bound = ' num2str(ErrBd,3) ...
   ', Using ' int2str(n) ' data'])
disp(' ')
xlabel('\(x\)')


%% Algorithm 3 Sample location and kernel are adaptive
n0 = 4;
plotn = [0 n0:8 nmax];
ploti = 2;
figure
plot(xeval,feval,'color',colorScheme(1,:));
hold on
kernelth = @(t,x,lnth) MaternKernel(t,x,exp(lnth));
thetaRange = (-5:0.5:5)';
for n = n0:nmax
   if n == n0
      xdata(1:n0) = seqFixedDes(1:n0);
      fdata(1:n0) = f(xdata(1:n0));
   else
      xdata(n) = xeval(whKX);
      fdata(n) = f(xdata(n));
   end
   lnthOptim = selectTheta(thetaRange,kernelth,xdata(1:n),fdata(1:n));
   thetaOptim = exp(lnthOptim)
   kernel = @(t,x) MaternKernel(t,x,thetaOptim);
   [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:n,:), xeval, kernel, kerneldiag);
   [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   AX = ABfun(errKX,errKNull,Ainf,B0);
   [Appx, fluctNorm, ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX );
   if n == plotn(ploti)
      plot(xeval,Appx,'color',colorScheme(mod(ploti-1,6)+1,:))
      nrange = plotn(ploti-1)+1:n;
      plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(ploti-1,6)+1,:))
      ploti = ploti+1;
   end
   if ErrBd < abstol
      plot(xeval,Appx,'color',colorScheme(mod(ploti-1,6)+1,:))
      nrange = plotn(ploti-1)+1:n;
      plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(ploti-1,6)+1,:))
      break
   end
end
trueErr = max(abs(feval - Appx));
disp('Algorithm 3')
disp(['True error = ' num2str(trueErr,3) ', Error bound = ' num2str(ErrBd,3) ...
   ', Using ' int2str(n) ' data'])
disp(' ')
xlabel('\(x\)')
   
   