%% Adaptive Multivariate Sampling Example 4
gail.InitializeWorkspaceDisplay
format short e
warning('off')

[~,~,xeval,neval,Ainf,B0,errFudge] = StdParam;
abstolVec = [0.05 0.02 0.01 0.005 0.002 0.001]';
ntol = size(abstolVec,1);
theta = 1;

f = @(x) exp(-6*x).*sin(8*x+0.1) - 0.1;
feval = f(xeval);
colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon;];
nmax = 500;
xdata(nmax,1) = 0;
fdata(nmax,1) = 0;
kerneldiag = @(x) ones(size(x,1),1);
errKNull = 1;

figure %Univariate function
plot(xeval,feval);
xlabel('\(x\)')
ylabel('\(f(x)\)');
print('-depsc','UniFunPlot.eps')

%% Algorithm 3 Sample location and kernel are adaptive
n0 = 3;
plotn = [0 n0 nmax];
ploti = 2;
h(ntol+length(plotn)+1,1) = 0; 
legendLabel = cell(ntol+length(plotn)+1,1);
figure
h(1) = plot(xeval,feval,'color',colorScheme(1,:));
legendLabel{1} = '\(f(x)\)';
hold on
itol = 1;
abstol = abstolVec(itol);
nNeed(ntol,1) = 0;
ErrBdVec(nmax,1) = 0;
trueErr(nmax,1) = 0;
InErrBars(nmax,1) = 0;
nstart = 0;
coli = ploti;
AXvec(nmax,1) = 0;
BXvec(nmax,1) = 0;
thOptimVec(nmax,1) = 0;
%kernelth = @(t,x,lnth) MaternKernel(t,x,exp(lnth));
%thetaRange = (-5:0.5:5)';
kernelth = @(t,x,lnth) MaternKernel(t,x,lnth);
thetaRange = (0.1:0.1:2)';
for n = n0:10
   if n == n0
      xdata(1:n0) = seqFixedDes(1:n0);
      fdata(1:n0) = f(xdata(1:n0));
   else
      xdata(n) = xeval(whKX);
      fdata(n) = f(xdata(n));
   end
   lnthOptim = selectTheta(thetaRange,kernelth,kerneldiag,xdata(1:n),fdata(1:n), ...
      xeval,errKNull,Ainf,B0);
   thetaOptim = exp(lnthOptim);
   thOptimVec(n,:) = thetaOptim;
   kernel = @(t,x) MaternKernel(t,x,thetaOptim);
   [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:n,:), xeval, kernel, kerneldiag);
   [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   [AX, BX] = ABfun(errKX,errKNull,Ainf,B0);
   AXvec(n) = AX;
   BXvec(n) = BX;   
   [Appx, fluctNorm, ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX );
   ErrBdVec(n) = ErrBd;
   trueErr(n) = max(abs(feval - Appx));
   errFudge = eps*sqrt(cond(Kmat));
   InErrBars(n) = sum(abs(feval - Appx) <= ErrBdx + errFudge)/neval;
   %if InErrBars(n) < 1, keyboard, end
   %if n == plotn(ploti)
      h(coli) = plot(xeval,Appx,'color',colorScheme(mod(coli,6)+1,:));
      legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
      %plot(xeval,Appx + trueErr(n),'-.','color',colorScheme(mod(coli,6)+1,:));
      %plot(xeval,Appx - trueErr(n),'-.','color',colorScheme(mod(coli,6)+1,:));
      nrange = nstart+1:n;
      plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(coli,6)+1,:))
      ploti = ploti+1;
      coli = coli+1;
      nstart = n;
   %end
   %if ErrBd < abstol
    %  if abstol >= 0.01
    %    h(coli) = plot(xeval,Appx,'color',colorScheme(mod(coli,6)+1,:));
    %     legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
    %     nrange = nstart+1:n;
    %     plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(coli,6)+1,:))
    %     nstart = n;
    %     coli = coli+1;
    %  end
    %  nNeed(itol) = n;
    %  itol = itol + 1;
    %  if itol > ntol, break, end
    %  abstol = abstolVec(itol);
   %end
end
disp('Algorithm 3 with Univariate Function')
xlabel('\(x\)')
ylabel('\(f(x), \ f(x_i), \ \)APP\((\mathsf{X},\textbf{\textit{y}})\)')
legend(h(1:coli-1),legendLabel{1:coli-1},'location','north','orientation','vertical','box','off')
print('-depsc','UniFunAlg3.eps')
Alg3SummaryData = [(n0:n)' ErrBdVec(n0:n) trueErr(n0:n) InErrBars(n0:n)];
hold off;
figure
plot(xdata(1:n),'.')
title('Sample Points');
print('-depsc','UniFunAlg3SamplePoints.eps')