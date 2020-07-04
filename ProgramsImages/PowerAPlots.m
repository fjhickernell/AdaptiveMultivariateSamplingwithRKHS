%% Plot of Power function vs theta and A(X) vs theta
gail.InitializeWorkspaceDisplay

warning('off','MATLAB:handle_graphics:exceptions:SceneNode');

[thetavec,nth,xplot,nxplot,Ainf,B0] = ...
   StdParam;

%[~,kerneldiag] = MaternKernel([],[],1);
colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];
cOrder = 1:6;

n = 4;

myFun = @(x) 2 * (exp(-60.*(x-1/2).^2) - 0.6);
xdata = seqFixedDes(1:n);
ydata = myFun(xdata);
yplot = myFun(xplot);


%% Plot of errK(X,.) AND error bound for different theta
figerrK = figure;
set(gca,'Yscale','log')
axis([0 1 1e-5 10])
hold on
figerrBd = figure;
set(gca,'Yscale','log')
hold on
legendLabel = cell(length(nth));

for ii = 1:nth
   theta = thetavec(ii);
   kernel = @(t,x) MaternKernel(t,x,theta);
   [Kmat, Kdateval, Kdiageval] = KMP(xdata, xplot, kernel);
   [errKXx,errKX] = powerfun(Kmat, Kdateval, Kdiageval);
   AX = ABfun(errKX,max(Kdiageval),Ainf,B0);
   [~,~,errBdx] = Approx(ydata, Kmat, Kdateval, errKXx, errKX, AX);
   figure(figerrK)
   plot(xplot,errKXx,'color',colorScheme(cOrder(ii),:))
   figure(figerrBd);
   plot(xplot,errKXx,'color',colorScheme(cOrder(ii),:));
   legendLabel{ii} = ['\(\theta = ' num2str(theta) ' \quad \)'];
end
figure(figerrK)
plot(xdata,zeros(n,1),'.k')
xlabel('\(x\)')
ylabel('ERRK\((\mathsf{X},x)\)')
legend(legendLabel,'box','off','Orientation','horizontal','location','north')
print('-depsc','errKplotth.eps')
figure(figerrBd);
plot(xdata,zeros(n,1),'.k')
xlabel('\(x\)')
ylabel('\(A(\mathsf{X})\)ERRK\((\mathsf{X},x) \sqrt{\textit{\textbf{y}}^T\mathsf{K}^{-1}\textit{\textbf{y}}}\)')
legend(legendLabel,'box','off','Orientation','horizontal','location','north')
print('-depsc','errBdplotth.eps')

%% Plot of error bound for many different theta
figerrBdLots = figure;
set(gca,'Xscale','log')
set(gca,'Yscale','log')
axis([1e-3 10 0.01 10])
hold on
lotsthvec = 10.^(-3:0.1:1)';
nthlots = length(lotsthvec);
errBd(nthlots,1) = 0;

for ii = 1:nthlots
   theta = lotsthvec(ii);
   kernel = @(t,x) MaternKernel(t,x,theta);
   [Kmat, Kdateval, Kdiageval] = KMP(xdata, xplot, kernel);
   [errKXx,errKX] = powerfun(Kmat, Kdateval, Kdiageval);
   AX = ABfun(errKX,max(Kdiageval),Ainf,B0);
   [~,~,~,errBd(ii)] = Approx(ydata, Kmat, Kdateval, errKXx, errKX, AX);
end
figure(figerrBdLots);
plot(lotsthvec,errBd,'.');
xlabel('\(\theta\)')
ylabel('\(A(\mathsf{X})\)ERRK\((\mathsf{X}) \sqrt{\textit{\textbf{y}}^T\mathsf{K}^{-1}\textit{\textbf{y}}}\)')
print('-depsc','errBdplotthlots.eps')

%% Plot function, approximation and error bars for theta with minimun error bound
kernelth = @(t,x,lnth) MaternKernel(t,x,exp(lnth));
thetaRange = (-5:0.5:5)';
lnthOptim = selectTheta(thetaRange,kernelth,xdata,ydata, ...
   xplot,Ainf,B0,'minErrBd');
theta = exp(lnthOptim)
kernel = @(t,x) MaternKernel(t,x,theta);
[Kmat, Kdateval, Kdiageval] = KMP(xdata, xplot, kernel);
[errKXx,errKX] = powerfun(Kmat, Kdateval, Kdiageval);
AX = ABfun(errKX,max(Kdiageval),Ainf,B0);
[Appx, fluctNorm, ErrBdx, ErrBd] = ...
   Approx(ydata, Kmat, Kdateval, errKXx, errKX, AX);
figure
hold on
h = plot(xplot,yplot,xdata,ydata,'.');
h = [h; plot(xplot,Appx,'color', MATLABGreen)];
h = [h; plot(xplot,Appx+ErrBdx*[-1,1],'color', MATLABPurple)];
legend(h(1:3), {'\(f(x)\)','APP\((\mathsf{X},\textit{\textbf{y}})\)', ...
   'APP\((\mathsf{X},\textit{\textbf{y}}) \pm \)ERRBD\((\mathsf{X},\textit{\textbf{y}},x)\)'}, ...
   'location','south','orientation','vertical','box','off')
trueErr = max(abs(yplot-Appx))
ErrBd
print('-depsc','minErrBdThAppx.eps')

%% Plot function, approximation and error bars for theta with our selection criterion
kernelth = @(t,x,lnth) MaternKernel(t,x,exp(lnth));
thetaRange = (-5:0.5:5)';
lnthOptim = selectTheta(thetaRange,kernelth,xdata,ydata, ...
   xplot,Ainf,B0,'EmpBayesAx');
theta = exp(lnthOptim)
kernel = @(t,x) MaternKernel(t,x,theta);
[Kmat, Kdateval, Kdiageval] = KMP(xdata, xplot, kernel);
[errKXx,errKX] = powerfun(Kmat, Kdateval, Kdiageval);
AX = ABfun(errKX,max(Kdiageval),Ainf,B0);
[Appx, fluctNorm, ErrBdx, ErrBd] = ...
   Approx(ydata, Kmat, Kdateval, errKXx, errKX, AX);
figure
hold on
h = plot(xplot,yplot,xdata,ydata,'.');
h = [h; plot(xplot,Appx,'color', MATLABGreen)];
h = [h; plot(xplot,Appx+ErrBdx*[-1,1],'color', MATLABPurple)];
legend(h(1:3), {'\(f(x)\)','APP\((\mathsf{X},\textit{\textbf{y}})\)', ...
   'APP\((\mathsf{X},\textit{\textbf{y}}) \pm \)ERRBD\((\mathsf{X},\textit{\textbf{y}},x)\)'}, ...
   'location','south','orientation','vertical','box','off')
trueErr = max(abs(yplot-Appx))
ErrBd
print('-depsc','OurSelectThAppx.eps')


%% Plot of A(X) for different theta and n
nmax = 100;
xdata = seqFixedDes(1:nmax);
figure
set(gca,'Xscale','log')
axis([1 nmax 0 0.8])
hold on
legendLabel = cell(length(nth));
AX(nth,nmax) = 0;
BX(nth,nmax) = 0;

for ii = 1:nth
   theta = thetavec(ii);
   kernel = @(t,x) MaternKernel(t,x,theta);
   for n = 1:nmax   
      if n == 1
         xdata(1) = seqFixedDes(1);
      else
         xdata(n) = xplot(whKX);
      end
      [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:n,:), xplot, kernel);
      [~,errKX,whKX] = powerfun(Kmat, Kdateval, Kdiageval);
      [AX(ii,n), BX(ii,n)] = ABfun(errKX,max(Kdiageval),Ainf,B0);
   end
   plot(1:nmax,AX(ii,:),'.','color',colorScheme(cOrder(ii),:))
   legendLabel{ii} = ['\(\theta = ' num2str(theta) ' \quad \)'];
end
xlabel('\(n\)')
ylabel('A\((\mathsf{X})\)')
legend(legendLabel,'box','off','Orientation','horizontal','location','north')
print('-depsc','Aplotth.eps')

%% Plot of designs for different theta
figure
nmax = 13;
axis([0 1 0.5 nmax+3])
hold on
legendLabel = cell(length(nth));
shift = [-0.3 0 0.3];
h = zeros(nth,1);

for ii = 1:nth
   theta = thetavec(ii);
   kernel = @(t,x) MaternKernel(t,x,theta);
   for n = 1:nmax   
      if n == 1
         xdata(1) = seqFixedDes(1);
      else
         xdata(n) = xplot(whKX);
      end
      [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:n,:), xplot, kernel);
      [~,errKX,whKX] = powerfun(Kmat, Kdateval, Kdiageval);
      h(ii) = plot(xdata(1:n),n*ones(n,1)+shift(ii),'.','color',colorScheme(cOrder(ii),:));
   end
   legendLabel{ii} = ['\(\theta = ' num2str(theta) ' \quad \)'];
end
ylabel('\(n\)')
xlabel('\(x\)')
legend(h,legendLabel,'box','off','Orientation','horizontal','location','north')
print('-depsc','Desplotth.eps')

%% Plot of orthonormal basis for different theta
thetaVec = [0.25 1 4 16]';
nth = size(thetaVec,1);

n = 3;
xdata = seqFixedDes(1:n);

figure
axis([0 1 -1 1.5])
hold on
legendLabel = cell(length(nth));
h = zeros(nth,1);

for ii = 1:nth
   theta = thetaVec(ii);
   kernel = @(t,x) MaternKernel(t,x,theta);
   [Kmat, Kdateval] = KMP(xdata, xplot, kernel);
   [eigvec,eigval] = eig(Kmat);
   Fbasis = (((diag(eigval).^(-1/2)).*eigvec')*Kdateval)';
   htemp = plot(xplot,Fbasis,'color',colorScheme(cOrder(ii),:));
   h(ii) = htemp(1);
   legendLabel{ii} = ['\(\theta = ' num2str(theta) ' \quad \)'];
end
plot([0 1],[0 0],'color','k','linewidth',1)
plot(xdata,zeros(n,1),'.k')
xlabel('\(x\)')
ylabel('Bases')
legend(h,legendLabel,'box','off','Orientation','vertical','location','north','NumColumns',2)
print('-depsc','Basisplotth.eps')
