%% Plot of Power function vs theta and A(X) vs theta
gail.InitializeWorkspaceDisplay

xplot = (0:0.001:1)';
nxplot = size(xplot,1);

thetavec = [1 4 16]';
nth = size(thetavec,1);
theta = 1;

kerneldiag = @(x) ones(size(x,1),1);
colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];
cOrder = 1:6;

%% Plot of errK(X,.) for different theta
n = 8;
xdata = seqFixedDes(1:n);
figpow = figure;
set(gca,'Yscale','log')
axis([0 1 1e-5 10])
hold on
legendLabel = cell(length(nth));

for ii = 1:nth
   theta = thetavec(ii);
   kernel = @(t,x) MaternKernel(t,x,theta);
   [Kmat, Kdateval, Kdiageval] = KMP(xdata, xplot, kernel, kerneldiag);
   errKXx = powerfun(Kmat, Kdateval, Kdiageval);
   plot(xplot,errKXx,'color',colorScheme(cOrder(ii),:))
   legendLabel{ii} = ['\(\theta = ' num2str(theta) ' \quad \)'];
end
plot(xdata,zeros(n,1),'.k')
xlabel('\(x\)')
ylabel('ERRK\((\mathsf{X},x)\)')
legend(legendLabel,'box','off','Orientation','horizontal')
print('-depsc','errKplotth.eps')

%% Plot of A(X) for different theta and n
Ainf = 0.2;
B0 = 0.5;
errKNull = 1;
nmax = 100;
xdata = seqFixedDes(1:nmax);
figpow = figure;
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
      [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:n,:), xplot, kernel, kerneldiag);
      [~,errKX,whKX] = powerfun(Kmat, Kdateval, Kdiageval);
      [AX(ii,n), BX(ii,n)] = ABfun(errKX,errKNull,Ainf,B0);
   end
   plot(1:nmax,AX(ii,:),'.','color',colorScheme(cOrder(ii),:))
   legendLabel{ii} = ['\(\theta = ' num2str(theta) ' \quad \)'];
end
xlabel('\(n\)')
ylabel('A\((\mathsf{X})\)')
legend(legendLabel,'box','off','Orientation','horizontal')
print('-depsc','Aplotth.eps')

   
