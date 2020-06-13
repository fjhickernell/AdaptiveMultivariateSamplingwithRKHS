%% Plot of Power function vs theta and A(X) vs theta
gail.InitializeWorkspaceDisplay

[thetavec,nth,xplot,nxplot,Ainf,B0] = ...
   StdParam;

[~,kerneldiag,errKNull] = MaternKernel([],[],1);
colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];
cOrder = 1:6;

%% Plot of errK(X,.) and error bound for different theta
n = 16;
xdata = seqFixedDes(1:n);
ydata = simpleFun(xdata);
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
   [Kmat, Kdateval, Kdiageval] = KMP(xdata, xplot, kernel, kerneldiag);
   [errKXx,errKX] = powerfun(Kmat, Kdateval, Kdiageval);
   [AX, BX] = ABfun(errKX,errKNull,Ainf,B0);
   normAppx = sqrt(ydata'*(Kmat\ydata));
   errBdx = errKXx*(AX*normAppx);
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
n = 8;
xdata = seqFixedDes(1:n);
ydata = simpleFun(xdata);
figerrBdLots = figure;
set(gca,'Xscale','log')
set(gca,'Yscale','log')
axis([0.01 10 0.01 10])
hold on
lotsthvec = 10.^(-2:0.1:1)';
nthlots = length(lotsthvec);
errBd(nthlots,1) = 0;

for ii = 1:nthlots
   theta = lotsthvec(ii);
   kernel = @(t,x) MaternKernel(t,x,theta);
   [Kmat, Kdateval, Kdiageval] = KMP(xdata, xplot, kernel, kerneldiag);
   [errKXx,errKX] = powerfun(Kmat, Kdateval, Kdiageval);
   [AX, BX] = ABfun(errKX,errKNull,Ainf,B0);
   normAppx = sqrt(ydata'*(Kmat\ydata));
   errBd(ii) = errKX*(AX*normAppx);
end
figure(figerrBdLots);
plot(lotsthvec,errBd,'.');
xlabel('\(\theta\)')
ylabel('\(A(\mathsf{X})\)ERRK\((\mathsf{X}) \sqrt{\textit{\textbf{y}}^T\mathsf{K}^{-1}\textit{\textbf{y}}}\)')
print('-depsc','errBdplotthlots.eps')

%%
[~,whth] = min(errBd);


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
      [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:n,:), xplot, kernel, kerneldiag);
      [~,errKX,whKX] = powerfun(Kmat, Kdateval, Kdiageval);
      [AX(ii,n), BX(ii,n)] = ABfun(errKX,errKNull,Ainf,B0);
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
h(nth,1) = 0;

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
      h(ii) = plot(xdata(1:n),n*ones(n,1)+shift(ii),'.','color',colorScheme(cOrder(ii),:));
   end
   legendLabel{ii} = ['\(\theta = ' num2str(theta) ' \quad \)'];
end
ylabel('\(n\)')
xlabel('\(x\)')
legend(h,legendLabel,'box','off','Orientation','horizontal','location','north')
print('-depsc','Desplotth.eps')

    
