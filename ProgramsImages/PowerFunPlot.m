%% Plotting the ellipses for different values of theta
gail.InitializeWorkspaceDisplay

xplot = (0:0.001:1)';
nxplot = size(xplot,1);
theta = 1;

kernel = @(t,x) MaternKernel(t,x,theta);
kerneldiag = @(x) ones(size(x,1),1);
colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];
cOrder = [1:6];


%% First no trend
nvec = [4  16  64];
nn = length(nvec);
ndata = max(nvec);
xdata = seqFixedDes(1:ndata);

figpow = figure;
set(gca,'Yscale','log')
hold on
legendLabel = cell(length(nn));
for ii = 1:nn
   [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:nvec(ii)), xplot, kernel, kerneldiag);
   errKX = powerfun(Kmat, Kdateval, Kdiageval);
   plot(xplot,errKX,'color',colorScheme(cOrder(ii),:))
   legendLabel{ii} = ['\(n = ' int2str(nvec(ii)) ' \quad \)'];
end
plot(xdata,zeros(ndata,1),'.k')
xlabel('\(x\)')
ylabel('ERRK\((\mathsf{X},x)\)')
legend(legendLabel,'box','off','Orientation','horizontal')
print('-depsc','errKplot.eps')
   
   
%% Add a trend
polyfun = @(x)[ones(size(x,1),1) x-1/2];
nvec = [4 16 64];
nn = length(nvec);
ndata = max(nvec);
xdata = seqFixedDes(1:ndata);
figpowpolyK = figure;
set(gca,'Yscale','log')
hold on
figpowpolyP = figure; 
set(gca,'Yscale','log')
hold on
legendLabel = cell(length(nn));
for ii = 1:nn
   [Kmat, Kdateval, Kdiageval, ~, Peval, PTKinv] = ...
      KMP(xdata(1:nvec(ii)), xplot, kernel, kerneldiag, polyfun);
   [errKX,errPX] = powerfun(Kmat, Kdateval, Kdiageval, Peval, PTKinv);
   figure(figpowpolyK); plot(xplot,errKX,'color',colorScheme(cOrder(ii),:))
   figure(figpowpolyP); plot(xplot,errPX,'color',colorScheme(cOrder(ii),:))
   legendLabel{ii} = ['\(n = ' int2str(nvec(ii)) ' \quad\)'];
end
figure(figpowpolyK); 
plot(xdata,zeros(ndata,1),'.k')
xlabel('\(x\)')
ylabel('ERRK\((\mathsf{X},x)\)')
legend(legendLabel,'box','off','Orientation','horizontal')
figure(figpowpolyP); 
plot(xdata,zeros(ndata,1),'.k')
xlabel('\(x\)')
ylabel('ERRP\((\mathsf{X},x)\)')
legend(legendLabel,'box','off','Orientation','horizontal')
