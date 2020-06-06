%% Plotting the ellipses for different values of theta
gail.InitializeWorkspaceDisplay

xplot = (0:0.001:1)';
nxplot = size(xplot,1);
theta = 1;

kernel = @(t,x) MaternKernel(t,x,theta);
kerneldiag = @(x) ones(size(x,1),1);

nvec = [4  16  64];
nn = length(nvec);
ndata = max(nvec);
xdata = seqFixedDes(1:ndata);

colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];
cOrder = [1:6];
figpow = figure;
set(gca,'Yscale','log')
hold on
legendLabel = cell(length(nn));
for ii = 1:nn
   powkeval = powerfun(xplot,xdata(1:nvec(ii)),kernel,kerneldiag);
   plot(xplot,powkeval,'color',colorScheme(cOrder(ii),:))
   legendLabel{ii} = ['\(n = ' int2str(nvec(ii)) ' \quad \)'];
end
plot(xdata,zeros(ndata,1),'.k')
xlabel('\(x\)')
ylabel('errK\((x)\)')
legend(legendLabel,'box','off','Orientation','horizontal')
print('-depsc','errKplot.eps')
   
   
   
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
   [powkeval,powpeval] = powerfun(xplot,xdata(1:nvec(ii)),kernel,kerneldiag,polyfun);
   figure(figpowpolyK); plot(xplot,powkeval,'color',colorScheme(cOrder(ii),:))
   figure(figpowpolyP); plot(xplot,powpeval,'color',colorScheme(cOrder(ii),:))
   legendLabel{ii} = ['\(n = ' int2str(nvec(ii)) ' \quad\)'];
end
figure(figpowpolyK); 
plot(xdata,zeros(ndata,1),'.k')
xlabel('\(x\)')
legend(legendLabel,'box','off','Orientation','horizontal')
figure(figpowpolyP); 
plot(xdata,zeros(ndata,1),'.k')
xlabel('\(x\)')
legend(legendLabel,'box','off','Orientation','horizontal')
