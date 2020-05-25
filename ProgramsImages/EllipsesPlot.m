%% Plotting the ellipses for different values of theta
gail.InitializeWorkspaceDisplay

f = @(x) exp(-6.*(x-1/2).^2) - 0.6;
xplot = (0:0.001:1)';
nxplot = size(xplot,1);
fplot = f(xplot);
plot(xplot,fplot)
xlabel('\(x\)')
ylabel('\(f(x)\)')
print -depsc samplef1D.eps

phiplot = xplot*(2*pi);
spherevec = [cos(phiplot) sin(phiplot)];

xdata = [0.3; 0.6];
ydata = f(xdata);

efigh = figure;
hold on
kfigh = figure;
hold on
colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple];
cOrder = [1 3 2 4];
thetavec = [0.5 8 2]';
thOrder = [1 3 2];
ntheta = size(thetavec,1);
h=zeros(ntheta,1);
hK=zeros(ntheta,1);
legendLabel = cell(length(thetavec));
for i = 1:ntheta
   theta = thetavec(i);
   Kmat = MaternKernel(xdata,xdata,theta);
   L = chol(Kmat,'lower');
   ellipsevec = sqrt(ydata'*(Kmat\ydata))*spherevec*L';
   figure(efigh);
   h(i)= fill(ellipsevec(:,1), ellipsevec(:,2), colorScheme(cOrder(i),:), ...
      'EdgeColor', colorScheme(cOrder(i),:));
   legendLabel{i} = ['\(\theta = ' num2str(theta) '\)'];
   Kplot = MaternKernel(xplot,0.5,theta);
   figure(kfigh);
   hK(i)=plot(xplot,Kplot,'-','color',colorScheme(cOrder(i),:));
end
figure(efigh);
plot(ydata(1),ydata(2),'k.');
xlabel('\(z_1\)')
ylabel('\(z_2\)')
legend(h(thOrder),legendLabel{thOrder},  ...
   'location','northwest','box','off');
print -depsc ellipsesPlot.eps

figure(kfigh);
xlabel('\(x\)')
ylabel('\(K_{\theta}(x,0.5)\)')
legend(hK(thOrder),legendLabel{thOrder},  ...
   'location','south','box','off');
print -depsc KthetaPlot.eps




% 