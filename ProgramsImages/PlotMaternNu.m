%% Plot Matern kernel
clearvars
nu = 40;
xplot = (0.001:0.001:3);
sc = sqrt(2*nu);
kernel = @ (x) (2^(1-nu)/gamma(nu)) .* (sc*x).^nu .* besselk(nu,sc*x);
figure
plot(xplot,kernel(xplot),xplot,exp(-xplot.^2/2),xplot,exp(-xplot))
