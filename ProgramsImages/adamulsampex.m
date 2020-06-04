%% Adaptive Multivariate Sampling Examples
%% 1-D example
f = @(x) exp(-6*x).*sin(8*x+0.1) - 0.1;
%f = @(x) f0(x);
axisBox = [0 1 -0.3 0.5];

%% Plot function
xData = [0:0.1:0.6 0.8:0.1:1]';
fData = f(xData);
xPlot = (0:0.002:1)';
%1D Common Matern Kernel
dist = @(x,y) abs(x - y');
kernel =@(dis) (1 + dis).*exp(-dis);

% To derive B0
KDataData = kernel(dist(xData,xData));
KPlotData = kernel(dist(xPlot,xData));
KPlotPlot = kernel(dist(xPlot,xPlot));
bot = norm(KPlotPlot,inf);
top = norm(KPlotPlot -KPlotData*(KDataData\KPlotData'),inf);
BData = sqrt(top/bot); %0.0063

B0 = 0.007;
Ainf = 10;
AData = Ainf*B0/(B0-BData);
