%%Explore designs
clearvars
InitializeDisplay
xeval = (0:0.0001:1)';
prm.nmax = 17;
prm.plotn = [1 2 4 8 16];
prm.isDiagnose = true;
prm.xLim = [0 1]';
prm.yLim = [1e-2 1e1]';
prm.legendPos = 'north';
prm.colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];
prm.kername = 'MaternOne';
[xdataOne, prm] = OptimalDesignXonly(@MaternKernelOne,xeval, prm);
Even = (0:prm.nmax-1)'/(prm.nmax-1);
compare = [sort(xdataOne) Even abs(sort(xdataOne) - Even)]


prm.nmax = 9;
prm.plotn = [1 2 4 8];
prm.yLim = [1e-7 1e2]';
prm.kername = 'GaussianKernel';
theta = 1;
[xdataGauss, prm] = OptimalDesignXonly(@(t,x) GaussKernel(t,x,theta),xeval, prm);
xdataGauss = sort(xdataGauss);
Even = (0:prm.nmax-1)'/(prm.nmax-1);
Cheby = (1+sin(pi*(-1/2 + Even)))/2
compare = [xdataGauss Cheby abs(xdataGauss - Cheby)]

figure
plot(Even,xdataGauss,'.','color',MATLABBlue)
hold on
plot(Even,Cheby,'.','color',MATLABOrange)