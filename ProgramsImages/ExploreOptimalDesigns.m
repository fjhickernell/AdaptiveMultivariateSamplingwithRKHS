%%Explore designs
clearvars
InitializeDisplay
xeval = (0:0.0001:1)';
prm.nmax = 17;
prm.isDiagnose = true;
prm.xLim = [0 1]';
prm.legendPos = 'north';
prm.colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];

prm.kername = 'MaternOne'
[xdata, prm] = OptimalDesignXonly(@MaternKernelOne,xeval, prm);

[xdata, prm] = OptimalDesignXonly(@GaussKernel,xeval, prm);
