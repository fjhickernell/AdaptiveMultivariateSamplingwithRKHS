%% Adaptive Multivariate Sampling Example 4
gail.InitializeWorkspaceDisplay
format short e
warning('off')

[~,~,xeval,neval,Ainf,B0,errFudge] = StdParam;
abstolVec = [0.05 0.02 0.01]';
ntol = size(abstolVec,1);
theta = 1;

f = @(x) exp(-6*x).*sin(8*x+0.1) - 0.1;
prm.fname = 'LeftPeakFun';
kernel = @(t,x) MaternKernel(t,x,theta);
prm.kername = 'Matern';
feval = f(xeval);
prm.colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];
prm.legendPos = 'northeast';
prm.isDiagnose = true;
prm.Ainf = Ainf;
prm.B0 = B0;
prm.n0 = 1;
prm.nmax = 500;
prm.whDes = 'uniform';

figure %Univariate function
plot(xeval,feval);
xlabel('\(x\)')
ylabel('\(f(x)\)');
print('-depsc',[prm.fname 'Plot.eps'])


%% Algorithm 3 Sample location and kernel are adaptive
disp('Algorithm 3')
kernelth = @(t,x,theta) MaternKernel(t,x,theta,true);
prm.n0 = 3;
prm.thetaRange = (-5:0.5:5)';
prm.whobj = 'EmpBayesAx';
prm.plotSites = true;
[Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm] = ...
   AdaptAlgo3(f, kernelth, xeval, feval, abstolVec, prm);
fprintf(1,'\n\n')


%% Algorithm 3 Sample location and kernel spatial dependence are adaptive
disp('Algorithm 3')
prm.kername = 'SpatialMatern';
kernelth = @(t,x,theta) MaternKernel(t,x,theta,true);
prm.n0 = 5;
[thaa,thbb] = meshgrid(-5:0.5:5,-5:0.5:1);
prm.thetaRange  = [thaa(:) thbb(:)];
[Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm] = ...
   AdaptAlgo3(f, kernelth, xeval, feval, abstolVec, prm);
fprintf(1,'\n\n')


