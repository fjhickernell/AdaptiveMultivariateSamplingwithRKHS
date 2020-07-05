%% Example of easy function

gail.InitializeWorkspaceDisplay
format short e
warning('off')

[~,~,xeval,neval,Ainf,B0] = StdParam;
abstolVec = [0.05 0.02 0.01 0.005 0.002 0.001]';
ntol = size(abstolVec,1);
theta = 1;

f = @simpleFun;
fname = 'SimpleFun';
kernel = @(t,x) MaternKernel(t,x,theta);
kername = 'Matern';
feval = f(xeval);
colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];
nmax = 500;
xdata(nmax,1) = 0;
fdata(nmax,1) = 0;
errKNull = 1;

figure %simple function
plot(xeval,feval);
xlabel('\(x\)')
ylabel('\(f(x)\)');
print('-depsc','SimpleFunPlot.eps')

%% Algorithm 1 Sample size is adaptive
disp('Algorithm 1')
[Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
   AdaptAlgo1(f, kernel, xdata, fdata, xeval, feval, ...
   abstolVec, Ainf, B0, true, colorScheme, fname, kername);
disp(['Necessary condition flag = ' int2str(NeccFlag(end))])

%% Algorithm 2 Sample location is adaptive
disp('Algorithm 2')
[Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
   AdaptAlgo2(f, kernel, xdata, fdata, xeval, feval, ...
   abstolVec, Ainf, B0, true, colorScheme, fname, kername);
disp(['Necessary condition flag = ' int2str(NeccFlag(end))])

%% Algorithm 3 Sample location and kernel are adaptive
disp('Algorithm 3')
kernelth = @(t,x,theta) MaternKernel(t,x,theta,true);
n0 = 3;
thetaRange = (-5:0.5:5)';
[Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm] = ...
   AdaptAlgo3(f, kernelth, xdata, fdata, xeval, feval, ...
   abstolVec, Ainf, B0, n0, thetaRange, ...
   true, colorScheme, fname, kername);

   