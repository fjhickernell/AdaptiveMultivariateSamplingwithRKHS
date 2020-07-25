%% Example of simple function

gail.InitializeWorkspaceDisplay
format short e
warning('off')

[~,~,xeval,neval,Ainf,B0] = StdParam;
abstolVec = [0.05 0.02 0.01 0.005]';
ntol = size(abstolVec,1);

f = @(x) sin(2*pi*(x-0.1));
prm.fname = 'CurrinSineFun';
prm.colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];
prm.legendPos = 'south';
prm.isDiagnose = true;
prm.Ainf = Ainf;
prm.B0 = B0;
prm.n0 = 1;
prm.nmax = 500;
prm.theta = 1;
%prm.whDes = 'unif_grid';
prm.whDes = 'adapt_th';
prm.whobj = 'EmpBayesAx';
prm.plotSites = false;
prm.xLim = [0;1];
prm.yLim = [-2;2];
kernelth = @(t,x,theta) MaternKernel(t,x,theta,true);
kernel = @(t,x) MaternKernel(t,x,prm.theta);
prm.kername = 'Matern';
feval = f(xeval);


figure %simple function
plot(xeval,feval);
xlabel('\(x\)')
ylabel('\(f(x)\)');
axis([prm.xLim', prm.yLim'])
print('-depsc',[prm.fname 'Plot.eps'])

%% Algorithm 1 Sample size is adaptive
disp('Algorithm 1')
[Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
   AdaptAlgo1(f, kernel, xeval, feval, abstolVec, prm);
disp(['Necessary condition flag = ' int2str(NeccFlag(end))])
fprintf(1,'\n\n')

%% Algorithm 2 Sample location is adaptive
disp('Algorithm 2')
[Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
   AdaptAlgo2(f, kernel, xeval, feval, abstolVec, prm);
disp(['Necessary condition flag = ' int2str(NeccFlag(end))])
fprintf(1,'\n\n')

%% Algorithm 3 Sample location and kernel are adaptive
disp('Algorithm 3')
prm.n0 = 5;
prm.thetaRange = (-5:0.5:5)';
prm.whobj = 'EmpBayesAx';
[Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm] = ...
   AdaptAlgo3(f, kernelth, xeval, feval, abstolVec, prm);
fprintf(1,'\n\n')
   