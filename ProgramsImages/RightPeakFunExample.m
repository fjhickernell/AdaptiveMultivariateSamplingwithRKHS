%% Example for 
gail.InitializeWorkspaceDisplay
format short e
warning('off')

[~,~,xeval,neval,Ainf,B0,errFudge] = StdParam;
%abstolVec = [0.05 0.02 0.01 0.005 0.002 0.001]';
abstolVec = [0.05 0.02 0.01]';
ntol = size(abstolVec,1);
theta = 1;

f = @(x) sin(pi*x.^4)-x;
prm.fname = 'RightPeakFun';
feval = f(xeval);
prm.colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon;];
prm.legendPos = 'southwest';
prm.isDiagnose = true;
prm.Ainf = Ainf;
prm.B0 = B0;
prm.n0 = 1;
prm.nmax = 500;
prm.whDes = 'uniform';
prm.plotSites = true;

figure %Univariate function
plot(xeval,feval);
xlabel('\(x\)')
ylabel('\(f(x)\)');
print('-depsc',[prm.fname 'Plot.eps'])

disp('Algorithm 3')
prm.kername = 'SpatialMatern';
kernelth = @(t,x,th) MaternKernel(t,x,th,true);
diary(sprintf('run_%s_%d.txt',datestr(now,'yyyy_mm_dd_HH_MM_SS'),randi([1,10000],1)));
prm.n0 = 9;
%prm.whobj = 'EmpBayesAx';
prm.whobj = 'minErrBd';
prm.isDiagnose = true;
[thaa, thbb] = meshgrid(-5:0.5:5, -5:0.5:5);
prm.thetaRange = [thaa(:) thbb(:)];
[Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm] = ...
   AdaptAlgo3(f, kernelth,xeval, feval, abstolVec, prm);
fprintf(1,'\n\n')
