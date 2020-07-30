%% Adaptive Multivariate Sampling Example 4
gail.InitializeWorkspaceDisplay
format short e
warning('off')

[~,~,~,neval,Ainf,B0,errFudge] = StdParam;
%abstolVec = [0.05 0.02 0.01]';
abstolVec = 0.05;
ntol = size(abstolVec,1);
theta = 1;

f = @(x) 1/6*((30+5*x(:,1).*sin(5*x(:,1))).*...
    (4+exp(-5*x(:,2)))-100);
prm.fname = 'LimFun';
prm.colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];
prm.legendPos = 'northeast';
prm.isDiagnose = true;
prm.Ainf = Ainf;
prm.B0 = B0;
prm.n0 = 1;
prm.nmax = 20;
prm.theta = [1 1];
%prm.whDes = 'unif_grid';
prm.whDes = 'adapt_th';
prm.plotSites = false;
kernelth = @(t,x,theta) MaternKernel(t,x,theta,true);
kernel = @(t,x) MaternKernel(t,x,prm.theta);
prm.kername = 'Matern';
prm.isDiagnose = false;
x = 0:0.002:1;
[xx, yy]= meshgrid(x,x);
xeval  = [xx(:),yy(:)];
feval = f(xeval);
fplot = reshape(feval,length(x),length(x));

figure 
mesh(xx,yy,fplot);
xlabel('\(x_1\)')
ylabel('\(x_2\)');
zlabel('\(f(x_1,x_2)\)');
print('-depsc',[prm.fname 'Plot.eps'])



%% Algorithm 3 Sample location and kernel spatial dependence are adaptive
disp('Algorithm 3')
prm.kername = 'SpatialMatern';
prm.n0 = 10;
prm.theta = [1 1 0 0];
prm.whObj = 'EmpBayesAx';
[thaa,thbb] = meshgrid(-5:0.5:5,-5:0.5:5);
prm.thetaRange  = [thaa(:) thaa(:) thbb(:) thbb(:)];
[Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm] = ...
   AdaptAlgo3(f, kernelth, xeval, feval, abstolVec, prm);
fprintf(1,'\n\n')


