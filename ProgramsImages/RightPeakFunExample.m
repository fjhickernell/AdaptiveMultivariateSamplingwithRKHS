%% Right peak function
S = struct('type','{}','subs',{{':'}}); prm = [];
[prm(1:3).AlgName] = subsref({'Algo2','Algo3','Algo3'},S);
[prm,kernelth] = parseFunAppxParam(prm);
nAlg = size(prm,2);
f = @(x) sin(pi*x.^4)-x;
[prm.fname] = subsref(repmat({'RightPeakFun'},1,nAlg),S);
[prm.kername] = subsref({'Matern','Matern','SpatialMatern'},S);
[prm.legendPos] = subsref(repmat({'southwest'},1,nAlg),S);
[prm.n0] = subsref({1,10,10},S);
[prm.theta] = subsref({1,1,[1 0]},S);
[prm.currentTheta] = subsref({0,0,[0 0]},S);
xRange = (-5:0.5:5)';
[thaa,thbb] = meshgrid(xRange,xRange);
[prm.thetaRange] = subsref({xRange,xRange,[thaa(:) thbb(:)]},S);
[prm.plotSites] = subsref({false,true,true},S);
[prm.yLim] = subsref(repmat({[-1.2;0.5]},1,nAlg),S);

%%
RunExample(f,prm,kernelth)

%% With uniform grid 
[prm.whDes] = subsref(repmat({'unif_grid'},1,nAlg),S);
RunExample(f,prm,kernelth)
