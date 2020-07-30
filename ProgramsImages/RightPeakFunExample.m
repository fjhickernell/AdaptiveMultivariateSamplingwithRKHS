%% Right peak function
f = @(x) sin(pi*x.^4)-x;
[prm.fname] = subsref(repmat({'RightPeakFun'},1,nAlg),S);
[prm.AlgName] = subsref({'Algo2','Algo3','Algo3'},S);
[prm.kername] = subsref({'Matern','Matern','SpatialMatern'},S);
[prm.legendPos] = subsref(repmat({'southwest'},1,nAlg),S);
[prm.n0] = subsref({1,5,5},S);
[prm.theta] = subsref({1,1,[1 0]},S);
xRange = (-5:0.5:5)';
[thaa,thbb] = meshgrid(xRange,xRange);
[prm.thetaRange] = subsref({xRange,xRange,[thaa(:) thbb(:)]},S);
[prm.plotSites] = subsref({false,true,true},S);
[prm.yLim] = subsref(repmat({[-1.2;0.5]},1,nAlg),S);
RunExample(f,prm,abstolVec,kernelth)
