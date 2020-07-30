%% Left Peak Function Example
f = @(x) exp(-6*x).*sin(8*x+0.1) - 0.1;
[prm.fname] = subsref(repmat({'LeftPeakFun'},1,nAlg),S);
[prm.AlgName] = subsref({'Algo2','Algo3','Algo3'},S);
[prm.kername] = subsref({'Matern','Matern','SpatialMatern'},S);
[prm.n0] = subsref({1,5,5},S);
[prm.theta] = subsref({1,1,[1 0]},S);
xRange = (-5:0.5:5)';
[thaa,thbb] = meshgrid(xRange,xRange);
[prm.thetaRange] = subsref({xRange,xRange,[thaa(:) thbb(:)]},S);
[prm.yLim] = subsref(repmat({[-0.2;0.5]},1,nAlg),S);
[prm.legendPos] = subsref(repmat({'northeast'},1,nAlg),S);
[prm.plotSites] = subsref({false,true,true},S);
RunExample(f,prm,abstolVec,kernelth)
