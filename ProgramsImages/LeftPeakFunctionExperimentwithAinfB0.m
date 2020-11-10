%% Left Peak Function Example
S = struct('type','{}','subs',{{':'}}); prm = [];
[prm.AlgName] = subsref({'Algo3'},S);
[prm,kernelth] = parseFunAppxParam(prm);
nAlg = size(prm,2);
f = @(x) exp(-6*x).*sin(8*x+0.1) - 0.1;
[prm.fname] = subsref(repmat({'LeftPeakFun'},1,nAlg),S);
[prm.kername] = subsref({'Matern'},S);
[prm.n0] = subsref({10},S);
[prm.theta] = subsref({1},S);
[prm.currentTheta] = subsref({1},S);
xRange = (-5:0.5:5)';
[thaa,thbb] = meshgrid(xRange,xRange);
[prm.thetaRange] = subsref({xRange,xRange,[thaa(:) thbb(:)]},S);
[prm.yLim] = subsref(repmat({[-0.2;0.5]},1,nAlg),S);
[prm.legendPos] = subsref(repmat({'northeast'},1,nAlg),S);
[prm.plotSites] = subsref({true},S);

%%
RunExample(f,prm,kernelth)
