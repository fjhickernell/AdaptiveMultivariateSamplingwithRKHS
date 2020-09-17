S = struct('type','{}','subs',{{':'}}); prm=[];
[prm(1).AlgName] = subsref({'Algo3'},S);
[prm.kername] = subsref({'SpatialMatern'},S);
[prm,kernelth] = parseFunAppxParam(prm);
nAlg = size(prm,2);
%f = @(x) 1/6*((30+5*x(:,1).*sin(5*x(:,1))).*...
%    (4+exp(-5*x(:,2)))-100);
f = @(x) cos(sum(x,2)).*exp(prod(x,2));
[prm.fname] = subsref(repmat({'ChengDuFun'},1,nAlg),S);
%[prm.kername] = subsref({'SpatialMatern'},S);
[prm.n0] = subsref({10},S);
[prm.theta] = subsref({[1 1 0 0]},S);
xRange = (-5:0.5:5)';
[thaa,thbb] = meshgrid(xRange,xRange);
[prm.thetaRange] = subsref({...
    [thaa(:) thaa(:) thbb(:) thbb(:)]},S);
[prm.yLim] = subsref(repmat({[-0.2;0.5]},1,nAlg),S);
[prm.legendPos] = subsref(repmat({'northeast'},1,nAlg),S);
[prm.plotSites] = subsref({false},S);
[prm.abstolVec] = subsref({[0.1 0.05 0.02 0.01]'},S);
[prm.nmax] = subsref({500},S);
[prm.canvasTheta]= subsref({false},S);
[prm.currentTheta]= subsref({[0 0 0 0]},S);
[prm.whDes] = subsref({'adapt_th'},S);
[prm.whObj] = subsref({'EmpBayes'},S);
highdimflag = 2;
%%
RunExample(f,prm,kernelth,highdimflag)
