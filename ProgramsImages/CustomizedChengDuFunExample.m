S = struct('type','{}','subs',{{':'}}); prm=[];
[prm(1).AlgName] = subsref({'Algo3'},S);
%[prm.kername] = subsref({'SpatialMatern'},S);
[prm.kername] = subsref({'Matern'},S);
[prm,kernelth] = parseFunAppxParam(prm);
nAlg = size(prm,2);
%alpha = [2*pi 0.0001]';
alpha = [0.0001 2*pi ]';
f = @(x) cos(x*alpha);
[prm.fname] = subsref(repmat({'CustChengDuFun'},1,nAlg),S);
%[prm.kername] = subsref({'SpatialMatern'},S);
[prm.n0] = subsref({10},S);
%[prm.theta] = subsref({[1 1 0 0]},S);
[prm.theta] = subsref({[1 1]},S);
xRange = (-5:0.5:5)';
[thaa,thbb] = meshgrid(xRange,xRange);
%[prm.thetaRange] = subsref({...
%   [thaa(:) thaa(:) thbb(:) thbb(:)]},S);
[prm.thetaRange] = subsref({[thaa(:) thaa(:)]},S);
[prm.yLim] = subsref(repmat({[-0.2;0.5]},1,nAlg),S);
[prm.legendPos] = subsref(repmat({''},1,nAlg),S);
[prm.plotSites] = subsref({false},S);
%[prm.abstolVec] = subsref({[0.1 0.05 0.02 ]'},S);
[prm.abstolVec] = subsref({[0.1 0.05]'},S);
[prm.nmax] = subsref({500},S);
[prm.canvasTheta]= subsref({false},S);
[prm.currentTheta]= subsref({[0 0]},S);
%[prm.currentTheta]= subsref({[0 0 0 0]},S);
[prm.whDes] = subsref({'adapt_th'},S);
[prm.whObj] = subsref({'EmpBayesAx'},S);
highdimflag = 2;
%%
tic
RunExample(f,prm,kernelth,highdimflag)
toc
