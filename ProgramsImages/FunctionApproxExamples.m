%% Function approximation examples

%% Open File
clearvars
Outfile = 'FunctionApproxExOut.txt';
fopen(Outfile,'w+');
diary(Outfile)
StartTime = tic;

%% Defaults
S = struct('type','{}','subs',{{':'}});
abstolVec = [0.05 0.02 0.01 0.005 0.002 0.001]';
AlgName = {'Algo1', 'Algo2', 'Algo3'};
nAlg = length(AlgName);
[prm(1:nAlg).AlgName] = subsref(AlgName,S);
kernelth = @(t,x,theta) MaternKernel(t,x,theta,true);
[prm.kername] = subsref(repmat({'Matern'},1,nAlg),S);
[prm.n0] = subsref(repmat({1},1,nAlg),S);
prm(3).n0 = 5;

%% Simple function
f = @simpleFun;
[prm.fname] = subsref(repmat({'SimpleFun'},1,nAlg),S);
RunExample(f,prm,abstolVec,kernelth)

%% With uniform grid 
[prm.whDes] = subsref(repmat({'unif_grid'},1,nAlg),S);
RunExample(f,prm,abstolVec,kernelth)

%% Currin sinusoidal
f = @(x) sin(2*pi*(x-0.1));
[prm.fname] = subsref(repmat({'CurrinSineFun'},1,nAlg),S);
[prm.whDes] = subsref(repmat({'unif_grid'},1,nAlg),S);
[prm.yLim] = subsref(repmat({[-2;1.5]},1,nAlg),S);
RunExample(f,prm,abstolVec,kernelth)

%% With adpatpive theta 
[prm.whDes] = subsref(repmat({'adapt_th'},1,nAlg),S);
RunExample(f,prm,abstolVec,kernelth)

%% Left peak function
f = @(x) exp(-6*x).*sin(8*x+0.1) - 0.1;
[prm.fname] = subsref(repmat({'LeftPeakFun'},1,nAlg),S);
%abstolVec = [0.05 0.02 0.01]';
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

%% Finish up
toc(StartTime)
diary off


