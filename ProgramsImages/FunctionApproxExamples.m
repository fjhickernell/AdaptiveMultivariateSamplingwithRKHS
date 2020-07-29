%% Function approximation examples

%% Open File
clearvars
Outfile = 'FunctionApproxExOut.txt';
fopen(Outfile,'w+');
diary(Outfile)

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
RunExample(f,prm,abstolVec,kernelth,AlgName)

%% With uniform grid 
[prm.whDes] = subsref(repmat({'unif_grid'},1,nAlg),S);
RunExample(f,prm,abstolVec,kernelth,AlgName)

%% Currin sinusoidal
f = @(x) sin(2*pi*(x-0.1));
[prm.fname] = subsref(repmat({'CurrinSineFun'},1,nAlg),S);
[prm.yLim] = subsref(repmat({[-2;1.5]},1,nAlg),S);
RunExample(f,prm,abstolVec,kernelth,AlgName)

%% With adpatpive theta 
[prm.whDes] = subsref(repmat({'adapt_th'},1,nAlg),S);
RunExample(f,prm,abstolVec,kernelth,AlgName)

diary off


