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
SimpleFunExample

%% Currin sinusoidal
currinSinFunExample

%% Left peak function
LeftPeakFunctionExample

%% Right peak function
RightPeakFunExample

%% Finish up
toc(StartTime)
diary off


