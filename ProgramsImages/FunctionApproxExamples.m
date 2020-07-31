%% Function approximation examples

%% Open File
clearvars
Outfile = 'FunctionApproxExOut.txt'; %open a diary file to capture command window print out
fopen(Outfile,'w+');
diary(Outfile)
StartTime = tic; %start time

%% Defaults
S = struct('type','{}','subs',{{':'}});
[prm(1:3).AlgName] = subsref({'Algo1', 'Algo2', 'Algo3'},S);
[prm,kernelth] = parseFunAppxParam(prm);
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


