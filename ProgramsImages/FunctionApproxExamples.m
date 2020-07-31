%% Function approximation examples

%% Open File
clearvars
Outfile = 'FunctionApproxExOut.txt'; %open a diary file to capture command window print out
fopen(Outfile,'w+');
diary(Outfile)
StartTime = tic; %start time

%% Defaults
S = struct('type','{}','subs',{{':'}});

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


