%% Function approximation examples

%% Start
clearvars
warning('off')

%% Open Diary File
Outfile = 'FunctionApproxExOut.txt'; %open a diary file to capture command window print out
fopen(Outfile,'w+');
diary(Outfile)
StartTime = tic; %start time

%% Simple function
SimpleFunExample

%% Currin sinusoidal
currinSinFunExample

%% Left peak function
LeftPeakFunctionExample

%% Right peak function
RightPeakFunExample

%%% Lim function
LimFunExample

%% Finish up
toc(StartTime)
diary off


