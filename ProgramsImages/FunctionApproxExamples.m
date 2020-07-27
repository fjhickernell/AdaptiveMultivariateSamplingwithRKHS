%% Function approximation examples

%% Default
clearvars
abstolVec = [0.05 0.02 0.01 0.005 0.002 0.001]';
AlgName = {'Algo1', 'Algo2', 'Algo3'};
kernelth = @(t,x,theta) MaternKernel(t,x,theta,true);
prm(3).kername = 'Matern';
for ii = 1:2; prm(ii).kername = 'Matern'; end

%% Simple function
f = @simpleFun;
prm(3).fname = 'SimpleFun';
for ii = 1:2; prm(ii).fname = 'SimpleFun'; end
prm(3).n0 = 5;
for ii = 1:2; prm(ii).n0 = 1; end
RunExample(f,prm,abstolVec,kernelth,AlgName)

%% With uniform grid design
prm.whdes = 'unif_grid';
RunExample(f,prm,abstolVec,kernelth,AlgName)

%% Currin sinusoidal
f = @(x) sin(2*pi*(x-0.1));
prm = [];
prm.yLim = [-2;2];
prm(3).fname = 'CurrinSineFun';
for ii = 1:2; prm(ii).fname = 'SimpleFun'; end
prm(3).n0 = 5;
for ii = 1:2; prm(ii).n0 = 1; end
prm.legendPos = 'south';
RunExample(f,prm,abstolVec,kernelth,AlgName)

