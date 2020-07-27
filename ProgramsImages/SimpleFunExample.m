%% Example of easy function

f = @simpleFun;
prm.fname = 'SimpleFun';
kernelth = @(t,x,theta) MaternKernel(t,x,theta,true);
prm.kername = 'Matern';
abstolVec = [0.05 0.02 0.01 0.005 0.002 0.001]';
AlgName = {'Algo1', 'Algo2', 'Algo3'};
RunExample(f,prm,abstolVec,kernelth,AlgName)

%%
prm.whdes = 'unif_grid';
RunExample(f,prm,abstolVec,kernelth,AlgName)


return

   