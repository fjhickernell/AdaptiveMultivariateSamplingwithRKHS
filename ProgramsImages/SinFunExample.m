%% Example of sinusodial function
tic
clearvars
sinFunEx = FunctionApproxProblem(repmat({@sinFun},1,3));
%simpleFunEx = set_prop(simpleFunEx,'kernelOrig',{@MaternKernel});
sinFunEx = set_prop(sinFunEx,'Algo',{@AdaptAlgo1,@AdaptAlgo2,@AdaptAlgo3});
sinFunEx = set_prop(sinFunEx,'theta',{0 , 0 , 1});
sinFunEx = set_prop(sinFunEx,'n0',repmat({1},1,3));
sinFunEx = set_prop(sinFunEx,'nmax',repmat({50},1,3));
sinFunEx = set_prop(sinFunEx,'yLim',repmat({[-1.5;1.5]},1,3));
sinFunEx = set_prop(sinFunEx,'abstolVec',repmat({[0.01 0.001 0.0001 0.00001]'},1,3));

%%
[OutsinFunGaussOptimTheta,sinFunEx] = RunFunAppxExample(sinFunEx);

%% With uniform grid 
%sinFunEx = set_prop(sinFunEx,'whDes',{'unif_grid'});
%simpleFunOutUnifGrid = RunFunAppxExample(sinFunEx);

toc


   