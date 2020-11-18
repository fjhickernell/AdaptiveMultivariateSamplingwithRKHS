%% Example of easy function
tic
clearvars
simpleFunEx = FunctionApproxProblem({@simpleFun});
simpleFunEx = set_prop(simpleFunEx,'kernelOrig',{@MaternKernelGeneral});
simpleFunEx = set_prop(simpleFunEx,'Algo',{@AdaptAlgo3});
simpleFunEx = set_prop(simpleFunEx,'theta',{[0.4,-1]});
simpleFunEx = set_prop(simpleFunEx,'n0',{5});

%%
simpleFunOutUnif = RunFunAppxExample(simpleFunEx);

%% With uniform grid 
simpleFunEx = set_prop(simpleFunEx,'whDes',{'unif_grid'});
simpleFunOutUnifGrid = RunFunAppxExample(simpleFunEx);

toc


   