%% Example of easy function
tic
clearvars
simpleFunEx = FunctionApproxProblem({@simpleFun,@simpleFun});
simpleFunEx(2).Algo = @AdaptAlgo2;
simpleFunEx = set_prop(simpleFunEx,'whDes',{'unifChebyshev'});
simpleFunEx = set_prop(simpleFunEx,'n0',{5});

%%
simpleFunOutUnifCheby = RunFunAppxExample(simpleFunEx);

%% With uniform grid 
simpleFunEx = set_prop(simpleFunEx,'whDes',{'unif_grid'});
simpleFunOutUnifGrid = RunFunAppxExample(simpleFunEx);

toc


   