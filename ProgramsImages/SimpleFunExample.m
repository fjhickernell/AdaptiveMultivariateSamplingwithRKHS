%% Example of easy function
tic
S = struct('type','{}','subs',{{':'}}); 
simpleFunEx = FunctionApproxProblem({@simpleFun,@simpleFun});
class_d = length(simpleFunEx);
simpleFunEx(2).Algo = @AdaptAlgo2;
[simpleFunEx.whDes] = subsref(repmat({'unifChebyshev'},1,2),S);
[simpleFunEx.n0] = subsref(repmat({5},1,2),S);

%%
simpleFunOutUnifCheby = RunFunAppxExample(simpleFunEx)

%% With uniform grid 
[simpleFunEx.whDes] = subsref(repmat({'unif_grid'},1,class_d),S);
simpleFunOutUnifGrid = RunFunAppxExample(simpleFunEx)

toc


   