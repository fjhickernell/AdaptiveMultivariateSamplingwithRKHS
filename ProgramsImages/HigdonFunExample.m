%% Example of FORRESTER ET AL. (2008) FUNCTION function
tic
clearvars
HigdonFunEx = FunctionApproxProblem(repmat({@HigdonFun},1,3));
%simpleFunEx = set_prop(simpleFunEx,'kernelOrig',{@MaternKernel});
HigdonFunEx = set_prop(HigdonFunEx,'Algo',{@AdaptAlgo1,@AdaptAlgo2,@AdaptAlgo3});
HigdonFunEx = set_prop(HigdonFunEx,'theta',{2, 2, 2});
HigdonFunEx = set_prop(HigdonFunEx,'n0',{10,10,10});
HigdonFunEx = set_prop(HigdonFunEx,'nmax',repmat({100},1,3));
HigdonFunEx = set_prop(HigdonFunEx,'yLim',repmat({[-1.5;1.5]},1,3));
HigdonFunEx = set_prop(HigdonFunEx,'abstolVec',repmat({[0.01 0.001 0.0001 0.00001]'},1,3));

%%
[OutHigdonFunGaussOptimTheta,HigdonFunEx] = RunFunAppxExample(HigdonFunEx);

toc

%%