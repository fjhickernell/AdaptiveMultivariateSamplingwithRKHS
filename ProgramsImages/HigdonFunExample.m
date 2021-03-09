%% Example of FORRESTER ET AL. (2008) FUNCTION function
tic
clearvars
HigdonFunEx = FunctionApproxProblem(repmat({@HigdonFun},1,3));
%simpleFunEx = set_prop(simpleFunEx,'kernelOrig',{@MaternKernel});
HigdonFunEx = set_prop(HigdonFunEx,'Algo',{@AdaptAlgo1,@AdaptAlgo2,@AdaptAlgo3});
HigdonFunEx = set_prop(HigdonFunEx,'theta',{0, 0, 1});
HigdonFunEx = set_prop(HigdonFunEx,'n0',repmat({5},1,3));
HigdonFunEx = set_prop(HigdonFunEx,'nmax',repmat({100},1,3));
HigdonFunEx = set_prop(HigdonFunEx,'yLim',repmat({[-1.5;1.5]},1,3));
HigdonFunEx = set_prop(HigdonFunEx,'abstolVec',repmat({[0.01 0.001 0.0001 0.00001]'},1,3));

%%
[OutHigdonFunGaussOptimTheta,HigdonFunEx] = RunFunAppxExample(HigdonFunEx);

toc

%%
clearvars
HigdonFunEx = FunctionApproxProblem({@HigdonFun});
%simpleFunEx = set_prop(simpleFunEx,'kernelOrig',{@MaternKernel});
HigdonFunEx = set_prop(HigdonFunEx,'Algo',{@AdaptAlgo1});
HigdonFunEx = set_prop(HigdonFunEx,'theta',{0});
HigdonFunEx = set_prop(HigdonFunEx,'n0',{1});
HigdonFunEx = set_prop(HigdonFunEx,'nmax',{100});
HigdonFunEx = set_prop(HigdonFunEx,'yLim',{[-1.5;1.5]});
HigdonFunEx = set_prop(HigdonFunEx,'abstolVec',{[0.01 0.001 0.0001 0.00001]'});

%%
[OutHigdonFunGaussOptimTheta,HigdonFunEx] = RunFunAppxExample(HigdonFunEx);

   