%% Exponential times Cosine of Sum Example
tic
clearvars
ExpCosSumFunEx = FunctionApproxProblem(repmat({@ExpCosSumFun},1,2));
ExpCosSumFunEx = set_prop(ExpCosSumFunEx,'fparam',{[1 0.1; 1 0.1]});
ExpCosSumFunEx = set_prop(ExpCosSumFunEx,'whDes',{'apdapt_th','unifChebyshev'});
ExpCosSumFunEx = set_prop(ExpCosSumFunEx,'nmax',{200});
ExpCosSumFunEx = set_prop(ExpCosSumFunEx,'theta',{[-1 -1]});
ExpCosSumFunEx = set_prop(ExpCosSumFunEx,'xLim',{[0 0; 1 1]});
ExpCosSumFunEx = set_prop(ExpCosSumFunEx,'Algo',{@AdaptAlgo3});
ExpCosSumFunEx = set_prop(ExpCosSumFunEx,'n0',{10});
ExpCosSumFunEx = set_prop(ExpCosSumFunEx,'yLim',{[-1;1]});
ExpCosSumFunEx = set_prop(ExpCosSumFunEx,'legendPos',{'eastoutside'});
%%
[OutExpCosSumGauss,ExpCosSumFunEx] = RunFunAppxExample(ExpCosSumFunEx);

ExpCosSumFunEx = set_prop(ExpCosSumFunEx,'kernelOrig',{@MaternKernel});
%[OutCosSumGaussMatern,cosSumFunEx] = RunFunAppxExample(cosSumFunEx);
