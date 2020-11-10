%% Cosine of Sum Example
tic
clearvars
cosSumFunEx = FunctionApproxProblem(repmat({@CosSumFun},1,2));
cosSumFunEx = set_prop(cosSumFunEx,'fparam',{[1 0.1]});
cosSumFunEx = set_prop(cosSumFunEx,'whDes',{'apdapt_th','unifChebyshev'});
cosSumFunEx = set_prop(cosSumFunEx,'nmax',{200});
cosSumFunEx = set_prop(cosSumFunEx,'theta',{[-1 -1]});
cosSumFunEx = set_prop(cosSumFunEx,'xLim',{[0 0; 1 1]});
cosSumFunEx = set_prop(cosSumFunEx,'Algo',{@AdaptAlgo3});
cosSumFunEx = set_prop(cosSumFunEx,'n0',{10});
cosSumFunEx = set_prop(cosSumFunEx,'yLim',{[-1;1]});
cosSumFunEx = set_prop(cosSumFunEx,'legendPos',{'eastoutside'});
%%
[OutCosSumGauss,cosSumFunEx] = RunFunAppxExample(cosSumFunEx);

cosSumFunEx = set_prop(cosSumFunEx,'kernelOrig',{@MaternKernel});
%[OutCosSumGaussMatern,cosSumFunEx] = RunFunAppxExample(cosSumFunEx);
