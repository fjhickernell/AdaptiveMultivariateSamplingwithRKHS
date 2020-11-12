%% Cosine of Sum Example
tic
clearvars
cosSumFunEx = FunctionApproxProblem(repmat({@CosSumFun},1,3));
d = 12;
fpar = 2.^(-1:-1:-d);
fpar = fpar(1,randperm(d));
cosSumFunEx = set_prop(cosSumFunEx,'fparam',{fpar});
cosSumFunEx = set_prop(cosSumFunEx,'whDes',{'apdapt_th','unifChebyshev','seqChebyshev'});
cosSumFunEx = set_prop(cosSumFunEx,'nmax',{200});
cosSumFunEx = set_prop(cosSumFunEx,'theta',{-ones(1,d)});
cosSumFunEx = set_prop(cosSumFunEx,'xLim',{[zeros(1,d); ones(1,d)]});
cosSumFunEx = set_prop(cosSumFunEx,'Algo',{@AdaptAlgo3});
cosSumFunEx = set_prop(cosSumFunEx,'n0',{10});
cosSumFunEx = set_prop(cosSumFunEx,'yLim',{[-1;1]});
cosSumFunEx = set_prop(cosSumFunEx,'legendPos',{'eastoutside'});
%%
[OutCosSumGauss,cosSumFunEx] = RunFunAppxExample(cosSumFunEx);

cosSumFunEx = set_prop(cosSumFunEx,'kernelOrig',{@MaternKernel});
%[OutCosSumGaussMatern,cosSumFunEx] = RunFunAppxExample(cosSumFunEx);
