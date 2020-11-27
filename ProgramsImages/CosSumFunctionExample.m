%% Cosine of Sum Example
tic
clearvars
cosSumFunEx = FunctionApproxProblem(repmat({@CosSumFun},1,1));
rng(47)
d = 12;
fpar = 2.^(-1:-1:-d);
fpar = fpar(1,randperm(d));
cosSumFunEx = set_prop(cosSumFunEx,'fparam',{fpar});
cosSumFunEx = set_prop(cosSumFunEx,'whDes',{'lattice'});
cosSumFunEx = set_prop(cosSumFunEx,'whObj',{'EmpBayes'});
%cosSumFunEx = set_prop(cosSumFunEx,'whDes',{'apdapt_th'});
cosSumFunEx = set_prop(cosSumFunEx,'nmax',{50});
cosSumFunEx = set_prop(cosSumFunEx,'theta',{-3*ones(1,d)});
cosSumFunEx = set_prop(cosSumFunEx,'xLim',{[zeros(1,d); ones(1,d)]});
cosSumFunEx = set_prop(cosSumFunEx,'Algo',{@AdaptAlgo3});
cosSumFunEx = set_prop(cosSumFunEx,'n0',{2^4});
cosSumFunEx = set_prop(cosSumFunEx,'yLim',{[-1;1]});
cosSumFunEx = set_prop(cosSumFunEx,'legendPos',{'eastoutside'});
cosSumFunEx = set_prop(cosSumFunEx,'abstolVec',{[0.1 0.05 0.02 0.01 0.005 0.002 0.001]'});
%%
cosSumFunEx = set_prop(cosSumFunEx,'neval',{2^14});
[OutCosSumGauss,cosSumFunEx] = RunFunAppxExample(cosSumFunEx);

cosSumFunEx = set_prop(cosSumFunEx,'Algo',{@AdaptAlgo2});
cosSumFunEx = set_prop(cosSumFunEx,'theta',{-3*ones(1,d)});
cosSumFunEx = set_prop(cosSumFunEx,'abstolVec',{[0.1 0.05 0.02 0.01]'});
cosSumFunEx = set_prop(cosSumFunEx,'nmax',{50});
[OutCosSumGaussMatern,cosSumFunEx] = RunFunAppxExample(cosSumFunEx);
%  cosSumFunEx = set_prop(cosSumFunEx,'kernelOrig',{@MaternKernelGeneral});
%  cosSumFunEx = set_prop(cosSumFunEx,'theta',{[0.4 -ones(1,d)]});
%  [OutCosSumFunGeneralMatern,cosSumFunEx]= RunFunAppxExample(cosSumFunEx);
