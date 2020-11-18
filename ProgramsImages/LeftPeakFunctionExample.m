%% Left Peak Function Example
tic
leftPeakEx = FunctionApproxProblem(repmat({@LeftPeakFun},2,1));
leftPeakEx = set_prop(leftPeakEx,'Algo',{@AdaptAlgo3});
leftPeakEx = set_prop(leftPeakEx,'whDes',{'unif'});
leftPeakEx = set_prop(leftPeakEx,'theta',{[0.4 -1]});
leftPeakEx = set_prop(leftPeakEx,'n0',{5,10});
leftPeakEx = set_prop(leftPeakEx,'nmax',{20});
leftPeakEx = set_prop(leftPeakEx,'yLim',{[-0.4; 0.5]});
leftPeakEx = set_prop(leftPeakEx,'legendPos',{'northeast'});
leftPeakEx = set_prop(leftPeakEx,'kernelOrig',{@MaternKernelGeneral});
%%
OutLeftPeakGauss = RunFunAppxExample(leftPeakEx);

%%
leftPeakEx = set_prop(leftPeakEx,'kernelOrig',{@MaternKernel});
%OutLeftPeakMatern = RunFunAppxExample(leftPeakEx);
toc