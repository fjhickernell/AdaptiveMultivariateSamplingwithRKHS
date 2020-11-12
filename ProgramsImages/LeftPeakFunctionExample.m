%% Left Peak Function Example
tic
leftPeakEx = FunctionApproxProblem(repmat({@LeftPeakFun},1,2));
leftPeakEx = set_prop(leftPeakEx,'Algo',{@AdaptAlgo3,@AdaptAlgo3});
leftPeakEx = set_prop(leftPeakEx,'theta',{-1,[-1 0]});
leftPeakEx = set_prop(leftPeakEx,'n0',{5});
leftPeakEx = set_prop(leftPeakEx,'nmax',{100});
leftPeakEx = set_prop(leftPeakEx,'yLim',{[-0.4;0.5]});
leftPeakEx = set_prop(leftPeakEx,'legendPos',{'northeast'});
%%
OutLeftPeakGauss = RunFunAppxExample(leftPeakEx);

%%
leftPeakEx = set_prop(leftPeakEx,'kernelOrig',{@MaternKernel});
%OutLeftPeakMatern = RunFunAppxExample(leftPeakEx);
toc