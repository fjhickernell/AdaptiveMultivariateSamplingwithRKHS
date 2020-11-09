%% Left Peak Function Example
leftPeakEx = FunctionApproxProblem(repmat({@LeftPeakFun},1,2));
leftPeakEx = set_prop(leftPeakEx,'Algo',{@AdaptAlgo2,@AdaptAlgo3});
leftPeakEx = set_prop(leftPeakEx,'n0',{1,10});
leftPeakEx = set_prop(leftPeakEx,'yLim',{[-0.5;0.5]});
leftPeakEx = set_prop(leftPeakEx,'legendPos',{'northeast'});
%%
OutLeftPeak = RunFunAppxExample(leftPeakEx);
