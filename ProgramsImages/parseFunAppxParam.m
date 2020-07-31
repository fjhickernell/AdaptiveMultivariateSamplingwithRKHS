function [prm,kernelth] = parseFunAppxParam(param)

kernelth = @(t,x,theta) MaternKernel(t,x,theta,true);

gail.InitializeDisplay
[~,~,~,~,Ainf,B0] = StdParam;
p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'AlgName','Algo1')
addParameter(p,'colorScheme',[MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon])
addParameter(p,'isDiagnose',true)
addParameter(p,'Ainf',Ainf)
addParameter(p,'B0',B0)
addParameter(p,'n0',1)
addParameter(p,'nmax',500)
addParameter(p,'theta',1)
addParameter(p,'fname','')
addParameter(p,'kername','Matern')
addParameter(p,'legendPos','south')
addParameter(p,'whDes','adapt_th')
addParameter(p,'whObj','EmpBayesAx')
addParameter(p,'plotSites',false)
addParameter(p,'xLim',[0;1])
addParameter(p,'yLim',[-1;1])
addParameter(p,'thetaRange',(-5:0.5:5)')
addParameter(p,'abstolVec',[0.05 0.02 0.01 0.005 0.002 0.001]')

dimParam = size(param,2);
parse(p,param(dimParam));
prm(dimParam) = p.Results;
for kk = 1:dimParam-1
   parse(p,param(kk));
   prm(kk) = p.Results;
end

