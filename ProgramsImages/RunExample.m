function RunExample(f,param,abstolVec,kernelth,AlgName)


%% Example for a function

gail.InitializeDisplay
format short e
warning('off')

[~,~,xeval,~,Ainf,B0] = StdParam;
p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'colorScheme',[MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon])
addParameter(p,'isDiagnose',true)
addParameter(p,'Ainf',Ainf)
addParameter(p,'B0',B0)
addParameter(p,'n0',1)
addParameter(p,'nmax',500)
addParameter(p,'theta',1)
addParameter(p,'fname','')
addParameter(p,'kername','')
addParameter(p,'legendPos','south')
addParameter(p,'whDes','adapt_th')
addParameter(p,'whobj','EmpBayesAx')
addParameter(p,'plotSites',false)
addParameter(p,'xLim',[0;1])
addParameter(p,'yLim',[-1;1])
addParameter(p,'thetaRange',(-5:0.5:5)')

parse(p,param(1));
prm = p.Results;

kernel = @(t,x) kernelth(t,x,prm.theta);
feval = f(xeval);

figure %plot function
plot(xeval,feval);
axis([prm.xLim', prm.yLim'])
xlabel('\(x\)')
ylabel('\(f(x)\)');
print('-depsc',[prm.fname 'Plot.eps'])

for kk = 1:length(AlgName)
   parse(p,param(kk));
   prm = p.Results;
   if strcmp(AlgName{kk},'Algo1')
      disp([prm.fname ' ' prm.kername ' Algorithm 1'])
      [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
         AdaptAlgo1(f, kernel, xeval, feval, abstolVec, prm);
   elseif strcmp(AlgName{kk},'Algo2')
      disp([prm.fname ' ' prm.kername ' Algorithm 2'])
      [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
         AdaptAlgo2(f, kernel, xeval, feval, abstolVec, prm);
   elseif strcmp(AlgName{kk},'Algo3')
      disp([prm.fname ' ' prm.kername ' Algorithm 3'])
      [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm] = ...
         AdaptAlgo3(f, kernelth, xeval, feval, abstolVec, prm);
   end
   fprintf(1,'\n\n')
   disp(['Necessary condition flag = ' int2str(NeccFlag(end))])
   fprintf(1,'\n\n')
end
   