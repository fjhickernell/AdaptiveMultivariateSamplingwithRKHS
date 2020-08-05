function RunExample(f,param,kernelth)

[~,~,xeval] = StdParam;

feval = f(xeval);
prm = param(1);
figure %plot function
plot(xeval,feval);
axis([prm.xLim', prm.yLim'])
xlabel('\(x\)')
ylabel('\(f(x)\)');
print('-depsc',[prm.fname 'Plot.eps'])

for kk = 1:size(param,2)
   prm = param(kk);
   kernel = @(t,x) kernelth(t,x,param(kk).theta);
   disp([prm.fname ' ' prm.kername ' ' prm.whDes ' ' prm.whObj ' ' prm.AlgName])
   if strcmp(prm.AlgName,'Algo1')
      [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
         AdaptAlgo1(f, kernel, xeval, feval, prm);
   elseif strcmp(prm.AlgName,'Algo2')
      [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
         AdaptAlgo2(f, kernel, xeval, feval, prm);
   elseif strcmp(prm.AlgName,'Algo3')
      [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, prm] = ...
         AdaptAlgo3(f, kernelth, xeval, feval, prm);
   end
   if any(strcmp(prm.AlgName,{'Algo1','Algo2'}))
      disp(['Necessary condition flag = ' int2str(NeccFlag(end))])
   end
   fprintf(1,'\n\n')
end
   