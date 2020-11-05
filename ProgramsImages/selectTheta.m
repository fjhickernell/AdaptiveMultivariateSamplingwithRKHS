function [realTheta,thetaOptim] = ...
   selectTheta(thetaTest,kernelth,xdata,ydata,xeval,prm)
% Here theta might be some transformation of theta
%objectKy = @(K,y) mean(log(max(eig(K),100*eps))) + log(y'*(K\y));
nth = size(thetaTest,1);
if prm.canvasTheta
   objectThTest(nth,1) = 0;
   for i = 1:nth
      objectThTest(i) = objectTh(thetaTest(i,:),kernelth,xdata,ydata,xeval,prm);
   end
   [~,wh] = min(objectThTest);
   theta0 = thetaTest(wh,:);
else
   theta0 = prm.currentTheta;
end
thetaOptim = fminsearch(@(th) objectTh(th,kernelth,xdata,ydata,xeval,prm),theta0);
d = size(xdata,2);
[~,~,~,realTheta] = kernelth(xdata(1,:),xdata(1,:),thetaOptim);

function obj = objectTh(theta,kernelth,xdata,ydata,xeval,prm)
kernel = @(t,x) kernelth(t,x,theta);
[Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata, xeval, kernel);
[~, errKX] = powerfun(Kmat, Kdateval, Kdiageval);
AXtheta = ABfun(errKX,errKNull,prm.Ainf,prm.B0);
switch prm.whObj
   case 'EmpBayes'
      obj = mean(log(max(eig(Kmat),100*eps))) + log(ydata'*(Kmat\ydata));
   case 'minErrBd'
      obj = 2*(log(AXtheta) + log(errKX)) + log(ydata'*(Kmat\ydata));
   otherwise %'EmpBayesAx'
      obj = mean(log(max(eig(Kmat),100*eps))) + log(ydata'*(Kmat\ydata)) ...
         + log(1+AXtheta.^2);
end
      

