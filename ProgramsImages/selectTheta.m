function [realTheta,prm] = selectTheta(prm,xdata,ydata,xeval)
% Here theta might be some transformation of theta
%objectKy = @(K,y) mean(log(max(eig(K),100*eps))) + log(y'*(K\y));
thetaTest = prm.thetaRange;
nth = size(thetaTest,1);
if prm.canvasTheta
   objectThTest(nth,1) = 0;
   for i = 1:nth
      objectThTest(i) = objectTh(thetaTest(i,:),xdata,ydata,xeval,prm);
   end
   [~,wh] = min(objectThTest);
   theta0 = thetaTest(wh,:);
else
   theta0 = prm.currentTheta;
end
thetaOptim = fminsearch(@(th) objectTh(th,xdata,ydata,xeval,prm),theta0);
prm.currentTheta = thetaOptim;
[~,~,~,realTheta] = prm.kernelth(xdata(1,:),xdata(1,:),thetaOptim);

function obj = objectTh(theta,xdata,ydata,xeval,prm)
kernel = @(t,x) prm.kernelth(t,x,theta);
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
      

