function [realTheta,currentTheta] = selectTheta(obj,xdata,ydata,xeval,currentTheta)
% Here theta might be some transformation of theta
%objectKy = @(K,y) mean(log(max(eig(K),100*eps))) + log(y'*(K\y));
thetaTest = obj.thetaRange;
nth = size(thetaTest,1);
if obj.canvasTheta
   objectThTest(nth,1) = 0;
   for i = 1:nth
      objectThTest(i) = objectTh(thetaTest(i,:),xdata,ydata,xeval,obj);
   end
   [~,wh] = min(objectThTest);
   theta0 = thetaTest(wh,:);
else
   theta0 = currentTheta;
end
thetaOptim = fminsearch(@(th) objectTh(th,xdata,ydata,xeval,obj),theta0, ...
   optimset('MaxIter',20,'TolFun',1e-3,'TolX',1e-3,'Display','none'));
currentTheta = thetaOptim;
[~,~,~,realTheta] = obj.kernelth(xdata(1,:),xdata(1,:),thetaOptim);

function object = objectTh(theta,xdata,ydata,xeval,obj)
kernel = @(t,x) obj.kernelth(t,x,theta);
[Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata, xeval, kernel);
[~, errKX] = powerfun(Kmat, Kdateval, Kdiageval);
AXtheta = ABfun(errKX,errKNull,obj.Ainf,obj.B0);
switch obj.whObj
   case 'EmpBayes'
      object = mean(log(max(eig(Kmat),100*eps))) + log(ydata'*(Kmat\ydata));
   case 'minErrBd'
      object = 2*(log(AXtheta) + log(errKX)) + log(ydata'*(Kmat\ydata));
   otherwise %'EmpBayesAx'
      %object = mean(log(max(eig(Kmat),100*eps))) + log(ydata'*(Kmat\ydata)) ...
      object = mean(log(max(eig(Kmat),100*eps))) + log(ydata'*KinvY(Kmat,ydata)) ...
         + log(1+AXtheta.^2);
end
      

