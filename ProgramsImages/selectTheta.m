function thetaOptim = selectTheta(thetaTest,kernelth,xdata,ydata,xeval,Ainf,B0,whobj)
% Here theta might be some transformation of theta
%objectKy = @(K,y) mean(log(max(eig(K),100*eps))) + log(y'*(K\y));
if nargin < 8; whobj = 'EmpBayesAx'; end
nth = size(thetaTest,1);
objectThTest(nth,1) = 0;
for i = 1:nth
   objectThTest(i) = objectTh(thetaTest(i,:),kernelth,xdata,ydata,xeval,Ainf,B0,whobj);
end
[smallObjTh,wh] = min(objectThTest);
theta0 = thetaTest(wh,:);
thetaOptim = fminsearch(@(th) objectTh(th,kernelth,xdata,ydata,xeval,Ainf,B0,whobj),theta0);

function obj = objectTh(theta,kernelth,xdata,ydata,xeval,Ainf,B0,whobj)
kernel = @(t,x) kernelth(t,x,theta);
[Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata, xeval, kernel);
[~, errKX] = powerfun(Kmat, Kdateval, Kdiageval);
AXtheta = ABfun(errKX,errKNull,Ainf,B0);
switch whobj
   case 'EmpBayes'
      obj = mean(log(max(eig(Kmat),100*eps))) + log(ydata'*(Kmat\ydata));
   case 'minErrBd'
      obj = 2*(log(AXtheta) + log(errKX)) + log(ydata'*(Kmat\ydata));
   otherwise %'EmpBayesAx'
      obj = mean(log(max(eig(Kmat),100*eps))) + log(ydata'*(Kmat\ydata)) ...
         + log(1+AXtheta.^2);
end
      

