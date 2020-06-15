function thetaOptim = selectTheta(thetaTest,kernelth,kerneldiag,xdata,ydata,xeval,errKNull,Ainf,B0)
% Here theta might be some transformation of theta
%objectKy = @(K,y) mean(log(max(eig(K),100*eps))) + log(y'*(K\y));
nth = size(thetaTest,1);
objectThTest(nth,1) = 0;
for i = 1:nth
   objectThTest(i) = objectTh(thetaTest(i,:),kernelth,kerneldiag,xdata,ydata,xeval,errKNull,Ainf,B0);
end
[smallObjTh,wh] = min(objectThTest);
theta0 = thetaTest(wh,:);
thetaOptim = fminsearch(@(th) objectTh(th,kernelth,kerneldiag,xdata,ydata,xeval,errKNull,Ainf,B0),theta0);

function obj = objectTh(theta,kernelth,kerneldiag,xdata,ydata,xeval,errKNull,Ainf,B0)
kernel = @(t,x) kernelth(t,x,theta);
[Kmat, Kdateval, Kdiageval] = KMP(xdata, xeval, kernel, kerneldiag);
[~, errKX] = powerfun(Kmat, Kdateval, Kdiageval);
AXtheta = ABfun(errKX,errKNull,Ainf,B0);
obj = mean(log(max(eig(Kmat),100*eps))) + log(ydata'*(Kmat\ydata)) ...
   + log(1+AXtheta.^2);

