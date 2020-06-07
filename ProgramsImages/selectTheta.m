function thetaOptim = selectTheta(thetaTest,kernelth,xdata,ydata)
% Here theta might be some transformation of theta
Ktheta = @(theta) kernelth(xdata,xdata,theta);
objectKy = @(K,y) mean(log(max(eig(K),100*eps))) + log(y'*(K\y));
objectTh = @(th) objectKy(Ktheta(th),ydata);
nth = size(thetaTest,1);
objectThTest(nth,1) = 0;
for i = 1:nth
   objectThTest(i) = objectTh(thetaTest(i,:));
end
[smallObjTh,wh] = min(objectThTest)
theta0 = thetaTest(wh,:);
thetaOptim = fminsearch(objectTh,theta0);
