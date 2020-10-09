%OTL Circuit Function
%
% \[
% V_m(x) = \frac{(V_{b_1} +0.74)\beta (R_{c_2}+9)}{\beta (R_{c_2}+9) +R_f}
%         + \frac{11.35 R_f}{\beta (R_{c_2}+9) +R_f 
%         + \frac{0.74R_f\beta (R_{c_2}+9)}{(\beta (R_{c_2}+9)
%         +R_f)R_{c_1}},
% \]
% where \(V_{b_1} = 12 R_{b_2}/(R_{b_1}+R_{b_2})\).
% \(R_{b_1} \in[50, 150], R_{b_2} \in [25, 70],
%   R_f \in [0.5, 3], R_{c_1} \in [1.2, 2.5],
%   R_{c_2} \in [0.25, 1.2], \beta \in[50, 300]\).
%
% f(x) = 

S = struct('type','{}','subs',{{':'}}); prm=[];
[prm(1).AlgName] = subsref({'Algo3'},S);
[prm.kername] = subsref({'SpatialMatern'},S);
[prm,kernelth] = parseFunAppxParam(prm);
nAlg = size(prm,2);
v = @(v) (250*v(:,2)+50).*(0.95*v(:,1)+9.25);
rf = @(x) 2.5*x+0.5;
rc1 = @(x) 1.3*x+1.2;
vb1 = @(b) (24*b(:,2)+12)./(2*b(:,2)+4*b(:,1)+3);
f = @(x) ((vb1(x(:,1:2))+0.74).*v(x(:,5:6)).*rc1(x(:,4))+...
    11.35*rf(x(:,3)).*rc1(x(:,4))+ 0.74*rf(x(:,3)).*v(x(:,5:6)))...
    ./(v(x(:,5:6))+rf(x(:,3)))./rc1(x(:,4));

[prm.fname] = subsref(repmat({'OTLCircuitFun'},1,nAlg),S);
%[prm.kername] = subsref({'SpatialMatern'},S);
[prm.n0] = subsref({50},S);
[prm.theta] = subsref({[1 1 1 1 1 1 0 0 0 0 0 0]},S);
xRange = (-5:0.5:5)';
[thaa,thbb] = meshgrid(xRange,xRange);
[prm.thetaRange] = subsref({...
    [thaa(:) thaa(:) thaa(:) thaa(:) thaa(:) thaa(:) thbb(:) thbb(:) thbb(:) thbb(:) thbb(:) thbb(:)]},S);
[prm.yLim] = subsref(repmat({[-0.2;0.5]},1,nAlg),S);
[prm.legendPos] = subsref(repmat({''},1,nAlg),S);
[prm.plotSites] = subsref({false},S);
%[prm.abstolVec] = subsref({[0.1 0.05 0.02 ]'},S);
[prm.abstolVec] = subsref({1},S);
[prm.nmax] = subsref({500},S);
[prm.canvasTheta]= subsref({false},S);
[prm.currentTheta]= subsref({[0 0 0 0 0 0 0 0 0 0 0 0]},S);
[prm.whDes] = subsref({'adapt_th'},S);
[prm.whObj] = subsref({'EmpBayes'},S);
[prm.isDiagnose] = subsref({false},S);
highdimflag = 6;
%%
tic
RunExample(f,prm,kernelth,highdimflag)
toc
