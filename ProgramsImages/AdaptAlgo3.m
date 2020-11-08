%% Algorithm 3 Sample location is adaptive and theta is adaptive
function [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag,err,xdata,fdata,prm] = ...
   AdaptAlgo3(f, xeval, feval, prm)
[neval,d] = size(xeval);
xdata(prm.nmax,d) = 0;
fdata(prm.nmax,1) = 0;
ntol = size(prm.abstolVec,1);
plotn = [0 prm.n0 prm.nmax];
if prm.isDiagnose
   [h,ploti,legendLabel] =  ...
      multiAppxDiagPrelim(plotn,ntol,xeval,feval,prm);
   coli = ploti;
end
itol = 1;
abstol = prm.abstolVec(itol);
nNeed(ntol,1) = 0;
ErrBdVec(prm.nmax,1) = 0;
trueErr(prm.nmax,1) = 0;
InErrBars(prm.nmax,1) = 0;
AppxNorm(prm.nmax,1) = 0;
NeccFlag = NaN;
nstart = 0;
AXvec(prm.nmax,1) = 0;
BXvec(prm.nmax,1) = 0;
dth = size(prm.thetaRange,2);
thOptimVec(prm.nmax,dth) = 0;
for n = prm.n0:prm.nmax
   gail.print_iterations(n,'n',true)
   if n == prm.n0
      if strcmp(prm.whDes,'unif_grid') && n > 1 && d == 1
         xdata(1:n,:) = (0:n-1)'/(n-1);
      elseif strcmp(prm.whDes,'unifChebyshev') && n > 1 %Chebyshev
         xdata(1:n,:) = (0:n-1)'/(n-1);
         xdata(1:n,:) = (1+sin(pi*(-1/2 + xdata(1:n,:))))/2;
      elseif strcmp(prm.whDes,'seqChebyshev') %sequential Chebyshev
         xdata(1:n,:) = seqFixedDes(1:n,d,1/3,'Chebyshev');
      elseif strcmp(prm.whDes,'adapt_th') && n > 1
         kernel = @(t,x) prm.kernelth(t,x,prm.theta);
         xdata(1,:) = seqFixedDes(1,d);
          for ii = 2:prm.n0
            [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:ii,:), xeval, kernel);
            [~, ~, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
            xdata(ii,:) = xeval(whKX,:);
         end
      else %sequential, the default
         xdata(1:n,:) = seqFixedDes(1:n,d);
      end
      fdata(1:n) = f(xdata(1:n,:));
   else
      xdata(n,:) = xeval(whKX,:);
      fdata(n) = f(xdata(n,:));
   end
   [thOptim,prm] = selectTheta(prm,xdata(1:n,:),fdata(1:n), ...
      xeval);
   thOptimVec(n,:) = thOptim;
   kernel = @(t,x) prm.kernelth(t,x,prm.currentTheta);
   [Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata(1:n,:), xeval, kernel);
   [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   [AX, BX] = ABfun(errKX,errKNull,prm.Ainf,prm.B0);
   AXvec(n) = AX;
   BXvec(n) = BX;
   [Appx, AppxNorm(n), ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX );
   ErrBdVec(n) = ErrBd;
   err = abs(feval-Appx);
   trueErr(n) = max(err);
   errFudge = eps*cond(Kmat);
   InErrBars(n) = sum(abs(feval - Appx) <= ErrBdx + errFudge)/neval;
   if prm.isDiagnose
      if n == plotn(ploti)
         [h,legendLabel,coli,nstart] =  ...
            multiAppxDiagAddData(h,legendLabel,coli,n,nstart, ...
            xdata,fdata,xeval,Appx,prm,NaN);
         ploti = ploti+1;
      end
   end
   if ErrBd < abstol
      if abstol >= 0.01
         if prm.isDiagnose
            [h,legendLabel,coli,nstart] =  ...
               multiAppxDiagAddData(h,legendLabel,coli,n,nstart, ...
               xdata,fdata,xeval,Appx,prm,abstol);
            if prm.plotSites
               plot(xdata(1:n),zeros(n,1),'k.')
            end
         end
      end
      nNeed(itol) = n;
      itol = itol + 1;
      if itol > ntol, break, end
      abstol = prm.abstolVec(itol);
   end
end
fprintf('\n')
ErrBdVec = ErrBdVec(1:n);
trueErr = trueErr(1:n);
InErrBars = InErrBars(1:n);
xdata = xdata(1:n,:);
fdata = fdata(1:n);
prm.final_theta = thOptim;

if prm.isDiagnose
   multiAppxDiagFinishPlotTable ...
      (h,legendLabel,ErrBdVec,trueErr,InErrBars, ...
      coli,n,ntol,nNeed,prm,xeval,'Alg3');
end
