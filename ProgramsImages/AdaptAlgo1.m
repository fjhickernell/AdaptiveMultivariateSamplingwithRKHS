%% Algorithm 1 Sample size is adaptive
function [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag, prm] = ...
   AdaptAlgo1(f, kernel, xeval, feval, prm)
xdata(prm.nmax,1) = 0;
fdata(prm.nmax,1) = 0;
ntol = size(prm.abstolVec,1);
neval = size(xeval,1);
plotn = [0 1 3 prm.nmax];
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
nstart = 0;
AXvec(prm.nmax,1) = 0;
minNormf(prm.nmax,1) = 0;
maxNormf = inf(prm.nmax,1);
NeccFlag(prm.nmax,1) = 0;
for n = 1:prm.nmax
   gail.print_iterations(n,'n',true)
   xdata(n) = seqFixedDes(n);
   fdata(n) = f(xdata(n));
   [Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata(1:n,:), xeval, kernel);
   [errKXx, errKX] = powerfun(Kmat, Kdateval, Kdiageval);
   AX = ABfun(errKX,errKNull,prm.Ainf,prm.B0);
   AXvec(n) = AX;
   [Appx, AppxNorm(n), ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX );
   minNormf(n+1) = max(minNormf(n),AppxNorm(n));
   maxNormf(n+1) = min(maxNormf(n),sqrt(1+AX.^2)*AppxNorm(n));
   if minNormf(n+1) > maxNormf(n+1)
      NeccFlag(n+1) = NeccFlag(n) + 1;
      minNormf(n+1) = AppxNorm(n);
      maxNormf(n+1) = sqrt(1+AX.^2)*AppxNorm(n);
   end
   ErrBdVec(n) = ErrBd;
   trueErr(n) = max(abs(feval - Appx));
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
AppxNorm = AppxNorm(1:n);
NeccFlag = NeccFlag(2:n+1);


if prm.isDiagnose
   multiAppxDiagFinishPlotTable ...
      (h,legendLabel,ErrBdVec,trueErr,InErrBars, ...
      coli,n,ntol,nNeed,prm,xeval,'Alg1');
end
