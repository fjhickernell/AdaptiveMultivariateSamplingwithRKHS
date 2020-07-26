%% Algorithm 2 Sample location is adaptive
function [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
   AdaptAlgo2(f, kernel, xeval, feval, abstolVec, prm)
xdata(prm.nmax,1) = 0;
fdata(prm.nmax,1) = 0;
ntol = size(abstolVec,1);
neval = size(xeval,1);
plotn = [0 1 3 prm.nmax];
if prm.isDiagnose
   [h,ploti,legendLabel] =  ...
      multiAppxDiagPrelim(plotn,ntol,xeval,feval,prm);
   coli = ploti;
end
itol = 1;
abstol = abstolVec(itol);
nNeed(ntol,1) = 0;
ErrBdVec(prm.nmax,1) = 0;
trueErr(prm.nmax,1) = 0;
InErrBars(prm.nmax,1) = 0;
AppxNorm(prm.nmax,1) = 0;
nstart = 0;
AXvec(prm.nmax,1) = 0;
BXvec(prm.nmax,1) = 0;
minNormf(prm.nmax,1) = 0;
maxNormf = inf(prm.nmax,1);
NeccFlag(prm.nmax,1) = 0;
for n = prm.n0:prm.nmax
   gail.print_iterations(n,'n',true)
   if n == prm.n0
      if strcmp(prm.whDes,'uniform') && n > 1
         xdata = (0:prm.n0-1)'/(prm.n0-1);
      else %sequential
         xdata(1:prm.n0) = seqFixedDes(1:prm.n0);
      end
   else
      xdata(n) = xeval(whKX);
   end
   fdata(n) = f(xdata(n));
   [Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata(1:n,:), xeval, kernel);
   [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   [AX, BX] = ABfun(errKX,errKNull,prm.Ainf,prm.B0);
   AXvec(n) = AX;
   BXvec(n) = BX;
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
      abstol = abstolVec(itol);
   end
end
fprintf('\n')
ErrBdVec = ErrBdVec(1:n);
trueErr = trueErr(1:n);
InErrBars = InErrBars(1:n);
NeccFlag = NeccFlag(2:min(prm.nmax,n+1));

if prm.isDiagnose
   multiAppxDiagFinishPlotTable ...
      (h,legendLabel,abstolVec,ErrBdVec,trueErr,InErrBars, ...
      coli,n,ntol,nNeed,prm,xeval,'Alg2');
end
