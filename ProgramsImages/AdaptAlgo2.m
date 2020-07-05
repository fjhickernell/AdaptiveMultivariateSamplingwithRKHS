%% Algorithm 2 Sample location is adaptive
function [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
   AdaptAlgo2(f,kernel,xdata,fdata,xeval,feval, ...
   abstolVec,Ainf,B0,isDiagnose,colorScheme,fname,kername)
if nargin < 10, isDiagnose = false; end
nmax = 500;
ntol = size(abstolVec,1);
neval = size(xeval,1);
n0 = 1;
plotn = [0 1 3 nmax];
if isDiagnose
   [h,ploti,legendLabel] =  ...
      multiAppxDiagPrelim(plotn,ntol,xeval,feval,colorScheme);
end
itol = 1;
abstol = abstolVec(itol);
nNeed(ntol,1) = 0;
ErrBdVec(nmax,1) = 0;
trueErr(nmax,1) = 0;
InErrBars(nmax,1) = 0;
AppxNorm(nmax,1) = 0;
nstart = 0;
coli = ploti;
AXvec(nmax,1) = 0;
BXvec(nmax,1) = 0;
minNormf(nmax,1) = 0;
maxNormf = inf(nmax,1);
NeccFlag(nmax,1) = 0;
for n = 1:nmax
   gail.print_iterations(n,'n',true)
   if n == 1
      xdata(1) = seqFixedDes(1);
   else
      xdata(n) = xeval(whKX);
   end
   fdata(n) = f(xdata(n));
   [Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata(1:n,:), xeval, kernel);
   [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   [AX, BX] = ABfun(errKX,errKNull,Ainf,B0);
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
   if isDiagnose
      if n == plotn(ploti)
         [h,legendLabel,coli,nstart] =  ...
            multiAppxDiagAddData(h,legendLabel,coli,n,nstart, ...
            xdata,fdata,xeval,Appx,colorScheme,NaN);
         ploti = ploti+1;
      end
   end
   if ErrBd < abstol
      if abstol >= 0.01
         if isDiagnose
            [h,legendLabel,coli,nstart] =  ...
               multiAppxDiagAddData(h,legendLabel,coli,n,nstart, ...
               xdata,fdata,xeval,Appx,colorScheme,abstol);
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
NeccFlag = NeccFlag(2:n+1);

if isDiagnose
   multiAppxDiagFinishPlotTable ...
      (h,legendLabel,abstolVec,ErrBdVec,trueErr,InErrBars, ...
      coli,n,n0,ntol,nNeed,fname,kername,'Alg2');
end
