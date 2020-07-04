%% Algorithm 3 Sample location is adaptive
function [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars] = ...
   AdaptAlgo3(f,kernelth,xdata,fdata,xeval,feval, ...
   abstolVec,Ainf,B0,n0,thetaRange,isDiagnose,colorScheme,fname,kername)
if nargin < 12, isDiagnose = false; end
nmax = 500;
ntol = size(abstolVec,1);
neval = size(xeval,1);
plotn = [0 n0 nmax];
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
nstart = 0;
coli = ploti;
AXvec(nmax,1) = 0;
BXvec(nmax,1) = 0;
thOptimVec(nmax,1) = 0;
for n = n0:nmax
   gail.print_iterations(n,'n',true)
   if n == n0
      xdata(1:n0) = seqFixedDes(1:n0);
      fdata(1:n0) = f(xdata(1:n0));
   else
      xdata(n) = xeval(whKX);
      fdata(n) = f(xdata(n));
   end
   lnthOptim = selectTheta(thetaRange,kernelth,xdata(1:n),fdata(1:n), ...
      xeval,Ainf,B0);
   thetaOptim = exp(lnthOptim);
   thOptimVec(n) = thetaOptim;
   kernel = @(t,x) kernelth(t,x,lnthOptim);
   [Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata(1:n,:), xeval, kernel);
   [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   [AX, BX] = ABfun(errKX,errKNull,Ainf,B0);
   AXvec(n) = AX;
   BXvec(n) = BX;
   [Appx, ~, ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX );
   ErrBdVec(n) = ErrBd;
   trueErr(n) = max(abs(feval - Appx));
   errFudge = eps*cond(Kmat);
   InErrBars(n) = sum(abs(feval - Appx) <= ErrBdx + errFudge)/neval;
   if isDiagnose
      if n == plotn(ploti)
         [h,legendLabel,coli,nstart] =  ...
            multiAppxDiagAddData(h,legendLabel,coli,n,nstart, ...
            xdata,fdata,xeval,Appx,colorScheme);
         ploti = ploti+1;
      end
   end
   if ErrBd < abstol
      if abstol >= 0.01
         if isDiagnose
            [h,legendLabel,coli,nstart] =  ...
               multiAppxDiagAddData(h,legendLabel,coli,n,nstart, ...
               xdata,fdata,xeval,Appx,colorScheme);
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

if isDiagnose
   multiAppxDiagFinishPlotTable ...
      (h,legendLabel,abstolVec,ErrBdVec,trueErr,InErrBars, ...
      coli,n,n0,ntol,nNeed,fname,kername,'Alg3');
end
