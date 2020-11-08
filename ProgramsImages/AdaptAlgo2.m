%% Algorithm 2 Sample location is adaptive, but theta is not
function [OutObj,AlgoName] = AdaptAlgo2(f, kernel, xeval, feval, obj)
if nargout > 1
   OutObj = [];
   AlgoName = 'AdaptAlgo2';
   return
end
OutObj = FunAppxOut(obj);
d = obj.dim;
xdata(obj.nmax,1) = 0;
fdata(obj.nmax,1) = 0;
ntol = size(obj.abstolVec,1);
neval = size(xeval,1);
plotn = [obj.n0-1:5 obj.nmax];
if obj.isDiagnose
   [h,ploti,legendLabel] =  ...
      multiAppxDiagPrelim(plotn,ntol,xeval,feval,obj);
   coli = ploti;
end
itol = 1;
abstol = obj.abstolVec(itol);
nNeed(ntol,1) = 0;
OutObj.ErrBdVec(obj.nmax,1) = 0;
OutObj.trueErr(obj.nmax,1) = 0;
OutObj.InErrBars(obj.nmax,1) = 0;
OutObj.AppxNorm(obj.nmax,1) = 0;
nstart = 0;
AXvec(obj.nmax,1) = 0;
BXvec(obj.nmax,1) = 0;
minNormf(obj.nmax,1) = 0;
maxNormf = inf(obj.nmax,1);
OutObj.NeccFlag(obj.nmax,1) = 0;
for n = obj.n0:obj.nmax
   print_iterations(n,'n',true)
   if n == obj.n0
      if strcmp(obj.whDes,'uniform') && n > 1
         xdata(1:n,:) = (0:n-1)'/(n-1);
      elseif strcmp(obj.whDes,'unifChebyshev') && n > 1 %Chebyshev
         xdata(1:n,:) = (0:n-1)'/(n-1);
         xdata(1:n,:) = (1+sin(pi*(-1/2 + xdata(1:n,:))))/2;
      elseif strcmp(obj.whDes,'seqChebyshev') %sequential Chebyshev
         xdata(1:n,:) = seqFixedDes(1:n,d,1/3,'Chebyshev');
      else %sequential
         xdata(1:n,:) = seqFixedDes(1:n);
      end
      fdata(1:n) = f(xdata(1:obj.n0,:));
   else
      xdata(n) = xeval(whKX);
      fdata(n) = f(xdata(n));
   end
   [Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata(1:n,:), xeval, kernel);
   [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   [AX, BX] = ABfun(errKX,errKNull,obj.Ainf,obj.B0);
   AXvec(n) = AX;
   BXvec(n) = BX;
   [Appx, OutObj.AppxNorm(n), ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX );
   minNormf(n+1) = max(minNormf(n),OutObj.AppxNorm(n));
   maxNormf(n+1) = min(maxNormf(n),sqrt(1+AX.^2)*OutObj.AppxNorm(n));
   if minNormf(n+1) > maxNormf(n+1)
      OutObj.NeccFlag(n+1) = OutObj.NeccFlag(n) + 1;
      minNormf(n+1) = OutObj.AppxNorm(n);
      maxNormf(n+1) = sqrt(1+AX.^2)*OutObj.AppxNorm(n);
   end
   OutObj.ErrBdVec(n) = ErrBd;
   OutObj.trueErr(n) = max(abs(feval - Appx));
   errFudge = eps*cond(Kmat);
   OutObj.InErrBars(n) = sum(abs(feval - Appx) <= ErrBdx + errFudge)/neval;
   if obj.isDiagnose
      if any(plotn==n)
         [h,legendLabel,coli,nstart] =  ...
            multiAppxDiagAddData(h,legendLabel,coli,n,nstart, ...
            xdata,fdata,xeval,Appx,obj,NaN);
         ploti = ploti+1;
      end
   end
   if ErrBd < abstol
      if abstol >= 0.01
         if obj.isDiagnose
            [h,legendLabel,coli,nstart] =  ...
               multiAppxDiagAddData(h,legendLabel,coli,n,nstart, ...
               xdata,fdata,xeval,Appx,obj,abstol);
            if obj.plotSites
               plot(xdata(1:n),zeros(n,1),'k.')
            end
         end
      end
      nNeed(itol) = n;
      itol = itol + 1;
      if itol > ntol, break, end
      abstol = obj.abstolVec(itol);
   end
end
fprintf('\n')
OutObj.Appx = Appx;
OutObj.ErrBdVec = OutObj.ErrBdVec(1:n);
OutObj.trueErr = OutObj.trueErr(1:n);
OutObj.InErrBars = OutObj.InErrBars(1:n);
OutObj.NeccFlag = OutObj.NeccFlag(2:min(obj.nmax,n+1));

if obj.isDiagnose
   multiAppxDiagFinishPlotTable ...
      (h,legendLabel,OutObj,coli,n,ntol,nNeed,obj,xeval);
end
