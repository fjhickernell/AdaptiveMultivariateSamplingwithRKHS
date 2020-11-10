%% Algorithm 1 
% Sample size is adaptive, but sample location and theta are not
function OutObj = AdaptAlgo1(f, kernel, xeval, feval, obj)
OutObj = FunAppxOut(obj);
d = obj.dim;
xdata(obj.nmax,d) = 0;
fdata(obj.nmax,1) = 0;
ntol = size(obj.abstolVec,1);
neval = size(xeval,1);
plotn = [0 1 3 obj.nmax];
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
minNormf(obj.nmax,1) = 0;
maxNormf = inf(obj.nmax,1);
OutObj.NeccFlag(obj.nmax+1,1) = 0;
for n = obj.n0:obj.nmax
   print_iterations(n,'n',true)
   if n == obj.n0
      xdata = fixedDesign(1:n,obj);
      fdata(1:n) = f(xdata(1:n,:));
   else
      xdata(n,:) = fixedDesign(n,obj);
      fdata(n) = f(xdata(n,:));
   end
   fdata(n) = f(xdata(n));
   [Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata(1:n,:), xeval, kernel);
   [errKXx, errKX] = powerfun(Kmat, Kdateval, Kdiageval);
   AX = ABfun(errKX,errKNull,obj.Ainf,obj.B0);
   AXvec(n) = AX;
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
      if n == plotn(ploti)
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
               plot(xdata(1:n,:),zeros(n,1),'k.')
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
OutObj.AppxNorm = OutObj.AppxNorm(1:n);
OutObj.NeccFlag = OutObj.NeccFlag(2:n+1);
OutObj.xdata = xdata(1:n,:);
OutObj.fdata = fdata(1:n);


if obj.isDiagnose
   multiAppxDiagFinishPlotTable ...
      (h,legendLabel,OutObj,coli,n,ntol,nNeed,obj,xeval);
end
