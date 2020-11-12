%% Algorithm 3 Sample location is adaptive and theta is adaptive
function OutObj = AdaptAlgo3(f, ~, xeval, feval, obj)
OutObj = FunAppxOut(obj);
d = obj.dim;
neval = size(xeval,1);
xdata(obj.nmax,d) = 0;
fdata(obj.nmax,1) = 0;
ntol = size(obj.abstolVec,1);
plotn = [0 obj.n0 obj.nmax];
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
OutObj.NeccFlag = NaN;
nstart = 0;
AXvec(obj.nmax,1) = 0;
BXvec(obj.nmax,1) = 0;
dth = size(obj.theta,2);
OutObj.currentTheta = obj.theta;
thOptimVec(obj.nmax,dth) = 0;
for n = obj.n0:obj.nmax
   print_iterations(n,'n',true)
   if n == obj.n0
      if strcmp(obj.whDes,'adapt_th')
         xdata(1,:) = (obj.xLim(1,:) + obj.xLim(2,:))/2;
         kernel = @(t,x) obj.kernelth(t,x,OutObj.currentTheta);
         for nn = 2:n
            [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:nn-1,:), xeval, kernel);
            [~, ~, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
            xdata(nn,:) = xeval(whKX,:);
         end
      else
         xdata = fixedDesign(1:n,obj);     
      end
      fdata(1:n) = f(xdata(1:n,:));
   else
      xdata(n,:) = xeval(whKX,:);
      fdata(n) = f(xdata(n,:));
   end
%    if n ==14, keyboard, end
   [thOptimVec(n,:),OutObj.currentTheta] = selectTheta(obj,xdata(1:n,:),fdata(1:n), ...
      xeval,OutObj.currentTheta);
   kernel = @(t,x) obj.kernelth(t,x,OutObj.currentTheta);
   [Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata(1:n,:), xeval, kernel);
   [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   [AX, BX] = ABfun(errKX,errKNull,obj.Ainf,obj.B0);
   AXvec(n) = AX;
   BXvec(n) = BX;
   [Appx, OutObj.AppxNorm(n), ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX );
   OutObj.ErrBdVec(n) = ErrBd;
   trueErrX = abs(feval-Appx);
   OutObj.trueErr(n) = max(trueErrX);
   errFudge = eps*cond(Kmat);
   OutObj.InErrBars(n) = sum(trueErrX <= ErrBdx + errFudge)/neval;
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
   if obj.isDiagnose && n > nstart
      if n >= plotn(ploti)
         [h,legendLabel,coli,nstart] =  ...
            multiAppxDiagAddData(h,legendLabel,coli,n,nstart, ...
            xdata,fdata,xeval,Appx,obj,NaN);
         ploti = ploti+1;
      end
   end
end
fprintf('\n')
OutObj.ErrBdVec = OutObj.ErrBdVec(1:n);
OutObj.trueErr = OutObj.trueErr(1:n);
OutObj.InErrBars = OutObj.InErrBars(1:n);
OutObj.xdata = xdata(1:n,:);
OutObj.fdata = fdata(1:n);
OutObj.thetaOptimalVec = thOptimVec(1:n,:);

if obj.isDiagnose
   multiAppxDiagFinishPlotTable ...
      (h,legendLabel,OutObj,coli,n,ntol,nNeed,obj,xeval,xdata,trueErrX,ErrBdx);
end
