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
dth = size(obj.thetaRange,2);
thOptimVec(obj.nmax,dth) = 0;
for n = obj.n0:obj.nmax
   gail.print_iterations(n,'n',true)
   if n == obj.n0
      if strcmp(obj.whDes,'unif_grid') && n > 1 && d == 1
         xdata(1:n,:) = (0:n-1)'/(n-1);
      elseif strcmp(obj.whDes,'unifChebyshev') && n > 1 %Chebyshev
         xdata(1:n,:) = (0:n-1)'/(n-1);
         xdata(1:n,:) = (1+sin(pi*(-1/2 + xdata(1:n,:))))/2;
      elseif strcmp(obj.whDes,'seqChebyshev') %sequential Chebyshev
         xdata(1:n,:) = seqFixedDes(1:n,d,1/3,'Chebyshev');
      elseif strcmp(obj.whDes,'adapt_th') && n > 1
         kernel = @(t,x) obj.kernelth(t,x,obj.theta);
         xdata(1,:) = seqFixedDes(1,d);
          for ii = 2:obj.n0
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
   [thOptim,obj] = selectTheta(obj,xdata(1:n,:),fdata(1:n), ...
      xeval);
   thOptimVec(n,:) = thOptim;
   kernel = @(t,x) obj.kernelth(t,x,obj.currentTheta);
   [Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata(1:n,:), xeval, kernel);
   [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   [AX, BX] = ABfun(errKX,errKNull,obj.Ainf,obj.B0);
   AXvec(n) = AX;
   BXvec(n) = BX;
   [Appx, OutObj.AppxNorm(n), ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX );
   OutObj.ErrBdVec(n) = ErrBd;
   err = abs(feval-Appx);
   OutObj.trueErr(n) = max(err);
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
OutObj.ErrBdVec = OutObj.ErrBdVec(1:n);
OutObj.trueErr = OutObj.trueErr(1:n);
OutObj.InErrBars = OutObj.InErrBars(1:n);
OutObj.xdata = xdata(1:n,:);
OutObj.fdata = fdata(1:n);
obj.final_theta = thOptim;

if obj.isDiagnose
   multiAppxDiagFinishPlotTable ...
      (h,legendLabel,OutObj,coli,n,ntol,nNeed,obj,xeval);
end
