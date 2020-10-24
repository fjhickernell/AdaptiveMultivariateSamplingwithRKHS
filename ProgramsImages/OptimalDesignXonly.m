%% Construct Optimal Design With Only x Values
function [xdata, prm] = OptimalDesignXonly(kernel, xeval, prm)
d = size(xeval,2);
xdata(prm.nmax,d) = 0;
plotn = [prm.plotn prm.nmax];
nstart = 0;
if prm.isDiagnose
   [h,ploti,legendLabel] =  ...
      multiDesignDiagPrelim(plotn,xeval,prm);
   coli = ploti;
end

for n = 1:prm.nmax
   gail.print_iterations(n,'n',true)
   if n == 1
      xdata(1,:) = 0.5;
      errKXx = [];
      errKX = [];
   else
      [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:n-1,:), xeval, kernel);
      [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
      xdata(n,:) = xeval(whKX,:);
   end
   if prm.isDiagnose
      if n == plotn(ploti)
         [h,legendLabel,coli,nstart] =  ...
            multiDesignDiagAddData(h,legendLabel,coli,n,nstart, ...
            xdata,xeval,errKXx,errKX,prm);
         ploti = ploti+1;
      end
   end
end
fprintf('\n')

if prm.isDiagnose
   multiDesignDiagFinishPlotTable ...
      (h,legendLabel,coli,prm,xeval);
end
