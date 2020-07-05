function [h,legendLabel,coli,nstart] = multiAppxDiagAddData ...
   (h,legendLabel,coli,n,nstart,xdata,fdata,xeval,Appx,colorScheme,abstol)

   h(coli) = plot(xeval,Appx,'color',colorScheme(mod(coli-1,6)+1,:));
   if isfinite(abstol)
      legendLabel{coli} = ['\(n = ' int2str(n) ',\ \varepsilon = ' num2str(abstol) '\)'];
   else
      legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
   end
   nrange = nstart+1:n;
   plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(coli-1,6)+1,:))
   coli = coli+1;
   nstart = n;
end
