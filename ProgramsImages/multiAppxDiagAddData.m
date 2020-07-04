function [h,legendLabel,coli,nstart] =  ...
   multiAppxDiagAddData(h,legendLabel,coli,n,nstart,xdata,fdata,xeval,Appx,colorScheme)

   h(coli) = plot(xeval,Appx,'color',colorScheme(mod(coli-1,6)+1,:));
   legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
   nrange = nstart+1:n;
   plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(coli-1,6)+1,:))
   coli = coli+1;
   nstart = n;
end
