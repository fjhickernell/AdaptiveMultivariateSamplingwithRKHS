function [h,legendLabel,coli,nstart] = multiDesignDiagAddData ...
   (h,legendLabel,coli,n,nstart,xdata,prm)
   [~,d] = size(xdata);
   if d==1
       legendLabel{coli} = ['\(n = ' int2str(n) '\)']; 
       nrange = nstart+1:n;
       h(coli) = plot(xdata(nrange),zeros(n-nstart,1),'.','color',prm.colorScheme(mod(coli-1,6)+1,:));
       coli = coli+1;
       nstart = n;
   else
       nrange = nstart+1:n;
       h(coli) = plot(xdata(nrange,1),xdata(nrange,2),'.','color',prm.colorScheme(mod(coli-1,6)+1,:));
       legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
       coli = coli+1;
       nstart = n;
   end
end
