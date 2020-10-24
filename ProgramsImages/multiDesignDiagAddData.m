function [h,legendLabel,coli,nstart] = multiDesignDiagAddData ...
   (h,legendLabel,coli,n,nstart,xdata,xeval,errKXx,errKX,prm)
   [~,d] = size(xdata);
   if d==1
       legendLabel{coli} = ['\(n = ' int2str(n) '\)']; 
       nrange = nstart+1:n;
       our_color = prm.colorScheme(mod(coli-1,6)+1,:);
       h(coli) = plot(xdata(nrange),1.2*prm.yLim(1)*ones(n-nstart,1),'.','color',our_color);
       if n > 1
          plot(xdata(n),errKX,'.','color',our_color);
          plot(xeval,errKXx,'-','color',prm.colorScheme(mod(coli-1,6)+1,:));
       end
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
