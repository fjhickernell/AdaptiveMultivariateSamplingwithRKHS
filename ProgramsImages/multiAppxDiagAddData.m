function [h,legendLabel,coli,nstart] = multiAppxDiagAddData ...
   (h,legendLabel,coli,n,nstart,xdata,fdata,xeval,Appx,obj,abstol,trueErrX,fignum)
   d = obj.dim;
   if d == 1 
       h(coli) = plot(xeval,Appx,'color',obj.colorScheme(mod(coli-1,6)+1,:));
       if isfinite(abstol)
           legendLabel{coli} = ['\(n = ' int2str(n) ',\ \varepsilon = ' num2str(abstol) '\)'];
       else
           legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
       end 
       nrange = nstart+1:n;
       plot(xdata(nrange),fdata(nrange),'.','color',obj.colorScheme(mod(coli-1,6)+1,:))
       coli = coli+1;
       nstart = n;
   elseif d == 2
       delta = 0.1;
       nrange = nstart+1:n;
       h(coli) = plot3(xdata(nrange,1),xdata(nrange,2),fdata(nrange) + delta,'.','color',obj.colorScheme(mod(coli-1,6)+1,:));
       if isfinite(abstol)
           legendLabel{coli} = ['\(n = ' int2str(n) ',\ \varepsilon = ' num2str(abstol) '\)'];     
       else
           legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
       end
       coli = coli+1;
       nstart = n;
   end
end
