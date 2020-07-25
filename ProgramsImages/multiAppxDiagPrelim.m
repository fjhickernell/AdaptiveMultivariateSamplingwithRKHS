function [h,ploti,legendLabel] =  ...
   multiAppxDiagPrelim(plotn,ntol,xeval,feval,prm)
   d = size(xeval,2);
   ploti = 2;
   h(ntol+length(plotn)+1,1) = 0; 
   legendLabel = cell(ntol+length(plotn)+1,1);
   figure
   h(1) = plot(xeval,feval,'color',prm.colorScheme(1,:));
   legendLabel{1} = '\(f(x)\)';
   if d ==1
      axis([prm.xLim', prm.yLim'])
   end
   hold on
end
