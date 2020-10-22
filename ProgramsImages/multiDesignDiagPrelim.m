function [h,ploti,legendLabel] =  ...
   multiDesignDiagPrelim(plotn,xeval,prm)
   d = size(xeval,2);
   ploti = 1;
   h(length(plotn)+1,1) = 0; 
   legendLabel = cell(length(plotn)+1,1); 
   if d == 1
      figure
      axis([prm.xLim', -0.1 0.1])
   else
      gail.RemovePlotAxes
   end
    hold on
end