function [h,ploti,legendLabel] =  ...
   multiAppxDiagPrelim(plotn,ntol,xeval,feval,prm)
   d = size(xeval,2);
   ploti = 2;
   h(ntol+length(plotn)+1,1) = 0; 
   legendLabel = cell(ntol+length(plotn)+1,1);
   if d == 1
      figure
      h(1) = plot(xeval,feval,'color',prm.colorScheme(1,:));
      legendLabel{1} = '\(f(x)\)';
      axis([prm.xLim', prm.yLim'])
   else
      tmp = sqrt(size(xeval,1));
      xx = reshape(xeval(:,1),tmp,tmp);
      yy = reshape(xeval(:,2),tmp,tmp);
      zz = reshape(feval,tmp,tmp);
      gail.RemovePlotAxes
      h(1) = surf(xx,yy,zz,'FaceColor','Interp','EdgeColor','None');
      colorbar
   end
    hold on
end