function [h,ploti,legendLabel] =  ...
   multiAppxDiagPrelim(plotn,ntol,xeval,feval,obj)
   d = obj.dim;
   ploti = 2;
   h(ntol+length(plotn)+1,1) = 0; 
   legendLabel = cell(ntol+length(plotn)+1,1); 
   if d == 1
      figure
      h(1) = plot(xeval,feval,'color',obj.colorScheme(1,:));
      legendLabel{1} = '\(f(x)\)';
      axis([obj.xLim', obj.yLim'])
      hold on
   elseif d == 2
      tmp = sqrt(size(xeval,1));
      xx = reshape(xeval(:,1),tmp,tmp);
      yy = reshape(xeval(:,2),tmp,tmp);
      zz = reshape(feval,tmp,tmp);
      RemovePlotAxes
      h(1) = surf(xx,yy,zz,'FaceColor','Interp','EdgeColor','None');
      legendLabel{1} = '\(f(x)\)';
      colorbar
      hold on
   end
end