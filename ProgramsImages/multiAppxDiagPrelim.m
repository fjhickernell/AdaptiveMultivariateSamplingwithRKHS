function [h,ploti,legendLabel] =  ...
   multiAppxDiagPrelim(plotn,ntol,xeval,feval,prm)

   ploti = 2;
   h(ntol+length(plotn)+1,1) = 0; 
   legendLabel = cell(ntol+length(plotn)+1,1);
   figure
   h(1) = plot(xeval,feval,'color',prm.colorScheme(1,:));
   legendLabel{1} = '\(f(x)\)';
   hold on
end
