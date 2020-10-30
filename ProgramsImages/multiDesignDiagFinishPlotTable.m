function multiDesignDiagFinishPlotTable ...
   (h,legendLabel,coli,prm,xeval)

   if size(xeval,2) < 2
       xlabel('\(x\)')
       legend(h(1:coli-1),legendLabel{1:coli-1},'location',prm.legendPos, ...
          'orientation','horizontal','box','off','NumColumns',ceil((coli-1)/2))
   else
       xlabel('\(x_1\)')
       ylabel('\(x_2\)')
       legend(h(1:coli-1),legendLabel{1:coli-1},'location',prm.legendPos,'orientation','vertical','box','off')
   end
   print('-depsc',['OptimalDesign_' prm.kername '_.eps'])
 
