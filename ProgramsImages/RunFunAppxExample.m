function [OutObject,object] = RunFunAppxExample(object)

warning('off','all')
nObj = length(object);
OutObject = FunAppxOut(object);
for kk = 1:nObj
   obj = object(kk);
   xeval = obj.xeval;
   f = obj.f;
   feval = f(xeval);
   d = obj.dim;
   if d == 1
      figure %plot function
      plot(xeval,feval);
      axis([obj.xLim', obj.yLim'])
      xlabel('\(x\)')
      ylabel('\(f(x)\)');
      print('-depsc',[obj.fname 'Plot.eps'])
   elseif d == 2
      figure
      neval = size(feval,1);
      xx = reshape(xeval(:,1), sqrt([neval,neval]));
      yy = reshape(xeval(:,2), sqrt([neval,neval]));
      fplot = reshape(feval,size(xx));
      surf(xx,yy,fplot,'FaceColor','Interp','EdgeColor','None')
      colorbar
      xlabel('\(x_1\)')
      ylabel('\(x_2\)');
      zlabel('\(f(x_1,x_2)\)')
      print('-depsc',[obj.fname 'Plot.eps'])
   end

   kernel = @(t,x) obj.kernelth(t,x,obj.theta);
   disp([obj.fname ' ' obj.kername ' ' obj.whDes ' ' obj.whObj ' ' obj.algoname])
   OutObject(kk) = obj.Algo(f, kernel, xeval, feval, obj);
   disp(['Necessary condition flag = ' int2str(OutObject(kk).NeccFlag(end))])
   fprintf(1,'\n\n')
end
   
