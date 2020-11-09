function OutObject = RunFunAppxExample(object)

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
   
return

dim == 2
    
    xeval = 0:0.005:1;
    obj = obj(1);
    [xx, yy]= meshgrid(xeval);
    xeval  = [xx(:),yy(:)];
    feval = f(xeval);
    fplot = reshape(feval,size(xx));
    surf(xx,yy,fplot,'FaceColor','Interp','EdgeColor','None')
    colorbar
    xlabel('\(x_1\)')
    ylabel('\(x_2\)');
    zlabel('\(f(x_1,x_2)\)')
    %title('\(f(\textbf{\textit{y}}) = \cos(x_1 + x_2) \exp(x_1 * x_2)\)')
    print('-depsc',[obj.fname 'Plot.eps'])
    
    for kk = 1:size(obj,2)
        obj = obj(kk);
        %kernel = @(t,x) kernelth(t,x,param(kk).theta);
        disp([obj.fname ' ' obj.kername ' ' obj.whDes ' ' obj.whObj ' ' obj.AlgName])
        [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag,err,xdata,fdata,obj] = ...
                AdaptAlgo3(f, kernelth, xeval, feval, obj);
        disp(['Necessary condition flag = ' int2str(NeccFlag(end))])
        fprintf(1,'\n\n')
        delta = 0.1;
        fplot = reshape(err,size(xx));
        gail.RemovePlotAxes
        surf(xx,yy,fplot,'FaceColor','Interp','EdgeColor','None')
        colorbar
        hold on
        plot3(xdata(:,1),xdata(:,2),fdata + 10*delta,'.','color','magenta');
        %title('\(f(\textbf{\textit{y}}) = \cos(x_1 + x_2) \exp(x_1 * x_2)\)')
        print('-depsc',[obj.fname 'ErrorPlot.eps'])
    end

    n = 2000;
    xeval = rand(n,dim);
    feval = f(xeval);
    for kk = 1:size(obj,2)
        obj = obj(kk);
        %kernel = @(t,x) kernelth(t,x,param(kk).theta);
        disp([obj.fname ' ' obj.kername ' ' obj.whDes ' ' obj.whObj ' ' obj.AlgName])
        [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag,err,xdata,fdata,obj] = ...
                AdaptAlgo3(f, kernelth, xeval, feval, obj);
        %disp(['Necessary condition flag = ' int2str(NeccFlag(end))])
        fprintf(1,'\n\n')
    end
end
   