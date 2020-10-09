function RunExample(f,param,kernelth, dim)

if nargin < 4
    [~,~,xeval] = StdParam;
    
    feval = f(xeval);
    prm = param(1);
    figure %plot function
    plot(xeval,feval);
    axis([prm.xLim', prm.yLim'])
    xlabel('\(x\)')
    ylabel('\(f(x)\)');
    print('-depsc',[prm.fname 'Plot.eps'])
    
    for kk = 1:size(param,2)
        prm = param(kk);
        kernel = @(t,x) kernelth(t,x,param(kk).theta);
        disp([prm.fname ' ' prm.kername ' ' prm.whDes ' ' prm.whObj ' ' prm.AlgName])
        if strcmp(prm.AlgName,'Algo1')
            [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
                AdaptAlgo1(f, kernel, xeval, feval, prm);
        elseif strcmp(prm.AlgName,'Algo2')
            [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
                AdaptAlgo2(f, kernel, xeval, feval, prm);
        elseif strcmp(prm.AlgName,'Algo3')
            [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag] = ...
                AdaptAlgo3(f, kernelth, xeval, feval, prm);
        end
        disp(['Necessary condition flag = ' int2str(NeccFlag(end))])
        fprintf(1,'\n\n')
    end
    
elseif dim == 2
    
    xeval = 0:0.005:1;
    prm = param(1);
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
    print('-depsc',[prm.fname 'Plot.eps'])
    
    for kk = 1:size(param,2)
        prm = param(kk);
        %kernel = @(t,x) kernelth(t,x,param(kk).theta);
        disp([prm.fname ' ' prm.kername ' ' prm.whDes ' ' prm.whObj ' ' prm.AlgName])
        [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag,err,xdata,fdata,prm] = ...
                AdaptAlgo3(f, kernelth, xeval, feval, prm);
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
        print('-',[prm.fname 'ErrorPlot.eps'])
    end
else
    n = 2000;
    xeval = rand(n,dim);
    feval = f(xeval);
    for kk = 1:size(param,2)
        prm = param(kk);
        %kernel = @(t,x) kernelth(t,x,param(kk).theta);
        disp([prm.fname ' ' prm.kername ' ' prm.whDes ' ' prm.whObj ' ' prm.AlgName])
        [Appx, ErrBdx, ErrBdVec, trueErr, InErrBars, AppxNorm, NeccFlag,err,xdata,fdata,prm] = ...
                AdaptAlgo3(f, kernelth, xeval, feval, prm);
        %disp(['Necessary condition flag = ' int2str(NeccFlag(end))])
        fprintf(1,'\n\n')
    end
end
   