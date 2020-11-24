function AlgSummaryData =  multiAppxDiagFinishPlotTable ...
   (h,legendLabel,OutObj,coli,n,ntol,nNeed,obj,xeval,xdata,trueErrX,ErrBdx)

   AlgSummaryData = [(1:n)' OutObj.ErrBdVec OutObj.trueErr OutObj.InErrBars];
   if obj.dim <=2
   if obj.dim == 1
       xlabel('\(x\)')
       ylabel('\(f(x), \ f(x_i), \ \)APP\((\mathsf{X},\textbf{\textit{y}})\)')
       legend(h(1:coli-1),legendLabel{1:coli-1},'location',obj.legendPos,'orientation','vertical','box','off')
   elseif obj.dim == 2
       xlabel('\(x_1\)')
       ylabel('\(x_2\)')
       legend(h(1:coli-1),legendLabel{1:coli-1},'location',obj.legendPos,'orientation','vertical','box','off')
 
       tmp = sqrt(size(xeval,1));
       xx = reshape(xeval(:,1),tmp,tmp);
       yy = reshape(xeval(:,2),tmp,tmp);
       zz = reshape(trueErrX,tmp,tmp);
       RemovePlotAxes
       surf(xx,yy,zz,'FaceColor','Interp','EdgeColor','None');
       delta = 0.1;
       plot3(xdata(1:n,1),xdata(1:n,2),trueErrX(1:n) + delta,'.','color',obj.colorScheme(1,:));
       xlabel('\(x_1\)')
       ylabel('\(x_2\)')
       colorbar
       legend('True Error','location',obj.legendPos,'orientation','vertical','box','off')
       
       RemovePlotAxes
       zz = reshape(ErrBdx,tmp,tmp);
       surf(xx,yy,zz,'FaceColor','Interp','EdgeColor','None');
       delta = 0.1;
       plot3(xdata(1:n,1),xdata(1:n,2),ErrBdx(1:n) + delta,'.','color',obj.colorScheme(1,:));
       xlabel('\(x_1\)')
       ylabel('\(x_2\)')
       colorbar
       legend('Error Bound','location',obj.legendPos,'orientation','vertical','box','off')
   end
   print('-depsc',[obj.algoname '_' obj.fname '_' obj.kername '_' obj.whDes '_' ...
      obj.whObj '_theta_' num2str(obj.theta(1)) '.eps'])
   end
   whEBfails = find(AlgSummaryData(obj.n0:n,2) < AlgSummaryData(obj.n0:n,3));
   disp(['Error bound fails ' int2str(length(whEBfails)) ' times'])
   disp('    for n = ')
   disp(whEBfails+obj.n0-1)
   whEBPointfails = find(AlgSummaryData(obj.n0:n,4) < 1);
   disp(['Pointwise error bound fails ' int2str(length(whEBPointfails)) ' times'])
   disp('    for n = ')
   disp(whEBPointfails+obj.n0-1)
   fprintf(1,'\n')

   fid = fopen([obj.algoname '_Out_' obj.fname '_' obj.kername '_' obj.whDes '_' ...
       obj.whObj '_theta_' num2str(obj.theta(1)) '.txt'],'w+');
   fprintf(fid,'\\[ \n \\begin{array}{rccccccc} \n');
   fprintf(fid,['\\multicolumn{' int2str(ntol+1) '}{l}{A_\\infty = ' cleanStringFJH(sprintf('%2.2f',obj.Ainf))]); 
   fprintf(fid,['\\qquad B_0 = ' cleanStringFJH(sprintf('%2.2g',obj.B0))]);
   fprintf(fid,['\\qquad \\mT = (' cleanStringFJH(sprintf('%2.2g, %2.2g,',xeval(1:2)))]);
   fprintf(fid,['\\ldots, ' cleanStringFJH(sprintf('%2.2g)',xeval(end)))]);
   fprintf(fid,['\\qquad n_0 = ' cleanStringFJH(sprintf('%2.0f',obj.n0))]);
   fprintf(fid,'} \\\\ \\hline \n');
   if any(strcmp(obj.algoname,'Alg3'))
      fprintf(fid,['\\multicolumn{' int2str(ntol+1) '}{l}{\\text{Final } \\theta = ']);
      dth = length(obj.final_theta);
      if dth == 1
         fprintf(fid,cleanStringFJH(sprintf('%4.1E',obj.final_theta)));
      else
         fprintf(fid,'(');
         for jj = 1:dth-1
            fprintf(fid,cleanStringFJH(sprintf('%4.1E, ',obj.final_theta(jj))));
         end
         fprintf(fid,cleanStringFJH(sprintf('%4.1E)',obj.final_theta(dth))));
      end
      fprintf(fid,'} \\\\ \\hline \n');
   end   
   [~,ntolReach] = find(nNeed <= obj.nmax,1,'last');
   fprintf(fid,'n & ');
   fprintf(fid,'%3.0f & ',nNeed(1:ntolReach-1));
   fprintf(fid,'%3.0f \\\\ \\hline \n', nNeed(ntolReach));
   fprintf(fid,'\\varepsilon & ');
   fprintf(fid,cleanStringFJH(sprintf('%3.1E & ',obj.abstolVec(1:ntolReach-1))));
   fprintf(fid,[cleanStringFJH(sprintf('%3.1E', obj.abstolVec(ntolReach))) ' \\\\ \\hline \n']);
   fprintf(fid,'\\errBd & ');
   fprintf(fid,cleanStringFJH(sprintf('%3.1E & ',OutObj.ErrBdVec(nNeed(1:ntolReach-1)))));
   fprintf(fid,[cleanStringFJH(sprintf('%3.1E', OutObj.ErrBdVec(nNeed(ntolReach)))) ' \\\\ \\hline \n']);
   fprintf(fid,'\\norm[\\infty]{f - \\APP(\\mX,\\by)} & ');
   fprintf(fid,cleanStringFJH(sprintf('%3.1E & ',OutObj.trueErr(nNeed(1:ntolReach-1)))));
   fprintf(fid,[cleanStringFJH(sprintf('%3.1E',OutObj.trueErr(nNeed(ntol)) )) ' \\\\ \\hline \n']);
   fprintf(fid,'\\text{Ptwise Error Bound} & ');
   fprintf(fid,'%3.0f\\%% & ',100*OutObj.InErrBars(nNeed(1:ntolReach-1)));
   fprintf(fid,'%3.0f\\%% \n', 100*OutObj.InErrBars(nNeed(ntolReach)));
   fprintf(fid,'\\end{array} \n \\]');
end
