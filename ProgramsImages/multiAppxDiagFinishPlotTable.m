function AlgSummaryData =  multiAppxDiagFinishPlotTable ...
   (h,legendLabel,abstolVec,ErrBdVec,trueErr,InErrBars, ...
   coli,n,ntol,nNeed,prm,xeval,algname)

   AlgSummaryData = [(1:n)' ErrBdVec trueErr InErrBars];
   xlabel('\(x\)')
   ylabel('\(f(x), \ f(x_i), \ \)APP\((\mathsf{X},\textbf{\textit{y}})\)')
   legend(h(1:coli-1),legendLabel{1:coli-1},'location',prm.legendPos,'orientation','vertical','box','off')
   print('-depsc',[algname '_' prm.fname '_' prm.kername '_' prm.whDes '_' ...
      prm.whobj '_theta_' prm.theta.eps'])
   whEBfails = find(AlgSummaryData(prm.n0:n,2) < AlgSummaryData(prm.n0:n,3));
   disp(['Error bound fails ' int2str(length(whEBfails)) ' times'])
   disp('    for n = ')
   disp(whEBfails+prm.n0-1)
   whEBPointfails = find(AlgSummaryData(prm.n0:n,4) < 1);
   disp(['Pointwise error bound fails ' int2str(length(whEBPointfails)) ' times'])
   disp('    for n = ')
   disp(whEBPointfails+prm.n0-1)
   fprintf(1,'\n')

   fid = fopen([algname '_Out_' prm.fname '_' prm.kername '_' prm.whDes '_' prm.whobj '.txt'],'w+');
   fprintf(fid,'\\[ \n \\begin{array}{rccccccc} \n');
   fprintf(fid,['\\multicolumn{7}{l}{A_\\infty = ' cleanStringFJH(sprintf('%2.2f',prm.Ainf))]); 
   fprintf(fid,['\\qquad B_0 = ' cleanStringFJH(sprintf('%2.2g',prm.B0))]);
   fprintf(fid,['\\qquad \\mT = (' cleanStringFJH(sprintf('%2.2g, %2.2g,',xeval(1:2)))]);
   fprintf(fid,['\\ldots, ' cleanStringFJH(sprintf('%2.2g)',xeval(end)))]);
   if any(strcmp(algname,{'Alg2','Alg3'}))
      fprintf(fid,['\\qquad n_0 = ' cleanStringFJH(sprintf('%2.0f',prm.n0))]);
   end
   fprintf(fid,'} \\\\ \\hline \n');
   if any(strcmp(algname,'Alg3'))
      fprintf(fid,'\\multicolumn{7}{l}{\text{Final } \\theta = ');
      dth = length(prm.final_theta);
      if dth == 1
         fprintf(fid,cleanStringFJH(sprintf('%4.1E',prm.final_theta)));
      else
         for jj ' ...
         for j
         cleanStringFJH(sprintf('%4.1E',prm.final_theta))]);
   end   
   fprintf(fid,'n & ');
   fprintf(fid,'%3.0f & ',nNeed(1:ntol-1));
   fprintf(fid,'%3.0f \\\\ \\hline \n', nNeed(ntol));
   fprintf(fid,'\\varepsilon & ');
   fprintf(fid,cleanStringFJH(sprintf('%3.1E & ',abstolVec(1:ntol-1))));
   fprintf(fid,[cleanStringFJH(sprintf('%3.1E', abstolVec(ntol))) ' \\\\ \\hline \n']);
   fprintf(fid,'\\errBd & ');
   fprintf(fid,cleanStringFJH(sprintf('%3.1E & ',ErrBdVec(nNeed(1:ntol-1)))));
   fprintf(fid,[cleanStringFJH(sprintf('%3.1E', ErrBdVec(nNeed(ntol)))) ' \\\\ \\hline \n']);
   fprintf(fid,'\\norm[\\infty]{f - \\APP(\\mX,\\by)} & ');
   fprintf(fid,cleanStringFJH(sprintf('%3.1E & ',trueErr(nNeed(1:ntol-1)))));
   fprintf(fid,[cleanStringFJH(sprintf('%3.1E',trueErr(nNeed(ntol)) )) ' \\\\ \\hline \n']);
   fprintf(fid,'\\text{Ptwise Error Bound} & ');
   fprintf(fid,'%3.0f\\%% & ',100*InErrBars(nNeed(1:ntol-1)));
   fprintf(fid,'%3.0f\\%% \n', 100*InErrBars(nNeed(ntol)));
   fprintf(fid,'\\end{array} \n \\]');
end
