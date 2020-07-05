function AlgSummaryData =  multiAppxDiagFinishPlotTable ...
   (h,legendLabel,abstolVec,ErrBdVec,trueErr,InErrBars, ...
   coli,n,ntol,nNeed,prm,algname)

   AlgSummaryData = [(1:n)' ErrBdVec trueErr InErrBars];
   xlabel('\(x\)')
   ylabel('\(f(x), \ f(x_i), \ \)APP\((\mathsf{X},\textbf{\textit{y}})\)')
   legend(h(1:coli-1),legendLabel{1:coli-1},'location',prm.legendPos,'orientation','vertical','box','off')
   print('-depsc',[algname '_' prm.fname '_' prm.kername '.eps'])
   whEBfails = find(AlgSummaryData(prm.n0:n,2) < AlgSummaryData(prm.n0:n,3));
   disp(['Error bound fails ' int2str(length(whEBfails)) ' times'])
   disp('    for n = ')
   disp(whEBfails+prm.n0-1)
   whEBPointfails = find(AlgSummaryData(prm.n0:n,4) < 1);
   disp(['Pointwise error bound fails ' int2str(length(whEBPointfails)) ' times'])
   disp('    for n = ')
   disp(whEBPointfails+prm.n0-1)
   fprintf(1,'\n')

   fid = fopen([algname '_Out' prm.fname '_' prm.kername '.txt'],'w+');
   fprintf(fid,'\\[ \n \\begin{array}{rccccccc} \n');
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
   fprintf(fid,'\\end{array} \n \\]');end
