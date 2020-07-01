%% Example of easy function

gail.InitializeWorkspaceDisplay
format short e
warning('off')

[~,~,xeval,neval,Ainf,B0] = StdParam;
abstolVec = [0.05 0.02 0.01 0.005 0.002 0.001]';
ntol = size(abstolVec,1);
theta = 1;

f = @simpleFun;
feval = f(xeval);
colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; MATLABCyan; MATLABMaroon];
nmax = 500;
xdata(nmax,1) = 0;
fdata(nmax,1) = 0;
errKNull = 1;

figure %simple function
plot(xeval,feval);
xlabel('\(x\)')
ylabel('\(f(x)\)');
print('-depsc','SimpleFunPlot.eps')

%% Algorithm 1 Sample size is adaptive
kernel = @(t,x) MaternKernel(t,x,theta);
plotn = [0 1 3 nmax];
ploti = 2;
h(ntol+length(plotn)+1,1) = 0; 
legendLabel = cell(ntol+length(plotn)+1,1);
figure
h(1) = plot(xeval,feval,'color',colorScheme(1,:));
legendLabel{1} = '\(f(x)\)';
hold on
itol = 1;
abstol = abstolVec(itol);
nNeed(ntol,1) = 0;
ErrBdVec(nmax,1) = 0;
trueErr(nmax,1) = 0;
InErrBars(nmax,1) = 0;
nstart = 0;
coli = ploti;
AXvec(nmax,1) = 0;
for n = 1:nmax
   xdata(n) = seqFixedDes(n);
   fdata(n) = f(xdata(n));
   [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:n,:), xeval, kernel);
   [errKXx, errKX] = powerfun(Kmat, Kdateval, Kdiageval);
   AX = ABfun(errKX,errKNull,Ainf,B0);
   AXvec(n) = AX;
   [Appx, fluctNorm, ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX );
   ErrBdVec(n) = ErrBd;
   trueErr(n) = max(abs(feval - Appx));
   errFudge = eps*cond(Kmat);
   InErrBars(n) = sum(abs(feval - Appx) <= ErrBdx + errFudge)/neval;
   if n == plotn(ploti)
      h(coli) = plot(xeval,Appx,'color',colorScheme(mod(coli-1,6)+1,:));
      legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
      nrange = nstart+1:n;
      plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(coli-1,6)+1,:))
      ploti = ploti+1;
      coli = coli+1;
      nstart = n;
   end
   if ErrBd < abstol
      if abstol >= 0.01
         h(coli) = plot(xeval,Appx,'color',colorScheme(mod(coli-1,6)+1,:));
         legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
         nrange = nstart+1:n;
         plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(coli-1,6)+1,:))
         nstart = n;
         coli = coli+1;
      end
      nNeed(itol) = n;
      itol = itol + 1;
      if itol > ntol, break, end
      abstol = abstolVec(itol);
   end
end
disp('Algorithm 1')
% disp(['True error = ' num2str(trueErr,3) ', Error bound = ' num2str(ErrBd,3) ...
%    ', Using ' int2str(n) ' data'])
% disp(' ')
xlabel('\(x\)')
ylabel('\(f(x), \ f(x_i), \ \)APP\((\mathsf{X},\textbf{\textit{y}})\)')
legend(h(1:coli-1),legendLabel{1:coli-1},'location','south','orientation','vertical','box','off')
print('-depsc','SimpleFunAlg1.eps')

Alg1SummaryData = [(1:n)' ErrBdVec(1:n) trueErr(1:n) InErrBars(1:n)];
whEBfails = find(Alg1SummaryData(:,2) < Alg1SummaryData(:,3));
disp(['Error bound fails ' int2str(length(whEBfails)) ' times'])
disp('    for n = ')
disp(whEBfails)
whEBPointfails = find(Alg1SummaryData(:,4) < 1);
disp(['Pointwise error bound fails ' int2str(length(whEBPointfails)) ' times'])
disp('    for n = ')
disp(whEBPointfails)
fprintf(1,'\n\n\n')

fid = fopen('Example1Table.txt','w+');
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
fprintf(fid,'\\abs{f(x) - \\APP(\\mX,\\by)(x)} \\le \\errBd(\\mX,\\by)(x) & ');
fprintf(fid,'%3.0f\\%% & ',100*InErrBars(nNeed(1:ntol-1)));
fprintf(fid,'%3.0f\\%% \n', 100*InErrBars(nNeed(ntol)));
fprintf(fid,'\\end{array} \n \\]');

%% Algorithm 2 Sample location is adaptive
plotn = [0 1 3 nmax];
ploti = 2;
h(ntol+length(plotn)+1,1) = 0; 
legendLabel = cell(ntol+length(plotn)+1,1);
figure
h(1) = plot(xeval,feval,'color',colorScheme(1,:));
legendLabel{1} = '\(f(x)\)';
hold on
itol = 1;
abstol = abstolVec(itol);
nNeed(ntol,1) = 0;
ErrBdVec(nmax,1) = 0;
trueErr(nmax,1) = 0;
InErrBars(nmax,1) = 0;
nstart = 0;
coli = ploti;
AXvec(nmax,1) = 0;
BXvec(nmax,1) = 0;
for n = 1:nmax
   if n == 1
      xdata(1) = seqFixedDes(1);
   else
      xdata(n) = xeval(whKX);
   end
   fdata(n) = f(xdata(n));
   [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:n,:), xeval, kernel);
   [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   [AX, BX] = ABfun(errKX,errKNull,Ainf,B0);
   AXvec(n) = AX;
   BXvec(n) = BX;
   [Appx, fluctNorm, ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX);
   ErrBdVec(n) = ErrBd;
   trueErr(n) = max(abs(feval - Appx));
   errFudge = eps*cond(Kmat);
   InErrBars(n) = sum(abs(feval - Appx) <= ErrBdx + errFudge)/neval;
   if n == plotn(ploti)
      h(coli) = plot(xeval,Appx,'color',colorScheme(mod(coli-1,6)+1,:));
      legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
      nrange = nstart+1:n;
      plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(coli-1,6)+1,:))
      ploti = ploti+1;
      coli = coli+1;
      nstart = n;
   end
   if ErrBd < abstol
      if abstol >= 0.01
         h(coli) = plot(xeval,Appx,'color',colorScheme(mod(coli-1,6)+1,:));
         legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
         nrange = nstart+1:n;
         plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(coli-1,6)+1,:))
         nstart = n;
         coli = coli+1;
      end
      nNeed(itol) = n;
      itol = itol + 1;
      if itol > ntol, break, end
      abstol = abstolVec(itol);
   end
end
disp('Algorithm 2')
xlabel('\(x\)')
ylabel('\(f(x), \ f(x_i), \ \)APP\((\mathsf{X},\textbf{\textit{y}})\)')
legend(h(1:coli-1),legendLabel{1:coli-1},'location','south','orientation','vertical','box','off')
print('-depsc','SimpleFunAlg2.eps')

Alg2SummaryData = [(1:n)' ErrBdVec(1:n) trueErr(1:n) InErrBars(1:n)];
whEBfails = find(Alg2SummaryData(:,2) < Alg2SummaryData(:,3));
disp(['Error bound fails ' int2str(length(whEBfails)) ' times'])
disp('    for n = ')
disp(whEBfails)
whEBPointfails = find(Alg2SummaryData(:,4) < 1);
disp(['Pointwise error bound fails ' int2str(length(whEBPointfails)) ' times'])
disp('    for n = ')
disp(whEBPointfails)
fprintf(1,'\n\n\n')

fid = fopen('Example2Table.txt','w+');
fprintf(fid,'\\[ \n \\begin{array}{rccccccc} \n');
fprintf(fid,'n & ');
fprintf(fid,'%3.0f & ',nNeed(1:ntol-1));
fprintf(fid,'%3.0f \\\\ \\hline \n', nNeed(ntol));
fprintf(fid,'\\varepsilon & ');
fprintf(fid,cleanStringFJH(sprintf('%3.1E & ',abstolVec(1:ntol-1))));
fprintf(fid,[cleanStringFJH(sprintf('%3.1E', abstolVec(ntol))) ' \\\\ \\hline \n']);
fprintf(fid,'\\errBd & ');
fprintf(fid,cleanStringFJH(sprintf('%3.1E & ',ErrBdVec(nNeed(1:ntol-1)))));
fprintf(fid,[cleanStringFJH(sprintf('%3.1E', ErrBdVec(nNeed(ntol)))) '\\\\ \\hline \n']);
fprintf(fid,'\\norm[\\infty]{f - \\APP(\\mX,\\by)} & ');
fprintf(fid,cleanStringFJH(sprintf('%3.1E & ',trueErr(nNeed(1:ntol-1)))));
fprintf(fid,[cleanStringFJH(sprintf('%3.1E',trueErr(nNeed(ntol)) )) '\\\\ \\hline \n']);
fprintf(fid,'\\abs{f(x) - \\APP(\\mX,\\by)(x)} \\le \\errBd(\\mX,\\by)(x) & ');
fprintf(fid,'%3.0f\\%% & ',100*InErrBars(nNeed(1:ntol-1)));
fprintf(fid,'%3.0f\\%% \n', 100*InErrBars(nNeed(ntol)));
fprintf(fid,'\\end{array} \n \\]');

%% Algorithm 3 Sample location and kernel are adaptive
n0 = 3;
plotn = [0 n0 nmax];
ploti = 2;
h(ntol+length(plotn)+1,1) = 0; 
legendLabel = cell(ntol+length(plotn)+1,1);
figure
h(1) = plot(xeval,feval,'color',colorScheme(1,:));
legendLabel{1} = '\(f(x)\)';
hold on
itol = 1;
abstol = abstolVec(itol);
nNeed(ntol,1) = 0;
ErrBdVec(nmax,1) = 0;
trueErr(nmax,1) = 0;
InErrBars(nmax,1) = 0;
nstart = 0;
coli = ploti;
AXvec(nmax,1) = 0;
BXvec(nmax,1) = 0;
thOptimVec(nmax,1) = 0;
kernelth = @(t,x,lnth) MaternKernel(t,x,exp(lnth));
thetaRange = (-5:0.5:5)';
for n = n0:nmax
   if n == n0
      xdata(1:n0) = seqFixedDes(1:n0);
      fdata(1:n0) = f(xdata(1:n0));
   else
      xdata(n) = xeval(whKX);
      fdata(n) = f(xdata(n));
   end
   lnthOptim = selectTheta(thetaRange,kernelth,xdata(1:n),fdata(1:n), ...
      xeval,errKNull,Ainf,B0);
   thetaOptim = exp(lnthOptim);
   thOptimVec(n) = thetaOptim;
   kernel = @(t,x) MaternKernel(t,x,thetaOptim);
   [Kmat, Kdateval, Kdiageval] = KMP(xdata(1:n,:), xeval, kernel);
   [errKXx, errKX, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   [AX, BX] = ABfun(errKX,errKNull,Ainf,B0);
   AXvec(n) = AX;
   BXvec(n) = BX;   
   [Appx, fluctNorm, ErrBdx, ErrBd] = Approx(fdata(1:n), Kmat, Kdateval, errKXx, errKX, AX );
   ErrBdVec(n) = ErrBd;
   trueErr(n) = max(abs(feval - Appx));
   errFudge = eps*cond(Kmat);
   InErrBars(n) = sum(abs(feval - Appx) <= ErrBdx + errFudge)/neval;
   %if InErrBars(n) < 1, keyboard, end
   if n == plotn(ploti)
      h(coli) = plot(xeval,Appx,'color',colorScheme(mod(coli,6)+1,:));
      legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
      nrange = nstart+1:n;
      plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(coli,6)+1,:))
      ploti = ploti+1;
      coli = coli+1;
      nstart = n;
   end
   if ErrBd < abstol
      if abstol >= 0.01
         h(coli) = plot(xeval,Appx,'color',colorScheme(mod(coli,6)+1,:));
         legendLabel{coli} = ['\(n = ' int2str(n) '\)'];
         nrange = nstart+1:n;
         plot(xdata(nrange),fdata(nrange),'.','color',colorScheme(mod(coli,6)+1,:))
         nstart = n;
         coli = coli+1;
      end
      nNeed(itol) = n;
      itol = itol + 1;
      if itol > ntol, break, end
      abstol = abstolVec(itol);
   end
end
disp('Algorithm 3')
xlabel('\(x\)')
ylabel('\(f(x), \ f(x_i), \ \)APP\((\mathsf{X},\textbf{\textit{y}})\)')
legend(h(1:coli-1),legendLabel{1:coli-1},'location','south','orientation','vertical','box','off')
print('-depsc','SimpleFunAlg3.eps')

Alg3SummaryData = [(n0:n)' ErrBdVec(n0:n) trueErr(n0:n) InErrBars(n0:n)];
whEBfails = find(Alg3SummaryData(:,2) < Alg3SummaryData(:,3));
disp(['Error bound fails ' int2str(length(whEBfails)) ' times'])
disp('    for n = ')
disp(whEBfails)
whEBPointfails = find(Alg3SummaryData(:,4) < 1);
disp(['Pointwise error bound fails ' int2str(length(whEBPointfails)) ' times'])
disp('    for n = ')
disp(whEBPointfails)

fid = fopen('Example3Table.txt','w+');
fprintf(fid,'\\[ \n \\begin{array}{rccccccc} \n');
fprintf(fid,'n & ');
fprintf(fid,'%3.0f & ',nNeed(1:ntol-1));
fprintf(fid,'%3.0f \\\\ \\hline \n', nNeed(ntol));
fprintf(fid,'\\varepsilon & ');
fprintf(fid,cleanStringFJH(sprintf('%3.1E & ',abstolVec(1:ntol-1))));
fprintf(fid,[cleanStringFJH(sprintf('%3.1E', abstolVec(ntol))) ' \\\\ \\hline \n']);
fprintf(fid,'\\errBd & ');
fprintf(fid,cleanStringFJH(sprintf('%3.1E & ',ErrBdVec(nNeed(1:ntol-1)))));
fprintf(fid,[cleanStringFJH(sprintf('%3.1E', ErrBdVec(nNeed(ntol)))) '\\\\ \\hline \n']);
fprintf(fid,'\\norm[\\infty]{f - \\APP(\\mX,\\by)} & ');
fprintf(fid,cleanStringFJH(sprintf('%3.1E & ',trueErr(nNeed(1:ntol-1)))));
fprintf(fid,[cleanStringFJH(sprintf('%3.1E',trueErr(nNeed(ntol)) )) '\\\\ \\hline \n']);
fprintf(fid,'\\abs{f(x) - \\APP(\\mX,\\by)(x)} \\le \\errBd(\\mX,\\by)(x) & ');
fprintf(fid,'%3.0f\\%% & ',100*InErrBars(nNeed(1:ntol-1)));
fprintf(fid,'%3.0f\\%% \n', 100*InErrBars(nNeed(ntol)));
fprintf(fid,'\\end{array} \n \\]');

   
   