%% Test ill conditoning

clearvars, close all
InitializeDisplay


%% Rank of Gram matrix with as n increases
thetavec = [0.01 0.1 1 10 100];
nth = size(thetavec,2);
mmax = 6;
nmax = 2^mmax;
nvec = 2.^(0:mmax)'+1;
xdata_all = (0:nmax)'/nmax;
xeval = 0;
rankKmat(nth,mmax+1) = 0;
rankKmatB(nth,mmax+1) = 0;
legtxt = cell(nth,1);
for ii = 1:nth
   th = thetavec(ii);
   kernel = @(t,x) GaussKernel(t,x,th); %without the 4th arg, no transform
   legtxt(ii) = {['\(\theta = ' num2str(th) '\)']};
   for m = 0:mmax
      nskip =2^(mmax-m);
      xdata = xdata_all(1:nskip:nmax+1);
      Kmat = KMP(xdata, xeval, kernel);
      rankKmat(ii,m+1) = rank(Kmat);
      [~,~,rankKmatB(ii,m+1)] = KinvY(Kmat,zeros(size(Kmat,1),1));
   end
end
figure;
h = loglog(nvec,rankKmatB');
legend(h,legtxt,'location','northwest','box','off')
xlabel('\(n\)')
ylabel('rank(\(\mathsf{K}\))')
axis([1 10^ceil(log10(nmax)) 1 10^ceil(log10(nmax))])
print('rankKernel.eps','-depsc')


%% Plot bases functions

thetavec = [0.01 0.1 1 10 100];
nth = size(thetavec,2);
nplot = 1e3;
xeval = (0:nplot)'/nplot;
n = 16;
y0 = zeros(n+1,1);
xdata = (0:n)'/n;
rankKmat(nth,mmax+1) = 0;
rankKmatB(nth,mmax+1) = 0;
legtxt = cell(nth,1);
basevec(n,2) = 0;
for ii = 1:nth
   th = thetavec(ii);
   kernel = @(t,x) GaussKernel(t,x,th); %without the 4th arg, no transform
   legtxt(ii) = {['\(\theta = ' num2str(th) '\)']};
   [Kmat,Kdateval] = KMP(xdata, xeval, kernel);
   [~,cutoff,nok,Vok,Sok,Vnot,Snot] = KinvY(Kmat,y0);
   rankKmatB(ii) = nok;
   basesVis = Kdateval'*Vok ./sqrt(Sok');
   basesInVis = Kdateval'*Vnot ./sqrt(cutoff);
   [~,wh] = max(abs(basesVis),[],1);
   basesVis = basesVis.*sign(basesVis(sub2ind([nplot+1 nok],wh,1:nok)));
   if nok < n
      [~,wh] = max(abs(basesInVis),[],1);
      basesInVis = basesInVis.*sign(basesInVis(sub2ind([nplot+1, n+1-nok],wh,1:n+1-nok)));
   end
   basevec(1,1) = 1;
   scale0 = max(abs(basesVis(:,1)));
   kk = 1;
   for jj = 2:nok
      scale = max(abs(basesVis(:,jj)));
      if scale < 0.1*scale0
         basevec(kk,2) = jj-1;
         kk = kk+1;
         basevec(kk,1) = jj;
      end
   end
   basevec(kk,2) = nok;
   subtitleTxt(1,:) = "\(n = "+int2str(n+1)+"\), rank \(="+int2str(nok)+"\)";
   for jj = 1:kk
      figure
      if basevec(jj,1) == basevec(jj,2)
         whbas = basevec(jj,1);
         plot(xeval,basesVis(:,whbas));
         subtitleTxt(2,:) = "basis "+num2str(whbas);
      else
         whbas = basevec(jj,1):basevec(jj,2);
         plot(xeval,basesVis(:,whbas));
         subtitleTxt(2,:) = "bases "+num2str(whbas(1))+ ...
            " to "+num2str(whbas(end));
      end
      xlabel('\(x\)')      
      title(['Visible, \(\theta = ' num2str(th) '\)'],subtitleTxt) 
      ax = gca;
      ax.TitleHorizontalAlignment = 'right';
      fudge = 0.05*max(max(abs(basesVis(:,whbas))));
      axis([0 1 min(min(basesVis(:,whbas)))-fudge ...
         max(max(basesVis(:,whbas)))+fudge])
   end
   figure
   if nok <n 
      plot(xeval,basesInVis);
      xlabel('\(x\)')
      subtitleTxt(2,:) = "bases "+num2str(nok+1)+ ...
            " to "+num2str(n);
      title(['Invisible, \(\theta = ' num2str(th) '\)'],subtitleTxt) 
      ax = gca;
      ax.TitleHorizontalAlignment = 'right';
      fudge = 0.05*max(max(abs(basesInVis)));
      axis([0 1 min(min(basesInVis))-fudge max(max(basesInVis))+fudge])
   end
end
return
      

      
      
   

%% Test function
f = @(x) exp(-5*x).*cos(10*x);
th = 1;
n = 64;
cutoff = 1e-10;
ntest = 1e4;
xdata = (0:(n-1))'/(n-1);
xeval = (0:(ntest-1))'/(ntest-1);
fdata = f(xdata);
feval = f(xeval);
kernel = @(t,x) GaussKernel(t,x,th); %without the 4th arg, no transform
[Kmat, Kdateval, Kdiageval, errKNull] = KMP(xdata, xeval, kernel);
[errKXx, errKX] = powerfun(Kmat, Kdateval, Kdiageval, cutoff);
AX = ABfun(errKX,errKNull,1,0.2);
[Appx, AppxNorm, ErrBdx, ErrBd] = ...
   Approx(fdata, Kmat, Kdateval, errKXx, errKX, AX, cutoff);
trueErrX = feval - Appx;
trueErr = max(abs(trueErrX))
ErrBd

figure(1)
plot(xeval,feval, ...
   xeval,Appx, ...
   xdata,fdata,'.')

figure(2)
semilogy(xeval,abs(trueErrX),...
   xeval,ErrBdx)
   