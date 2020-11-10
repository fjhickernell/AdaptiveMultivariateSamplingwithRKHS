function xdata = fixedDesign(idx,obj)
%FIXEDDESIGN provides the points for a fixed, non-adpative in position
%design
n = max(idx);
nidx = length(idx);
if strncmp(obj.whDes,'unif',1) && n > 1 && nidx == n && obj.dim ==1
   temp = (idx-1)'/(n-1);
   if strcmp(obj.whDes,'unifChebyshev') %Chebyshev
      temp = (1+sin(pi*(-1/2 + temp)))/2;
   end
elseif strcmp(obj.whDes,'seqChebyshev') %sequential Chebyshev
   temp = seqFixedDes(idx,obj.dim,1/3,'Chebyshev');
else %sequential
   temp = seqFixedDes(idx,obj.dim);
end
xdata = obj.xLim(1,:) + (obj.xLim(2,:) - obj.xLim(1,:)).*temp;


