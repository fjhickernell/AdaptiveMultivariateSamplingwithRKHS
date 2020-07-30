%% Currin sinusoidal
f = @(x) sin(2*pi*(x-0.1));
[prm.fname] = subsref(repmat({'CurrinSineFun'},1,nAlg),S);
[prm.yLim] = subsref(repmat({[-2;1.5]},1,nAlg),S);
RunExample(f,prm,abstolVec,kernelth)

%% With adpatpive theta 
[prm.whDes] = subsref(repmat({'adapt_th'},1,nAlg),S);
RunExample(f,prm,abstolVec,kernelth)

   