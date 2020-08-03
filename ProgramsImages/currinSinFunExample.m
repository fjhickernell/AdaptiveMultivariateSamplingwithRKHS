%% Currin sinusoidal
S = struct('type','{}','subs',{{':'}}); prm = [];
[prm(1:3).AlgName] = subsref({'Algo1', 'Algo2', 'Algo3'},S);
[prm,kernelth] = parseFunAppxParam(prm);
prm(3).n0 = 10;
nAlg = size(prm,2);
f = @(x) sin(2*pi*(x-0.1));
[prm.fname] = subsref(repmat({'CurrinSineFun'},1,nAlg),S);
[prm.yLim] = subsref(repmat({[-2;1.5]},1,nAlg),S);

%%
RunExample(f,prm,kernelth)

%% With adpatpive theta 
[prm.whDes] = subsref(repmat({'unif_grid'},1,nAlg),S);
RunExample(f,prm,kernelth)

   