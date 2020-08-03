%% Example of easy function
tic
S = struct('type','{}','subs',{{':'}}); prm=[];
[prm(1:3).AlgName] = subsref({'Algo1', 'Algo2', 'Algo3'},S);
%[prm.AlgName] = 'Algo3';
[prm,kernelth] = parseFunAppxParam(prm);
prm(3).n0 = 10;
nAlg = size(prm,2);
f = @simpleFun;
[prm.fname] = subsref(repmat({'SimpleFun'},1,nAlg),S);

%%
RunExample(f,prm,kernelth)

%% With uniform grid 
[prm.whDes] = subsref(repmat({'unif_grid'},1,nAlg),S);
RunExample(f,prm,kernelth)

toc


   