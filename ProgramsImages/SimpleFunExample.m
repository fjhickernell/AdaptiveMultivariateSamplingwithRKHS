%% Example of easy function

f = @simpleFun;
[prm.fname] = subsref(repmat({'SimpleFun'},1,nAlg),S);
RunExample(f,prm,abstolVec,kernelth)

%% With uniform grid 
[prm.whDes] = subsref(repmat({'unif_grid'},1,nAlg),S);
RunExample(f,prm,abstolVec,kernelth)


   