classdef FunctionApproxProblem
   %FUNCTIONAPPROXIMATIONPROBLEM sets up a function approximation problem
   %to be solved using our methods
   
   properties
      forig = @simpleFun
      fparam
      kernelOrig = @GaussKernel
      Algo = @AdaptAlgo1
      isDiagnose = true
      Ainf = 0.5
      B0 =  0.2
      n0 = 1
      nmax = 500
      neval = 20000
      theta = 1
      legendPos = 'south'
      whDes = 'adapt_th'
      whObj = 'EmpBayesAx'
      plotSites = false
      xLim = [0;1]
      yLim = [-1;1]
      thetaRange = (-5:0.5:5)'
      abstolVec = [0.05 0.02 0.01 0.005 0.002 0.001]'
      canvasTheta = false
      currentTheta = 0
      final_theta
      colorScheme
   end
   
   properties (Dependent)
      f
      dim
      fname
      xeval
      kernelth
      kername
      algoname
   end
   
   methods
      function obj = FunctionApproxProblem(inputFun)
         %FUNCTIONAPPROXPROBLEM Construct an instance of this class
         %   Detailed explanation goes here
         InitializeDisplay
         colorScheme = [MATLABBlue; MATLABOrange; MATLABGreen; MATLABPurple; ...
            MATLABCyan; MATLABMaroon];
         obj.colorScheme = colorScheme;
         if nargin > 0
            d_obj = length(inputFun);
            obj(1,d_obj) = obj;
            for kk = 1:d_obj
               obj(kk).forig = inputFun{kk}; %our nice color scheme
            end
         end
      end
      
      function obj = set_prop(obj,propname,val,idx)
         S = struct('type','{}','subs',{{':'}});
         if nargin < 4
            idx = 1:length(obj);
         end
         nval = length(val);
         nidx = length(idx);
         [obj(idx).(propname)] = subsref(repmat(val,1,nidx/nval),S);
      end
      
      function val = get.dim(obj)
         val = size(obj.xLim,2);
      end
      
      function val = get.f(obj)
         if isempty(obj.fparam)
            val = obj.forig;
         else
            val = @(x) obj.forig(x,obj.fparam);
         end
      end
      
      function val = get.fname(obj)
         val = functions(obj.forig).function; %get the name of the function
      end
      
      function val = get.kernelth(obj)
         val = @(t,x,th) obj.kernelOrig(t,x,th,true);  %set up kernel
      end
      
      function val = get.kername(obj)
         val = functions(obj.kernelOrig).function;  %and its name
      end
      
      function val = get.algoname(obj)
         val = functions(obj.Algo).function;  %and the algorithm name
      end
      
      function val = get.xeval(obj)
         switch obj.dim  %set the places where we will we look for the next data
            case 1
               val = obj.xLim(1) + (obj.xLim(2) - obj.xLim(1))* ...
                  (0:obj.neval-1)'/(obj.neval-1);
            case 2
               ngrid = floor(sqrt(obj.neval));
               nngrid = (0:ngrid-1)'/(ngrid-1);
               x = obj.xLim(1,1) + (obj.xLim(2,1) - obj.xLim(1,1))*nngrid;
               y = obj.xLim(1,2) + (obj.xLim(2,2) - obj.xLim(1,2))*nngrid;
               [xx, yy]= meshgrid(x,y);
               val = [xx(:),yy(:)];
            otherwise
               meval = floor(log2(obj.neval));
               val = seqFixedDes([1 2^meval], obj.dim);
         end
      end


   end
end

