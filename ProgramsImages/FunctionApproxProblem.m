classdef FunctionApproxProblem
   %FUNCTIONAPPROXIMATIONPROBLEM sets up a function approximation problem
   %to be solved using our methods
   
   properties
      f = @simpleFun
      kernelOrig = @GaussKernel
      Algo = @AdaptAlgo1
      isDiagnose = true
      Ainf = 0.5
      B0 =  0.2
      n0 = 1
      nmax = 500
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
      colorScheme
   end
   
   properties (Dependent)
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
               obj(kk).f = inputFun{kk}; %our nice color scheme
            end
         end
      end
      
      function val = get.dim(obj)
         val = size(obj.xLim,2);
      end
      
      function val = get.fname(obj)
         [~,val] = obj.f(obj.xLim(1,:)); %get the name of the function
      end
      
      function val = get.kernelth(obj)
         val = @(t,x,th) obj.kernelOrig(t,x,th,true);  %set up kernel
      end
      
      function val = get.kername(obj)
         [~,~,~,~,val] = ...
            obj.kernelth(obj.xLim(1,:),obj.xLim(1,:),obj.theta);  %and its name
      end
      
      function val = get.algoname(obj)
         [~,val] = obj.Algo();  %and its name
      end
     function val = get.xeval(obj)
         switch obj.dim  %set the places where we will we look for the next data
            case 1
               val = (0:0.0002:1)';
            case 2
               x = 0:0.005:1;
               [xx, yy]= meshgrid(x);
               val = [xx(:),yy(:)];
            otherwise
               val = seqFixedDes([1 2^13], obj.dim);
         end
      end


   end
end

