classdef FunAppxOut
   %FUNAPPXOUT gives output from the function approximation process
   
   properties
      Appx
      ErrBdx
      ErrBdVec
      trueErr
      InErrBars
      AppxNorm
      NeccFlag
      xdata
      fdata
      thetaOptimalVec
      currentTheta
   end
   
   methods
      function obj = FunAppxOut(FunAppxInpObj)
         %FUNAPPXOUT Construct an instance of this class
         %   Detailed explanation goes here
         if nargin > 0
            d_obj = length(FunAppxInpObj);
            obj(1,d_obj) = obj;
         end
      end
      
   end
end

