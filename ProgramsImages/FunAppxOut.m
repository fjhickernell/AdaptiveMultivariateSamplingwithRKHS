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
   end
   
   methods
      function obj = FunAppxOut(FunAppxInpObj)
         %FUNAPPXOUT Construct an instance of this class
         %   Detailed explanation goes here
         d_obj = length(FunAppxInpObj);
         obj(1,d_obj) = obj;
      end
      
   end
end

