function x = seqFixedDes(nrange, d, shift)
%SEQFIXEDDES gives an extensible lattice design that is shifted
if ~exist('d','var')
   if ~exist('shift','var')
      d = 1;
      shift = 1/3;
   end
else
   if ~exist('shift','var')
      shift = 1/3*ones(1,d);
   end
end
if ~exist('nrange','var')
   nrange = 1;
end


