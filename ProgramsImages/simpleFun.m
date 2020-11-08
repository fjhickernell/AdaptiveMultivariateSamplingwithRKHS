function [y,fun_name] = simpleFun(x)
y = 2 * (exp(-6.*(x-1/2).^2) - 0.6);
if nargout > 1
   fun_name = 'SimpleFunction';
end
