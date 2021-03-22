function [c,nok,Vok,Sok,Vnot,Snot] = KinvY(Kmat,y,cutoff)
%KINVY multiplies a vector times the inverse of the Gram matrix
n = size(y,1);
if nargin <= 2
   cutoff = n*eps;
elseif cutoff <=0
   cutoff = n*eps;
end
% Kmat =  Kmat + n*eps*eye(n);
% rcondKnug = rcond(Kmat)
% rankKnug = rank(Kmat)
[V,S,~] = svd(Kmat,'econ');
Sdiag = diag(S);
nok = find(Sdiag > cutoff,1,'last');
Vok = V(:,1:nok);
Sok = Sdiag(1:nok);
c = Vok*(Vok'*y ./ Sok);
if nargout > 4
   Vnot = V(:,nok+1:n);
   Snot = Sdiag(nok+1:n);
end

