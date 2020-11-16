function c = KinvY(Kmat,y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[V,S,~] = svd(Kmat,'econ');
Sdiag = diag(S);
nok = find(Sdiag > 1e-12,1,'last');
Vok = V(:,1:nok);
c = Vok*(Vok'*y ./ Sdiag(1:nok));
end

