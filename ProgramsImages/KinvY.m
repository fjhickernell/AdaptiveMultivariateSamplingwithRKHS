function c = KinvY(Kmat,y)
%KINVY multiplies a vector times the inverse of the Gram matrix
[V,S,~] = svd(Kmat,'econ');
Sdiag = diag(S);
nok = find(Sdiag > 1e-12,1,'last');
Vok = V(:,1:nok);
c = Vok*(Vok'*y ./ Sdiag(1:nok));
end

