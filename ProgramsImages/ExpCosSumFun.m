function y = ExpCosSumFun(x,a)
y = exp(sum(a(2,:).*x,2)).*cos(sum(a(1,:).* x,2));
