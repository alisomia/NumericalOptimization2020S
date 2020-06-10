function r = gradtest(func)
x = randn(size(func.init,1),1);
f = numgrad(func.f, x) ./ func.g(x);
r =  max(abs(f-1));