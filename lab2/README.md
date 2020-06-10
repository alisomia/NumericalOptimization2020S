### Numreical Optimization Lab 2

Ting Lin 1700010644

----------------------

This lab implements some methods about BB and CG.

`./report` : the report file



`./BB.m` : the BB program

`./conjGrad.m`: the nonlinear CG program



`./linesearch.m`:  inexact line search 

`./linsearch618.m`: exact line search



`./loadTestFunction`: load the test function

`./numgrad.m`: numerical gradient

`./gradtest.m` : test the gradient is valid













-----------------------

Usage:

```matlab
LoadTestFunction
[x,info] = conjGrad(trig4, [], [], [], "PRP+", 10000, 0, 1e-6,0)
%[x, info] = conjGrad(func, x0 = func.init, LineSearchRule = [], LineSerachOption = [], DirectionUpdateOption = "FR", maxiter, maxsubiter, tol, restart)

[x,info] = BB(trig4, [], [],10000, 1e-6)
% [x,info] = BB(func, x0 = func.init, LineSearchRule = [], LineSearchOption = [], maxiter, tol)
```

