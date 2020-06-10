### README

optimization method: lab 1

author: Ting Lin



-------

**Report** is in `./report`

`./loadTestFunction.m` is a script file for loading all the test function

`./linsearch.m` implements the inexact line search method

`./linsearch618.m` implements the exact line search method

`./dampedNewton.m` implements the damped Newton method

`./quasiNewton.m` implements the quasi Newton method (Broyden's family)

`./numgrad.m` implements the numerical gradient

`./do.m` is an auxiliary file for numerical experiment. 



----------

DEMO

`[x, info] = quasiNewton(die20, dieinit(20), [],"quad",[],10000,[],"mixed",0);`

`[x, info] = dampedNewton(die20, dieinit(20), [],"quad",[],10000,[],"mixed");`