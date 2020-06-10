%This SCRIPT file stores some test functions for numerical experiment
%
% STRUCTURE
% -----------------------------------------------
% func [struct]
%       func.f : The Objective fuction.
%       func.g: The gradient function.
%       func.h: The Hessian.
%
% OUTPUT
% -----------------------------------------------
% bd{n} [`func`]: brown and dennis function
% die{n} [`func`]: discrete integral equation
% dieinit [function int -> vector]: generate the initial value of die problem
% mse{n} [`func`]: minimal surface equation
% mseinit [function int -> vector] : generate the initial value of mse problem
%
%
% Last Update: 2020/04/01
% Author: Linting

clear all %#ok
simple.f = @(x) (x(1)^1 + 100*x(2)^2);
simple.g = @(x) [2*x(1) + 2*(x(1) - sin(x(2))); 200*x(2)];
simple.h = @(x) [2 0; 0 200];

%% main part 
% toy example  -- rosenbrock function
rosenbrock.f = @(x) (1-x(1))^2 + 100*(x(2) - x(1)^2)^2;
rosenbrock.g = @(x) [2*(x(1)-1) + 200*(x(1)^2-x(2))*(2*x(1)); 200*(x(2)-x(1)^2)];
rosenbrock.h = @(x) [2+200*(2*x(1))*(2*x(1)) + 400*(x(1)^2-x(2)), -400*x(1);-400*x(1), 200];

% brown and dennis function
bd4.f = @(x)BrownDennisFunc(x,4);
bd10.f = @(x)BrownDennisFunc(x,10);
bd20.f = @(x)BrownDennisFunc(x,20);
bd30.f = @(x)BrownDennisFunc(x,30);
bd40.f = @(x)BrownDennisFunc(x,40);
bd50.f = @(x)BrownDennisFunc(x,50);

% discrete integral equation
die2.f = @(x) DIE(x,2);
die10.f = @(x) DIE(x,10);
die20.f = @(x) DIE(x,20);
die30.f = @(x) DIE(x,30);
die40.f = @(x) DIE(x,40);
die50.f = @(x) DIE(x,50);
dieinit = @DIEinit;

% minimal surface equation
mse3.f = @(x) MSE(x,3);
mse5.f = @(x) MSE(x,5);
mse7.f = @(x) MSE(x,7);
mseinit = @MSEinit;

%% Sub-program
function fx = BrownDennisFunc(x,m)
global call
call = call + 1;
fx = 0;
for i = 1:m
    t = i/5;
    fx = fx + ((x(1) + t*x(2) -exp(t))^2 + (x(3) + x(4)*sin(t) - cos(t))^2)^2;
end
end

function fx = DIE(x,n)
global call
call = call + 1;
h = 1/(n+1);
fx = 0;
i = (1:n); j = (1:n)';
fxx = h*(1-i*h).*(j*h).*(x + j*h + 1).^3/2 .* (j<=i) + h*(i*h).*(1 - j*h).*(x + j*h + 1).^3/2 .* (j>i);
fxx = sum(fxx,1)' + x;
fx = sum(fxx.^2);
end

function x0 = DIEinit(n)
t = (1:n)'/n+1;
x0 = t.*(t-1);
end


function fx = MSE(x,n)   % mse with boundary value sin
global call
call = call + 1;
    A = ((0:n)'/n).^2 + ((0:n)/n).^2;
    h = 1/n;
    A(2:n,2:n) = reshape(x,n-1,n-1);
    B = h^2 * (1 + (A(2:n,2:n) - A(2:n,1:n-1)).^2/h^2 + (A(2:n,2:n) - A(1:n-1,2:n)).^2/h^2).^(1/2);
    B = B + h^2 * (1 + (A(2:n,2:n) - A(2:n,3:n+1)).^2/h^2 + (A(2:n,2:n) - A(3:n+1,2:n)).^2/h^2).^(1/2);
    fx = sum(sum(B));
end
    

function x0 = MSEinit(n)
    x0 = rand((n-1)*(n-1),1);   % randomized initial value to avoid starting from zero.
end
