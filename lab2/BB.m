function [xk,info] = BB(func, x0, option,maxiter, tol)
% This function implements BB method with nonmonotone line search.
% ref : The Barzilai and Borwein Gradient Method for 
% the Large Scale Unconstrained Minimization Problem
%
%
% USAGE
%-------------------------------------------------
% [xk, info] = BB(func, x0, option, maxiter, tol)
%   
%
% INPUT
% ------------------------------------------------
% func [struct]:
%       func.f [function vector -> double] : target function
%       func.g [function vector -> vector]: gradient of target function, use numerical gradient if is
%                   lacked 
% x0 [double]:  Initial Point (default = func.init)
% option [struct]: see description
%       option.delta [double] (default = 1)
%       option.sigma [double] (default = 0.6)
%       option.M [int] (default = 5)
%       option.eps [double] (default = 1e-6)
%       option.gamma [double] (default = 0.3)
% maxiter [int] : max iteration number (default = 100)
% tol [double] : tolerance (default = 1e-6)
%
%
% OUTPUT
%---------------------------------------------
% alpha [double]: step size 
% info [struct]:  info
%       info.success [bool]: -0: success -1: fail
%       info.reason [string]: the reason why program is stopped.
%       info.x [vec] : x
%       info.f [double] : f(x)
%       info.g [vec] : g(x)
%
% DESCRIPTION
% ---------------------------------------------
% ak = 1
% WHILE not converge 
%   IF ak>eps OR ak<1/eps THEN ak=delta ENDIF
%   lambda = 1/ak
%   WHILE f(xk - lambda*gk) >= max_{j<min(k,M)}f_{k-j} - gamma*lambda*(gk'*gk)
%       lambda = lambda * sigma
%   END
%   xkp1 = xp - gk*lambda
%   sk = xkp1 - xk;
%   yk = gkp1 - gk;
%   ak = (sk'*yk)/(sk'*sk);
% ENDWHILE
%
%   Linting@PKU
%   20200418
%
if isempty(tol); tol = 1e-6; end
if isempty(maxiter); maxiter = 100; end
if isempty(x0); x0 = func.init; end
if isempty(option)
    option.delta = 1;
    option.sigma = 0.4;
    option.eps = 1e-10;
    option.M = 10;
    option.gamma = 1e-4;
end

delta = option.delta;
sigma = option.sigma;
eps = option.eps;
M = option.M;
gamma = option.gamma;

if isfield(func, 'f'); f = func.f; else; error("No objective!"); end
if ~isfield(func, 'g');  func.g = @(x)numgrad(f,x);end
g = func.g;

tic
flist = ones(1,M) * (-inf);
xk = x0;
gk = g(xk);
fk = f(xk);
ak = 1;
dk = -gk;
subiter = 0;
for iter = 1:maxiter
    if norm(gk,inf) < tol*(1 + abs(fk)); break; end
    flist(mod(iter,M)+1) = fk;
    if norm(gk,2)>1; delta = 1; elseif norm(gk,2)<1e-5; delta = 1e5; else; delta = 1/(norm(gk,2)); end
    if ak<eps || ak > 1/eps; ak = delta;end
    lambda = 1/ak;
    
    % non monotone line search
    for subiterlocal = 1:5
        fk = f(xk - lambda *gk);
        if  fk < max(flist) - gamma*lambda*(gk'*gk); break; end
        lambda = sigma*lambda;
    end
    subiter = subiter + subiterlocal;
    xkp1 = xk-lambda*gk;
%     if lambda<tol
%         [~, subinfo] = linesearch618(func, xk, dk);
%             %subiter = subiter + subinfo.iter;
%             xkp1 = subinfo.x;
%     else
%         xkp1 = xk - lambda*gk;
%     end
    
    gkp1 = g(xkp1);
    
   
    sk = xkp1 - xk;
    yk = gkp1 - gk;
    ak = (sk'*yk)/(sk'*sk);
    xk = xkp1;
    gk = gkp1;
end

if norm(gk,inf) < tol*(1 + abs(fk))
    info.success = 1;
    info.reason = "|g| <= tol * (1 + |f|)";
elseif ~all(isfinite(gk))
    info.success = 0;
    info.reason = "CG failed!";
else
    info.success = 0;
    info.reason = "Max Iter Reached!";
end


info.iter = iter;
info.subiter = subiter;
info.x = xk;
info.f = fk;
info.g = gk;
info.time =toc;

