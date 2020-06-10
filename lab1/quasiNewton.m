function [xout, info] = quasiNewton(func, x0,rule, opt,varargin)
% this file implements the quasi newton method with modification.
%
% USAGE
% -----------------------------------------------
% [xout,info] = dampedNewton(func, x0, rule, opt)
% [xout,info] = dampedNewton(func, x0, rule, opt, tol)
% [xout,info] = dampedNewton(func, x0, rule, opt, tol, maxiter)
% [xout,info] = dampedNewton(func, x0, rule, opt, tol, maxiter, maxsubiter)
% [xout,info] = dampedNewton(func, x0, rule, opt, tol, maxiter, maxsubiter,
%                       modification)
% [xout,info] = dampedNewton(func, x0, rule, opt, tol, maxiter, maxsubiter,
%                       modification, phi)
% 
% INPUT
% -----------------------------------------------
% func [struct]:
%       func.f [function vector -> double] : target function
%       func.g [function vector -> vector]: gradient of target function, use numerical gradient if is
%                   lacked 
%       func.h [function vector -> matrix]: hessian of target function, use
%               numerical gradient if it is lacked.
% d [vector]: direction, can be replaced by [] when dim = 1
% rule [struct]: rule for line search 
%       rule.name [string]: 
%              - "goldstein"  (default)
%              - "wolfe"
%              - "armijo"
%              - "s-wolfe"
%       rule.rho [double]: (default = 0.3)
%       rule.sig [double]: used for `wolfe` and `s-wolfe`
% opt [string]: methods for getting new stepsize
%       - "none": simple newton update without line search
%       - "bt": backtracking 
%       - "quad": quadratic interpolation (default)
%       - "cubic": cubic interpolation
% tol [double]: tolerence, used in stop criterion: |x_{k+1} - x_k| \le tol
%                    (default = 1e-6)
% maxiter [int] : max iteration number (default = 100)
% maxsubiter [int] : max iteration number in line search sub-routine
%                            (default = 100)
% modification [string] : modification technique
%       - "none" (default)
%       - "LM"
%       - "mixed"
% phi [double] : the parameter of broyden's family
%
% OUTPUT
% -----------------------------------------------
% xout [vector]: the output x 
% info [struct]:  info for damped newton
%       info.success [bool]: -0: success -1: fail
%       info.reason [string]: the reason why program is stopped.
%       info.xn [cell{vec}] : the sequence of x_k 
%       info.f [double] : f(x)
%       info.g [vec] : g(x)
%       info.h [matrix] : h(x)
%       info.iter [int] : outer loop iteration time
%       info.subiter [int] : total iteration time in line search
%                                   sub-routine
%       info.time : run time
%
%
% DESCRIPTION
% -------------------------------------------
% Basic Routine:
%       While |x_{k+1} - x_k| >= tol
%           update hk and quasi newton newton direction  dk = -hk\gk
%           modify dk and normalize dk = dk/|dk|
%           linesearch with direction dk to get x_{k+1}
%       end while
%
% Modification
%       none
%       mixed : if dk'*gk >= 0.3 |dk||gk| then dk =-dk
%                    if |dk'*gk| < 0.3 |dk||gk| then dk = -gk
%      LM:  to find a proper nu >1 such that dk = -(hk + nu*I)\gk such that
%              dk'*gk < -0.3|dk||gk|
%
% Update rule of hk:
%     h_dfp = hk + (sk*sk')/(sk'*yk) - (hk*(yk*yk')*hk)/(yk'*hk*yk);
%     h_bfgs = hk + (1 + (yk'*hk*yk)/(yk'*sk))*(sk*sk')/(sk'*yk) - (sk*yk'*hk + hk*yk*sk')/(yk'*sk);
%     hk = h_dfp +phi*(h_bfgs - h_dfp);
%
%
% Latest Update: 2020/04/01
% Author: Linting
tic
global call
call  = 0;


if isfield(func, 'f'); f = func.f; else; error("No objective!"); end
if ~isfield(func, 'g');  func.g = @(x)numgrad(f,x);end
if ~isfield(func, 'h'); func.h = @(x)numgrad(func.g,x); end
g = func.g; h = func.h;

if nargin < 5; tol = 1e-8; else; tol = varargin{1}; end
if nargin < 6; maxiter = 100; else; maxiter = varargin{2}; end
if nargin < 7; maxsubiter = 100; else; maxsubiter = varargin{3}; end
if nargin <8; modification = "none"; else; modification = varargin{4}; end
if nargin <9; phi = 0; else; phi = varargin{5}; end

if isempty(tol); tol = 1e-8; end
if isempty(maxiter); maxiter = 100; end
if isempty(maxsubiter); maxsubiter = 100; end
if isempty(phi); phi = 0 ;end
if isempty(opt); opt = "quad"; end
iter = 1;
x{iter} = x0;
subiter = 0;
gk = g(x{iter});
hk = h(x{iter});
hk = (hk + hk')/2;
for iter = 1:maxiter
    if iter>1  && norm(x{iter} - x{iter-1}) < tol; break; end
    if all(isfinite(x{iter})) ~= 1; break; end
    %fk = f(x{iter});
    d = -hk\gk;
    
    if modification == "mixed"
        if d'*gk > 0.3*norm(d)*norm(gk)
            d = -d;
        elseif all(isfinite(d)) ~= 1 || d'*gk > -0.3*norm(d)*norm(gk)
            d = -gk;
        end
    end
    
    if modification == "LM"
        nu = 1;
        while all(isfinite(d)) ~= 1 || d'*gk > -0.3*norm(d)*norm(gk)
            d = -(hk+nu*eye(size(hk,1)))\gk;
            nu = nu * 2;
        end
    end
    if isempty(opt) || opt ~= "none"
        d = d/norm(d);
        [~, subinfo] = linesearch(func, x{iter}, d, rule, opt, 1, maxsubiter);
        if subinfo.success == 1
            subiter = subiter + subinfo.iter;
            x{iter+1} = subinfo.x;
        else
            %fprintf("Line Search Failed! Turn to exact line search\n");
            [~, subinfo] = linesearch618(func, x{iter}, d);
            subiter = subiter + subinfo.iter;
            x{iter+1} = subinfo.x;
        end
    else
        x{iter+1} = x{iter}+d;
    end
    gkp1 =  g(x{iter+1});
    if(gkp1==0); break; end
    yk = gkp1 - gk;
    gk = gkp1;
    sk = x{iter+1} - x{iter};
    % update hk by broyden's method.
    if ~isnumeric(phi)  % speical case for SR1
        hk = hk + ((sk-hk*yk)*(sk-hk*yk)')/((sk-hk*yk)'*yk);
        continue
    end
    h_dfp = hk + (sk*sk')/(sk'*yk) - (hk*(yk*yk')*hk)/(yk'*hk*yk);
    h_bfgs = hk + (1 + (yk'*hk*yk)/(yk'*sk))*(sk*sk')/(sk'*yk) - (sk*yk'*hk + hk*yk*sk')/(yk'*sk);
    hk = h_dfp +phi*(h_bfgs - h_dfp);
end


xout = x{iter};
if iter == maxiter
    info.success = 0;
    info.reason = "max iteration reached!";
elseif  all(isfinite(xout)) ~= 1% || norm(x{iter} - x{iter-1}) > 1
    info.success = 0;
    info.reason = "quasi Newton¡®s method failed!";
else
    info.success = 1;
    info.reason = "\|x_{k+1} - x_k\| < "+ string(tol);
    
end

info.iter = iter;
info.subiter = subiter;
info.xn = x;
info.f = func.f(x{iter});
info.g = func.g(x{iter});
info.h = func.h(x{iter});
info.time = toc;

info.call = call;
