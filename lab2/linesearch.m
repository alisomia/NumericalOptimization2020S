function [alpha, info] = linesearch(func, x0,  d, rule, opt,varargin)
% This function implements the line search methods with various rule
%
% USAGE
%-------------------------------------------------
% [alpha, info] = linesearch(func, x0, d, rule, opt)
% [alpha, info] = linesearch(func, x0, d, rule, opt, alpha0)
% [alpha, info] = linesearch(func, x0, d, rule, opt, alpha0, maxiter)
%
%
% INPUT
% ------------------------------------------------
% func [struct]:
%       func.f [function vector -> double] : target function
%       func.g [function vector -> vector]: gradient of target function, use numerical gradient if is
%                   lacked 
% x0 [double]:  Point
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
%       - "bt": backtracking 
%       - "quad": quadratic interpolation (default)
%       - "cubic": cubic interpolation
% alpha0 [double]: initial step size (default = 1)
% maxiter [int] : max iteration number (default = 100)
%
% OUTPUT
%---------------------------------------------
% alpha [double]: step size 
% info [struct]:  info for line search
%       info.success [bool]: -0: success -1: fail
%       info.reason [string]: the reason why program is stopped.
%       info.x [vec] : x0 + alpha*dir
%       info.fval [double] : f(x)
%       info.gval [vec] : g(x)
%
% DESCRPITION
%------------------------------------------------
% This Program implements the simple linesearch method, 
% Basic Routine:
%       while alpha does not satisfies rule
%           gen new alpha from opt
%       end while
% Rules:
%       Goldstein: f(alpha) - f0 < rho*g0'*d*alpha && f(alpha) - f0 > (1-rho)*g0'*d*alpha 
%       Armijo: f(alpha) - f0 < rho*g0'*d*alpha
%       Wolfe : f(alpha) - f0 < rho*g0'*d*alpha && g(alpha)'d > sig * g0'*d
%       Strong Wolfe: f(alpha) - f0 < rho*g0'*d*alpha && |g(alpha)'d| < - sig * g0'*d
%
% Options: to generate new stepsize
%       BackTracking: alpha <- alpha *0.9
%       Quadratic Interpolation: Using phi(0), phi'(0), phi(alpha) to
%                                               interpolate a quad
%                                               function, with its minimum being new alpha 
%       Cubic Interpolation: Using phi(0), phi'(0), phi(alpha), phi(alpha1)
%                                       to interpolate a cubic function,
%                                       with its minimum alpha2
 %                                      Then [alpha alpha1] <- [alpha1
 %                                      alpha2]
% Ver: 0.0
% Author: Linting


%% START INPUT CHECKING

% checking if `func` is valid
if isfield(func,'f')
    f = @(a) func.f(x0 + a*d) - func.f(x0);
else
    error('LineSearch InputError: No Objective is given!');
end

if isfield(func,'g')
    g = @(a) func.g(x0 + a*d);
else
    error('LineSearch InputError: No Gradient is given!');
end

if isempty(rule)
    rule.name = "s-wolfe";
    rule.sigma = 0.6;
    rule.rho = 0.3;
end

switch rule.name
    case "goldstein"
       checkAlpha = @checkAlpha_goldstein;
       rule.sigma = 0;
    case "armijo" 
        checkAlpha = @checkAlpha_armijo;
        rule.sigma = 0;
    case "wolfe"
        checkAlpha = @checkAlpha_wolfe;
    case "s-wolfe"
        checkAlpha = @checkAlpha_s_wolfe;
    otherwise
        error('LineSearch InputError: No Rules matched!');
end

if nargin <6; alpha = 1; else; alpha = varargin{1}; end % set default initial stepsize

if nargin <7; maxiter = 100; else; maxiter = varargin{2};end % set default max iter

% END INPUT CHECKING

%% MAIN CODE

    g0d = g(0)'*d;
    if opt == "cubic"; alpha1 = g0d*alpha^2*.5/(f(alpha) - g0d*alpha); end
    
    for iter = 1: maxiter
        fval = f(alpha);
        if checkAlpha(fval,g0d,alpha,d,rule.rho, rule.sigma,g); break; end
        switch opt
            case "bt"
                alpha = alpha * 0.5;
            case "quad"
                alpha = -g0d*alpha^2*.5/(fval - g0d*alpha);
            case "cubic"
                fval1 = f(alpha1);
                p = [alpha^2 -alpha1^2; -alpha^3 alpha1^3] * [fval1 - g0d*alpha1;fval - g0d*alpha] ...
        / ( alpha^2 * alpha1^2 *(alpha1-alpha) );
                p = roots([p.*[3; 2]; g0d]);
                if isreal(p); alpha2 = max(p); end
                alpha = alpha1;
                alpha1 = alpha2;
        end
    end
    


    
%% OUTPUT
if iter==maxiter
    info.success = 0;
    info.reason = "Failed after Max Iteration Reached!";
else
    info.success = 1;
    info.reason =  " Success after "+ string(iter) + " iteration: " ;
end
    info.rule = rule.name;
    info.opt = opt;
    info.iter = iter;
    info.x = x0 + alpha * d;
    info.fval = fval + func.f(x0); 
    info.gval = g(alpha);

end

%% SUB PROGRAM
% various rules to check if the stepsize is valid
function ok = checkAlpha_goldstein(fk,g0d,alphak,~,rho,~,~)
    gdalpha = g0d*alphak;
    ok = (fk < gdalpha*rho) &&(fk > gdalpha*(1-rho));
end

function ok = checkAlpha_armijo(fk,g0d,alphak,~,rho,~,~)
    gdalpha = g0d*alphak;
    ok = (fk < gdalpha*rho);
end

function ok = checkAlpha_wolfe(fk,g0d,alphak,d,rho,sigma,g)
gdalpha = g0d*alphak;
    ok = (fk < gdalpha*rho) && (g(alphak)'*d > sigma*g0d);
end

function ok = checkAlpha_s_wolfe(fk,g0d,alphak,d,rho,sigma,g)
gdalpha = g0d*alphak;
    ok = (fk < gdalpha*rho) && (abs(g(alphak)'*d) < - sigma*g0d);
end
