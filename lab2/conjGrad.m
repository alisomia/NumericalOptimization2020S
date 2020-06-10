function [xk, info] = conjGrad(func, x0, LineSearchRule, LineSearchOption, DirectionUpdateOption, maxiter, maxsubiter, tol, RestartLevel)
% This function implements nonlinear conjugate gradient method.
%
%
% USAGE
%-------------------------------------------------
% [xk, info] = conjGrad(func, x0, LineSearchRule, LineSearchOption, DirectionUpdateOption, maxiter, maxsubiter, tol，RestartLevel)
%   
%
% INPUT
% ------------------------------------------------
% func [struct]:
%       func.f [function vector -> double] : target function
%       func.g [function vector -> vector]: gradient of target function, use numerical gradient if is
%                   lacked 
% x0 [double]:  Initial Point (default = func.init)
% LineSearchRule [struct]: rule for line search 
%       rule.name [string]: 
%              - "goldstein"  (default)
%              - "wolfe"
%              - "armijo"
%              - "s-wolfe"
%       rule.rho [double]: (default = 0.3)
%       rule.sig [double]: used for `wolfe` and `s-wolfe`
% LineSearchOption [string]: methods for getting new stepsize
%       - "bt": backtracking 
%       - "quad": quadratic interpolation (default)
%       - "cubic": cubic interpolation
%  DirectionUpdateOption [string/int]: methods of updating direction
%       - "FR"
%       - "PRP"
%       - "PRP+"
%       - "FR-PRP"
%       - [n] : Hu-Storey method with restart thereshold n.
%
%
%
% maxiter [int] : max iteration number (default = 100)
% maxsubiter [int] : max iteration number for line search (default = 5)
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
% The algorithm start with dk = -gk
% WHILE not converge DO
%       LineSearch to find x_{k+1}, compute gkp1 = g_{k+1}
%       Update d_{k+1}
% ENDWHILE
%
%   The LINESEARCH step is used s-wolfe or exact line search method
%   Update rule is controlled by DirectionUpdateOption
%       all but HS algorithm follows $$d_{k+1} = -g_{k+1} +beta * dk$$
%           case "FR"
%                beta = (gkp1'*gkp1)/(gk'*gk);
%            case "PRP"
%                beta =  (gkp1'*gkp1 - gkp1'*gk)/(gk'*gk);
%            case "PRP+"
%                beta =  max((gkp1'*gkp1 - gkp1'*gk)/(gk'*gk),0);
%            case "FR-PRP"
%                beta_FR =  (gkp1'*gkp1)/(gk'*gk);
%                beta_PRP = (gkp1'*gkp1 - gkp1'*gk)/(gk'*gk);
%                if abs(beta_PRP) < beta_FR; beta = beta_PRP; end
%                if  beta_PRP < - beta_FR; beta = - beta_FR; end
%                if beta_PRP >beta_FR; beta = beta_FR; end
%   The HS algorithm is based on solving the following ssystem:
%
%
%
%   For the details and improvements, see [Hu and Storey, 1991]
%
%
%  Linting @ PKU 
% 20200418


%% set default args
if isempty(tol); tol = 1e-6; end
if isempty(maxiter); maxiter = 100; end
if isempty(maxsubiter); maxsubiter = 6; end
if isempty(LineSearchOption); LineSearchOption = "quad"; end
if isempty(LineSearchRule); LineSearchRule.name = "s-wolfe"; LineSearchRule.rho = 0.3; LineSearchRule.sigma = 0.6; end
if isempty(DirectionUpdateOption); DirectionUpdateOption = "FR"; end
if isempty(x0); x0 = func.init; end 

if isfield(func, 'f'); f = func.f; else; error("No objective!"); end
if ~isfield(func, 'g');  func.g = @(x)numgrad(f,x);end
%if DirectionUpdateOption == "HS"; DirectionUpdateOption  = 20; end
g = func.g;

tic
xk = x0;
gk = g(xk);
fk = f(xk);
dk = -gk;
subiter = 0;
HScounter = 0;
for iter = 1:maxiter
    if norm(gk,inf) < tol*(1 + abs(fk)); break; end
    if ~all(isfinite(gk)); break; end
    if gk'*dk>0; dkk = -dk/norm(dk);else; dkk = dk/norm(dk); end
    if maxsubiter > 0
        [alpha, subinfo] = linesearch(func, xk, dkk, LineSearchRule, LineSearchOption, 1, maxsubiter);
    else 
        [alpha, subinfo] = linesearch618(func, xk, dkk);
    end
%         if subinfo.success == 1
%             subiter = subiter + subinfo.iter;
%             xk = subinfo.x;
%         else
%             [alpha, subinfo] = linesearch618(func, xk, dk/norm(dk)); %#ok
%             subiter = subiter + subinfo.iter;
%             xk = subinfo.x;
%         end
    subiter = subiter + subinfo.iter;
    xk = subinfo.x;
       fkp1 = f(xk);
       gkp1 = g(xk);
       if(isstring(DirectionUpdateOption))
       switch DirectionUpdateOption
           case "FR"
               beta = (gkp1'*gkp1)/(gk'*gk);
           case "PRP"
               beta =  (gkp1'*gkp1 - gkp1'*gk)/(gk'*gk);
           case "PRP+"
               beta =  max((gkp1'*gkp1 - gkp1'*gk)/(gk'*gk),0);
           case "FR-PRP"
               beta_FR =  (gkp1'*gkp1)/(gk'*gk);
               beta_PRP = (gkp1'*gkp1 - gkp1'*gk)/(gk'*gk);
               if abs(beta_PRP) < beta_FR; beta = beta_PRP; end
               if  beta_PRP < - beta_FR; beta = - beta_FR; end
               if beta_PRP >beta_FR; beta = beta_FR; end
       end
       else
               HScounter = HScounter + 1;
               if HScounter <= DirectionUpdateOption
               eps = 1e-6;
               Hdk = (g(xk+eps*dk) - g(xk -eps*dk))/(2*eps);
               Hgk = (g(xk+eps*gkp1) - g(xk -eps*gkp1))/(2*eps);
               tk = dk'*Hdk;
               vk = gkp1'*Hgk;
               uk = gkp1'*Hdk;
               gkp1sq = gkp1'*gkp1;
               dksq = dk'*dk;
               gkp1dk = gkp1'*dk;
               end
               if HScounter <= DirectionUpdateOption %&& tk>0 && vk>0 && 1-uk^2/(tk*vk) > 1/4 - eps && (vk/gkp1sq)< tk/dksq
                   wk = tk*vk-uk^2;
                   dk = 1/wk*((uk*gkp1dk - tk*gkp1sq)*gk + (uk*gkp1sq - vk*gkp1dk)*dk);
               else
                   HScounter = 1;
                   dk = -gkp1;
               
               end
       end
       fk = fkp1;
       gk = gkp1;
       if isstring(DirectionUpdateOption)
           if RestartLevel~=0 && mod(iter, RestartLevel) ==0
               dk = -gk;
           else
               dk = -gk +beta * dk;
           end
        end
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
