function [alpha, info] = linesearch618(func,  x0,  d, varargin)
% an exact line search using .618 method
%
% USAGE
% -----------------------------------------------
% [alpha, info] = linesearch618(func, x0, d)
% [alpha, info] = linesearch618(func, x0, d, tol)
%
%
% INPUT
% ------------------------------------------------
% func [struct]:
%       func.f [function vector -> double] : target function, must be
%       coerovice, otherwise the function will not terminate.
% x0 [double]:  Point
% d [vector]: direction, can be replaced by [] when dim = 1
% tol [double]: the tolerence of .618 method 
% 
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
% DESCRIPTION
%----------------------------------------------
% Basic Routine: 
% The .618 method consists two parts, 
%     £¨I) determine a initial search interval
%       (II) apply .618 method
    % (I): OUTPUT -- (l,r) such that there exists m such that phi(l), phi(m),
    %                         phi(r) are high-low-high.
    %      ALGORITHM:
    %           1. if phi(1)<phi(0) find  x>1 such that phi(x)>phi(0) return
    %           (0,x)
    %           2. if phi(-1)<phi(0) find x<-1 such that phi(x)>phi(0) return (-x,0)
    %           3. otherwise, return £¨-1,1£©
    % (II): ALGORITHM:
    %           while (r-l)<tol
    %               ml = .618*l+.382*r,   mr = .382*l+.618*r
    %               if ml>mr then l = ml else r = mr
    %           end while
%

% Version 0.0
% Author Linting

if nargin == 3; tol = 1e-8; else; tol = varargin{1}; end

if isfield(func,'f')
    f = @(a) func.f(x0 + a*d) - func.f(x0);
else
    error('LineSearch InputError: No Function is given!');
end

%% Find the initial search interval
iter = 0;
x  =1; 
if f(x) < 0
    while f(x) <0; x = x*2; iter = iter+1; end
    l = 0; r = x;
elseif f(-x)<0
    while f(-x)<0; x = x*2;  iter = iter+1; end
    l = -x; r = 0;
else
    l = -x; r = x;
end
%% Apply .618 method
while r - l >tol
    m1 = .382*r + .618*l; m2 = .618*r + .382*l;
    if f(m1)>f(m2); l = m1; else; r = m2; end
     iter = iter+1; 
end
alpha = (r+l)/2;


info.success = 1; % if phi is coerocive the method must sucess, otherwise the loop wont end.
info.reason =  " Success after "+ string(iter) + " iteration: " ;
info.iter = iter;
info.x = x0 + alpha * d;
info.fval = func.f(info.x);
