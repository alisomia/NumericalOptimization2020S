function gval = numgrad(f, x)
%   A simple file implements the numerical gradient.
%   
% USAGE
% -----------------------------------------------
%   gval = numgrad(f,x)
%
% INPUT
% ------------------------------------------------
% f [function vector -> double] : the target (scalar) function
% f [function vector -> vector]* : the target (vector) function
% x [vector] : the point which the gradient is evaluated.
% 
% OUTPUT
%------------------------------------------------
% gval [vector]: the grad f(x)
%       if f is [function vector -> double]
% gval [matrix]* : the Jacobian of f(x)
%       if f is [function vector -> vector]
%
% Latest Update: 2020/04/01
% author: Linting
    dx = zeros(size(x,1),1); 
    gval =zeros(size(x,1),size(f(x)',2)); 
    for idx = 1:size(x,1)
        dx(idx)= 1e-7*max(abs(x));
        gval(idx,:) = (f(x+dx) - f(x-dx))'./(2*dx(idx));
        dx(idx) = 0;
    end
end