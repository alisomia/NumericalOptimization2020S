function [W,H,info] = nmf_ao_admm_euc(V,r,opt)
% use ALS to solve the least square problem, 
% the subproblem is solved by (projected) ADMM
%
%   
%
% INPUT
% --------------------------------------------------------
% V [m*n double]: the original nonnegative matrix, need to factorize
% r [int]: the rank of NMF
% opt [struct]: additional options
%   - admm_iter [int]: admm iteration, default = 10
%   - rho, mu [double]: admm param, defult = 5,1
%   - eps [double]: precision, default = 1e-8
%   - maxiter [int] : max iteration, default = 1000
% OUTPUT
% ----------------------------------------------------------
% W [m*r double], V [r*n double]: the factor
% info [struct]
%   - fvalue [double]: fvalue
%   - epoch [int]: iteration time
%   - loss [array]: the loss at each epoch
%
% Author: Ting Lin @ PKU

if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
if ~isfield(opt,'tol'); tol = 1e-2; else; tol = opt.tol; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
if ~isfield(opt, 'admm_iter'); admm_iter = 10; else; admm_iter = opt.admm_iter; end


[m,n] = size(V);
if isfield(opt,'init')
    W = opt.init.W;
    H = opt.init.H;
else
    W = rand(m,r); H = rand(r,n);
end
dual_W = zeros(size(W));
dual_H = zeros(size(H));
epoch = 0;
loss = zeros(1, maxiter);
t0 = cputime;
for epoch = 1:maxiter
    [H, dual_H] = admm_ls_update(V,W,H,dual_H, r, admm_iter, tol,eps);
    [W,dual_W] = admm_ls_update(V', H', W', dual_W', r, admm_iter, tol, eps);
    W = W';
    dual_W = dual_W';
    loss(epoch) = metric_euc(V,W,H);
    if epoch>1 && abs(loss(epoch)-loss(epoch-1))<tol*loss(epoch-1)
        break;
    end
end

info.name = 'AO-ADMM';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch);
info.epoch = epoch;
if opt.print
    info %#ok
end
end





