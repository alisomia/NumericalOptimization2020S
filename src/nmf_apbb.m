function [W,H,info] = nmf_apbb(V,r, opt)
% Alternative Projected Barzilai Borwein method (APBB) for non-negative matrix factorization (NMF).
%   
%
% INPUT
% --------------------------------------------------------
% V [m*n double]: the original nonnegative matrix, need to factorize
% r [int]: the rank of NMF
% opt [struct]: additional options
%   - metric [string]: "euc" or "KL"
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

[m,n] = size(V);
if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
%if ~isfield(opt,'alpha'); alpha = 1; else; alpha = opt.alpha; end
%if ~isfield(opt,'delta'); delta = 0.1; else; delta = opt.delta; end
if ~isfield(opt,'gamma'); gamma = 1e-4; else; gamma = opt.gamma; end
if ~isfield(opt, 'rho'); rho = 1e-4; else; rho = opt.rho;end
if ~isfield(opt,'M'); M = 5; else; M = opt.M; end
if ~isfield(opt, 'lambda_min'); lambda_min = 1e-4; else; lambda_min = opt.lambda_min; end
if ~isfield(opt, 'lambda_max'); lambda_max = 1e4; else; lambda_max = opt.lambda_max; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
if ~isfield(opt,'tol_grad_ratio'); tol_grad_ratio = 1e-4; else; tol_grad_ratio = opt.tol_grad_ratio; end
W = rand(m,r); H = rand(r,n);
epoch = 0;
loss = zeros(1, maxiter);
t0   = cputime;
gradW = W*(H*H') - V*H';
gradH = (W'*W)*H - W'*V;
init_grad = norm([gradW; gradH'],'fro');

for epoch = 1:maxiter
    projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)]);
    if projnorm < tol_grad_ratio*init_grad
        break;
    end
    
    [W, gradW] = pbbnls(H*V', H*H', W', gradW',rho, gamma, M, lambda_min, lambda_max, eps, 20);
    W = W';
    gradW = gradW';

    
    [H, gradH] = pbbnls(W'*V, W'*W, H, gradH, rho, gamma, M, lambda_min, lambda_max, eps, 20);

    loss(epoch) = metric_euc(V,W,H);
end
info.name = 'APBB';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch-1);
info.epoch = epoch;
if opt.print
    info %#ok
end
end