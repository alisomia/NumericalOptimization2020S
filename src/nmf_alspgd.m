function [W,H,info] = nmf_alspgd(V,r, opt)
% Projected Gradient Descent (PGD) for non-negative matrix factorization (NMF).
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
%if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'alpha'); alpha = 1; else; alpha = opt.alpha; end
%if ~isfield(opt,'delta'); delta = 0.1; else; delta = opt.delta; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
if ~isfield(opt,'tol_grad_ratio'); tol_grad_ratio = 1e-4; else; tol_grad_ratio = opt.tol_grad_ratio; end
if isfield(opt,'init')
    W = opt.init.W;
    H = opt.init.H;
else
    W = rand(m,r); H = rand(r,n);
end
epoch = 0;
loss = zeros(1, maxiter);
t0   = cputime;
tol_grad_ratio = 0.00001; % tol = [0.001; 0.0001; 0.00001];
gradW = W*(H*H') - V*H';
gradH = (W'*W)*H - W'*V;
init_grad = norm([gradW; gradH'],'fro');
tolW = max(0.001,tol_grad_ratio)*init_grad;
tolH = tolW;

for epoch = 1:maxiter
    projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)],'fro');
    if projnorm < tol_grad_ratio*init_grad
        break;
    end
    
    [W, gradW, iterW] = nlssubprob(V', H', W', tolW, 1000);
    W = W';
    gradW = gradW';
    
    if iterW == 1
        tolW = 0.1 * tolW;
    end
    
    [H, gradH, iterH] = nlssubprob(V, W, H, tolH, 1000);
    if iterH == 1
        tolH = 0.1 * tolH;
    end
    loss(epoch) = metric_euc(V,W,H);
end
info.name = 'ALSPGD';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch-1);
info.epoch = epoch;
if opt.print
    info %#ok
end
end