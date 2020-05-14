function [W,H,info] = nmf_anls_blockpivot(V,r, opt)
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
if isfield(opt,'init')
    W = opt.init.W;
    H = opt.init.H;
else
    W = rand(m,r); H = rand(r,n);
end
epoch = 0;
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end

loss = zeros(1, maxiter);
t0   = cputime;

for epoch = 1:maxiter
    H = nnlsm_blockpivot(W'*W, W'*V, 1, H);
    W = nnlsm_blockpivot(H*H', H*V', 1, W');
    W = W';
    loss(epoch) = metric_euc(V,W,H);
    if epoch>1 && abs(loss(epoch)-loss(epoch-1))<tol*loss(epoch-1)
        break;
    end
end
info.name = 'ANLS active set givens';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch-1);
info.epoch = epoch;
if opt.print
    info %#ok
end
end