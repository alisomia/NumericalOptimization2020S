function [W,H,info] = nmf_als(V,r,opt)
% Alternative Least Square (ALS) for non-negative matrix factorization (NMF).
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
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
W = rand(m,r); H = rand(r,n);
epoch = 0;
loss = zeros(1, maxiter);
t0   = cputime;
for epoch = 1:maxiter
    H = (W'*W) \ W' * V;        % H = inv(W'*W) * W' * V;
    H = H .* (H>0);
    W = V * H' / (H*H');        % W = V * H' * inv(H*H');
    W = (W>0) .* W;
    W = W ./ (repmat(sum(W), m, 1)+eps);
    loss(epoch) = metric_euc(V,W,H);
end
info.name = 'ALS';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch);
info.epoch = epoch;
if opt.print
    info %#ok
end
end