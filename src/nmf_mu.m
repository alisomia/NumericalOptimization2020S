function [W,H,info] = nmf_mu(V, r, opt)
% Multiplicative upates (MU) for non-negative matrix factorization (NMF).
%
% Ref:
%   Daniel D. Lee and H. Sebastian Seung,
%       "Algorithms for non-negative matrix factorization," NIPS 2000.
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

% Log
% 2020 05 12 Init

[m,n] = size(V);
if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
if ~isfield(opt,'tol'); tol = 1e-5; else; tol = opt.tol; end
if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end

% random initialization
if isfield(opt,'init')
    W = opt.init.W;
    H = opt.init.H;
else
    W = rand(m,r); H = rand(r,n);
end
epoch = 0;
loss = zeros(1, maxiter);
t0   = cputime;
if metric == "euc"
    
    for epoch = 1:maxiter
        % update H
        H = H .* (W' * V) ./ (W' * W * H);
        H = H + (H<eps) .* eps;
        
        % update W
        W = W .* (V * H') ./ (W * (H * H'));
        W = W + (W<eps) .* eps;
        loss(epoch) = metric_euc(V,W,H);
        if epoch>1 && abs(loss(epoch)-loss(epoch-1))<tol*loss(epoch-1)
            break;
        end
    end
    
elseif metric == "KL"
    for epoch = 1: maxiter
        W = W.*((V./(W*H+eps))*H')./repmat(sum(H,2)'+eps,size(V,1),1);
        H = H.*(W'*(V./(W*H+eps)))./repmat(sum(W,1)'+eps,1,size(V,2));
        loss(epoch) = sum(sum(-V.*(log((W*H+eps)./(V+eps))+1)+W*H));
        if epoch>1 && abs(loss(epoch)-loss(epoch-1))<tol*loss(epoch-1)
            break;
        end
    end
else
    fprintf('No metric matched!');
end
info.name = 'MU';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch);
info.epoch = epoch;
if opt.print
    info %#ok
end
end