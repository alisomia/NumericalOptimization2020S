function [W,H,info] = nmf_hals(V,r,opt)
% Hierarchial Alternative Least Square (HALS)  for non-negative matrix factorization (NMF).
%
%   Ref:
%                           Andrzej Cichocki and PHAN Anh-Huy,
%                           "Fast local algorithms for large scale nonnegative matrix and tensor factorizations,"
%                           IEICE Transactions on Fundamentals of Electronics, Communications and Computer Sciences, 
%                           vol. 92, no. 3, pp. 708-721, 2009.
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
             % update H
            VtW = V'*W;
            WtW = W'*W;            
            for k=1:r
                tmp = (VtW(:,k)' - (WtW(:,k)' * H) + (WtW(k,k) * H(k,:))) / WtW(k,k);
                tmp(tmp<=eps) = eps;
                H(k,:) = tmp;
            end 
            
            % update W
            VHt = V*H';
            HHt = H*H';
            for k=1:r
                tmp = (VHt(:,k) - (W * HHt(:,k)) + (W(:,k) * HHt(k,k))) / HHt(k,k);
                tmp(tmp<=eps) = eps;
                W(:,k) = tmp;
            end

            loss(epoch) = metric_euc(V,W,H);
end
info.name = 'HALS';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch);
info.epoch = epoch;
if opt.print
    info %#ok
end
end