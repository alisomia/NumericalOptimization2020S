function [W,H,info] = nmf_halsacc(V,r,opt)
% Acceralted Hierarchial Alternative Least Square (acc HALS)  for non-negative matrix factorization (NMF).
%
%   Ref:
%                           N. Gillis and F. Glineur,
%                           "Accelerated Multiplicative Updates and hierarchical ALS Algorithms for Nonnegative
%                           Matrix Factorization,",
%                           Neural Computation 24 (4), pp. 1085-1105, 2012.
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
if ~isfield(opt,'tol'); tol = 1e-5; else; tol = opt.tol; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
if ~isfield(opt,'alpha'); alpha = 2; else; alpha = opt.alpha; end
if ~isfield(opt,'delta'); delta = 0.1; else; delta = opt.delta; end

if isfield(opt,'init')
    W = opt.init.W;
    H = opt.init.H;
else
    W = rand(m,r); H = rand(r,n);
end
epoch = 0;
loss = zeros(1, maxiter);
t0   = cputime;

eit1 = cputime;
VHt = V*H';
HHt = H*H';
eit1 = cputime-eit1;

scaling = sum(sum(VHt.*W))/sum(sum( HHt.*(W'*W) ));
W = W * scaling;

for epoch = 1:maxiter
    % Update of W
    if epoch > 0 % Do not recompute A and B at first pass
        % Use actual computational time instead of estimates rhoU
        eit1 = cputime;
        VHt = V*H';
        HHt = H*H';
        eit1 = cputime-eit1;
    end
    W = HALSupdt(W',HHt',VHt', eit1, alpha, delta, eps);
    W = W';
    
    % Update of H
    eit1 = cputime;
    WtV = W'*V;
    WtW = W'*W;
    eit1 = cputime-eit1;
    H = HALSupdt(H, WtW, WtV, eit1, alpha, delta, eps);
    
    loss(epoch) = metric_euc(V,W,H);
    if epoch>1 && abs(loss(epoch)-loss(epoch-1))<tol*loss(epoch-1)
        break;
    end
end
info.name = 'Accelerated HALS';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch);
info.epoch = epoch;
if opt.print
    info %#ok
end
end

function V = HALSupdt(V,UtU,UtM,eit1,alpha,delta, ~)
[r, ~] = size(V);
eit2 = cputime; % Use actual computational time instead of estimates rhoU
cnt = 1; % Enter the loop at least once
eps = 1;
eps0 = 1;
eit3 = 0;
while cnt == 1 || (cputime-eit2 < (eit1+eit3)*alpha && eps >= (delta)^2*eps0)
    nodelta = 0; if cnt == 1, eit3 = cputime; end
    for k = 1 : r
        deltaV = max((UtM(k,:)-UtU(k,:)*V)/UtU(k,k),-V(k,:));
        V(k,:) = V(k,:) + deltaV;
        nodelta = nodelta + deltaV*deltaV'; % used to compute norm(V0-V,'fro')^2;
        if V(k,:) == 0, V(k,:) = eps*max(V(:)); end % safety procedure
    end
    if cnt == 1
        eps0 = nodelta;
        eit3 = cputime-eit3;
    end
    eps = nodelta; cnt = 0;
end
end
