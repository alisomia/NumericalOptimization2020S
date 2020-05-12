function [W,H,info] = nmf_mumod(V, r, opt)
%Modified multiplicative upates (MU) for NMF
%
% Ref:
%                           C.-J. Lin,
%                           "On the convergence of multiplicative update algorithms for nonnegative matrix factorization,"
%                           IEEE Trans. Neural Netw. vol.18, no.6, pp.1589?1596, 2007. 
%
[m,n] = size(V);
if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end

% random initialization
W = rand(m,r); H = rand(r,n);
epoch = 0;
loss = zeros(1, maxiter);
t0   = cputime;
delta = eps;
for epoch = 1:maxiter
            WtW = W' * W;
            WtV = W' * V;
            gradH = WtW * H - WtV;
            Hb = max(H, (gradH < 0)* eps);
            H = H - Hb ./ (WtW * Hb + delta) .* gradH;
            
            % update W
            HHt = H * H';
            VHt = V * H';
            gradW = W * HHt - VHt;
            Wb = max(W, (gradW < 0)* eps);
            W = W - Wb ./ (Wb * HHt + delta) .* gradW;
            
            S = sum(W,1);
            W = W ./ repmat(S,m,1);
            H = H .* repmat(S',1,n);
            loss(epoch) = metric_euc(V,W,H);
end
info.name = 'modified MU';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch);
info.epoch = epoch;
if opt.print
    info %#ok
end
end