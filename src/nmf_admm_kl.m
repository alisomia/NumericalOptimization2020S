function [W, H, info] = nmf_admm_kl(V, r, opt)
%   D.L. Sun and C. F?votte, "Alternating direction method of multipliers
%    for non-negative matrix factorization with the beta divergence",
%      ICASSP 2014.

[m,n] = size(V);
if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
if ~isfield(opt, 'rho'); rho = 5; else; rho = opt.rho; end

% random initialization
W = rand(m,r); H = rand(r,n);
epoch = 0;
loss = zeros(1, maxiter);
t0   = cputime;
X = W*H;
Wplus = W;
Hplus = H;
alphaX = zeros(size(X));
alphaW = zeros(size(W));
alphaH = zeros(size(H));
fixed=[];
free = setdiff(1:r, fixed);
for epoch = 1:maxiter
    % update for H
    H = (W'*W + eye(r)) \ (W'*X + Hplus + 1/rho*(W'*alphaX - alphaH));
    
    % update for W
    P = H*H' + eye(r);
    Q = H*X' + Wplus' + 1/rho*(H*alphaX' - alphaW');
    W(:,free) = ( P(:,free) \ (Q - P(:,fixed)*W(:,fixed)') )';
    
    % update for X (this is the only step that depends on beta)
    X_ap = W*H;
    b = rho*X_ap - alphaX - 1;
    X = (b + sqrt(b.^2 + 4*rho*V))/(2*rho);
    
    % update for H_+ and W_+
    Hplus = max(H + 1/rho*alphaH, 0);
    Wplus = max(W + 1/rho*alphaW, 0);
    
    % update for dual variables
    alphaX = alphaX + rho*(X - X_ap);
    alphaH = alphaH + rho*(H - Hplus);
    alphaW = alphaW + rho*(W - Wplus);
    
    loss(epoch) = metric_kl(V,W,H);
end

info.name = 'ADMM';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch);
info.epoch = epoch;
if opt.print
    info %#ok
end
end