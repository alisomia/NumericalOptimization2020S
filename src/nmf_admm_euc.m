function [W, H, info] = nmf_admm_euc(V, r, opt)
% Alternative Direction Multiplier method (ADMM) for NMF with metic
% euclidean
% % % %  ref:
%  DongjinSong?,DavidA.Meyer?,MartinRenqiangMin?,  "FastNonnegativeMatrixFactorizationwith Rank-oneADMM
%
%
%

[m,n] = size(V);
if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
if ~isfield(opt, 'tol'); tol = 1e-5; else; tol = opt.tol; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
if ~isfield(opt, 'rho'); rho = 5; else; rho = opt.rho; end

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
S = W;
T = H;
Lambda = zeros(size(W));
Pi = zeros(size(H));

for epoch = 1:maxiter
    W = (V*H'+rho*S-Lambda)/(H*H'+rho*eye(r));
    H = (W'*W+rho*eye(r))\(W'*V+rho*T-Pi);
    S = W + Lambda/rho;
    S(S<eps) = eps;
    T = H + Pi/rho;
    T(T<eps) = eps;
    Lambda = Lambda + rho*(W-S);
    Pi = Pi + rho*(H-T);
    loss(epoch) = metric_euc(V,W,H);
    if epoch>1 && abs(loss(epoch)-loss(epoch-1))<tol*loss(epoch-1)
        break;
    end
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