function [U, V, info] = nmf_admm_euc(X, r, opt)
% Alternative Direction Multiplier method (ADMM) for NMF with metic
% euclidean
%  ref:
%  DongjinSong?,DavidA.Meyer?,MartinRenqiangMin?,  "FastNonnegativeMatrixFactorizationwith Rank-oneADMM
%
%
%

[m,n] = size(X);
if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
if ~isfield(opt, 'tol'); tol = 1e-5; else; tol = opt.tol; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
if ~isfield(opt, 'rho'); rho = 5; else; rho = opt.rho; end

% random initialization
if isfield(opt,'init')
    U = opt.init.W;
    V = opt.init.H;
else
    U = rand(m,r); V = rand(r,n);
end
epoch = 0;
loss = zeros(1, maxiter);
t0   = cputime;
S = U;
T = V;
Lambda = zeros(size(U));
Pi = zeros(size(V));

for epoch = 1:maxiter
  U = (X*V'+rho*S-Lambda)/(V*V'+rho*eye(r));
  V = (U'*U+rho*eye(r))\(U'*X+rho*T-Pi);
  S = U + Lambda/rho;
  S(S<eps) = eps;
  T = V + Pi/rho;
  T(T<eps) = eps;
  Lambda = Lambda + rho*(U-S);
  Pi = Pi + rho*(V-T);
  loss(epoch) = metric_euc(X,U,V);
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