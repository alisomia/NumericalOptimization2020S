function [W,H,info] = nmf_muacc(V,r, opt)
%                   : mu_acc: Accelerated multiplicative updates (Accelerated MU)
%                       Reference:
%                           N. Gillis and F. Glineur,
%                           "Accelerated Multiplicative Updates and hierarchical ALS Algorithms for Nonnegative
%                           Matrix Factorization,",
%Neural Computation 24 (4), pp. 1085-1105, 2012.

[m,n] = size(V);
if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
if ~isfield(opt,'tol'); tol = 1e-5; else; tol = opt.tol; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
if ~isfield(opt,'alpha'); alpha = 2; else; alpha = opt.alpha; end
if ~isfield(opt,'delta'); delta = 0.1; else; delta = opt.delta; end
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
K = m*n;
rhoW = 1+(K+n*r)/(m*(r+1));
rhoH = 1+(K+m*r)/(n*(r+1));

% Update of W
for epoch=1:maxiter
    gamma = 1;
    eps0 = 1;
    j = 1;
    VHt = V * H';
    HHt = H * H.';
    while j <= floor(1+rhoW*alpha) && gamma >= delta*eps0
        W0 = W;
        W = max(eps, W .* (VHt ./ (W * HHt)));
        if j == 1;eps0 = norm(W0-W, 'fro');  end
        gamma = norm(W0-W, 'fro');
        j = j+1;
    end
    % Update of H
    gamma = 1;
    eps0 = 1;
    j = 1;
    WtV = W'*V;
    WtW = W.' * W;
    while j <= floor(1+rhoH*alpha) &&  gamma >= delta*eps0
        H0 = H;
        H = max(eps, H .* (WtV ./ (WtW * H)));
        if j == 1; eps0 = norm(H0-H, 'fro'); end
        gamma = norm(H0-H, 'fro');
        j = j+1;
    end
    loss(epoch) = metric_euc(V,W,H);
    if epoch>1 && abs(loss(epoch)-loss(epoch-1))<tol*loss(epoch-1)
        break;
    end
end

info.name = 'Accelerated MU';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch);
info.epoch = epoch;
if opt.print
    info %#ok
end
end