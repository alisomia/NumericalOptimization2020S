function [W,H,info] = nmf_fpa(V,r,opt)

%if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
%if ~isfield(opt,'tol'); tol = 1e-2; else tol = opt.tol; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
if ~isfield(opt, 'maxsubiter'); maxsubiter = 10; else; maxsubiter= opt.maxsubiter; end
%if ~isfield(opt, 'rho'); rho = 1; else; rho = opt.rho; end

t0   = cputime;
[m,n] = size(V);
W = rand(m,r); H = rand(r,n);
loss = zeros(1,maxiter);

% Set parameters
chi   = -V./(W*H);
chi   = bsxfun(@times, chi, 1./max(bsxfun(@times, -W'*chi, 1./sum(W,1)')));
Wbar  = W;
Wold  = W;
Hbar  = H;
Hold  = H;


for epoch = 1:maxiter
    
    % Computation of H
    sigma = sqrt(n/r) * sum(W(:)) ./ sum(V,1)  / norm(W, 'fro');
    tau   = sqrt(r/n) * sum(V,1)  ./ sum(W(:)) / norm(W, 'fro');

    for j = maxsubiter
        
        chi  = chi + bsxfun(@times, W*Hbar, sigma);
        chi  = (chi - sqrt(chi.^2 + bsxfun(@times, V, 4*sigma)))/2;
        H    = max(H - bsxfun(@times, W'*(chi + 1), tau), 0);
        Hbar = 2*H - Hold;
        Hold = H;
        
    end
    
    % Computation of W
    sigma = sqrt(m/r) * sum(H(:)) ./ sum(V,2)  / norm(H,'fro');
    tau   = sqrt(r/m) * sum(V,2)  ./ sum(H(:)) / norm(H,'fro');

    for j = 1:maxsubiter
        
        chi  = chi + bsxfun(@times, Wbar*H, sigma);
        chi  = (chi - sqrt(chi.^2 + bsxfun(@times, V, 4*sigma)))/2;
        W    = max(W - bsxfun(@times, (chi + 1)*H', tau), 0);
        Wbar = 2*W - Wold;
        Wold = W;
        
    end
    
    loss(epoch) = metric_kl(V,W,H);
        
end


info.name = 'FPA';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch);
info.epoch = epoch;
if opt.print
    info %#ok
end
end

