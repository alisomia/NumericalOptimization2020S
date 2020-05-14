function [W,H,info] = nmf_ao_admm_kl(V,r,opt)

if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
if ~isfield(opt,'tol'); tol = 1e-2; else tol = opt.tol; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
if ~isfield(opt, 'admm_iter'); admm_iter = 10; else; admm_iter = opt.admm_iter; end
if ~isfield(opt, 'rho'); rho = 1; else; rho = opt.rho; end

[m,n] = size(V);
W = rand(m,r); H = rand(r,n);
dual_W = zeros(size(W));
dual_H = zeros(size(H));
V_aux = zeros(size(V));
dual_V = zeros(size(V));
epoch = 0;
loss = zeros(1, maxiter);
t0 = cputime;
for epoch = 1:maxiter
    [H, dual_H, V_aux, dual_V] = admm_kl_update(V, V_aux, dual_V, W,H,dual_H, r, admm_iter, tol,eps);
    [W,dual_W, V_aux, dual_V] = admm_kl_update(V',V_aux', dual_V', H', W', dual_W', r, admm_iter, tol, eps);
    W = W';
    dual_W = dual_W';
    V_aux = V_aux';
    dual_V = dual_V';
    loss(epoch) = metric_euc(V,W,H);
end

info.name = 'AO-ADMM';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch);
info.epoch = epoch;
if opt.print
    info %#ok
end
end

