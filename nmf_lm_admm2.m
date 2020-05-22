% use LM to solve the least square problem, 
% the subproblem is solved by (projected) ADMM

function [W,H,info] = nmf_lm_admm(V,r,opt)

if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
if ~isfield(opt,'tol'); tol = 1e-2; else; tol = opt.tol; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
if ~isfield(opt, 'admm_iter'); admm_iter = 10; else; admm_iter = opt.admm_iter; end
if ~isfield(opt, 'rho'); rho = 5; else; rho = opt.rho; end
if ~isfield(opt, 'mu'); mu = 1; else; mu = opt.mu; end

[m,n] = size(V);
if isfield(opt,'init')
    W = opt.init.W;
    H = opt.init.H;
else
    W = rand(m,r); H = rand(r,n);
end
Waux = W;
Haux = H;
alphaW = zeros(size(W));
alphaH = zeros(size(H));
loss = zeros(1, maxiter);
t0 = cputime;
nu = 1;
R = V - W*H;
gradW = W*(H*H') - V*H';
gradH = (W'*W)*H - W'*V;
pong = 0;
% gamma = 0;
for epoch = 1: maxiter
    if nu >1e10; epoch= epoch - 1;break; end
    [Wnew,Hnew, Wauxnew, Hauxnew,~, ~] = ...
        lm_admm_update(W, H, Waux, Haux, R, W, H, alphaW, alphaH, ...
        nu, rho, admm_iter, m, n, r, eps, mu);
    Rnew = V - Wnew*Hnew;
    dW = Wnew - W;
    dH = Hnew - H;
    %delta_g = sum(sum(dW.*(nu*dW - gradW))) + sum(sum(dH.*(nu*dH - gradH)));
    %delta_g = delta_g / 2;
    dvec = WH2vec(dW,dH,m,n,r);
    delta_g = dvec'*JtJprho(dvec, W,H, nu, m,n,r) - 2*sum(sum(gradW.*dW)) - 2*sum(sum(gradH.*dH));
    delta_f = sum(sum(R.^2)) - sum(sum(Rnew.^2));
    if delta_g<0 && delta_f <0
    [H, ~] = admm_ls_update(V,W,H,alphaH, r, admm_iter, tol,eps);
    [W,~] = admm_ls_update(V', H', W', alphaW', r, admm_iter, tol, eps);
    W = W';
    R = V - W*H;
    gradW = W*(H*H') - V*H';
    gradH = (W'*W)*H - W'*V;
    loss(epoch) = metric_euc(V,W,H);
    mu = mu*2;
    pong = pong + 1;
    continue;
    end
    
    gamma = delta_f/delta_g;
%     gamma = delta_f;
    if gamma>0.75; nu = max(1e-4,nu/2); end
    if gamma<0.25; nu = nu*4; end
    if gamma > 0 
        W = Wnew;
        H = Hnew;
        Waux = Wauxnew;
        Haux = Hauxnew;
       
        R = Rnew;
        mu = mu/2;
        gradW = W*(H*H') - V*H';
    gradH = (W'*W)*H - W'*V;
    %disp('accept');
    %else
    %    mu = mu*2;
    end
    loss(epoch) = metric_euc(V,W,H);
end

info.name = 'LM-ADMM';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch);
info.epoch = epoch;
info.pong = pong;
if opt.print
    info %#ok
end
end









