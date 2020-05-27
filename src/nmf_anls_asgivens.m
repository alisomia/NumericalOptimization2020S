function [W,H,info] = nmf_anls_asgivens(V,r, opt)
% Alternative Projected Barzilai Borwein method (APBB) for non-negative matrix factorization (NMF).
%
%
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
if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
if ~isfield(opt,'tol'); tol = 1e-5; else; tol = opt.tol; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end

[m,n] = size(V);
if isfield(opt,'init')
    W = opt.init.W;
    H = opt.init.H;
else
    W = rand(m,r); H = rand(r,n);
end
epoch = 0;
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end

loss = zeros(1, maxiter);
t0   = cputime;

for epoch = 1:maxiter
    % disp(epoch);
    %     projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)]);
    %     if projnorm < tol_grad_ratio*init_grad
    %         break;
    %     end
    ow = 0;
    WtV = W' * V;
    for i=1:size(H,2)
        H(:,i) = nnls1_asgivens(W'*W, WtV(:,i), ow, 1, H(:,i));
    end
    
    HAt = H*V';
    Wt = W';
    for i=1:size(W,1)
        Wt(:,i) = nnls1_asgivens(H*H', HAt(:,i), ow, 1, Wt(:,i));
    end
    W = Wt';
    loss(epoch) = metric_euc(V,W,H);
    if epoch>1 && abs(loss(epoch)-loss(epoch-1))<tol*loss(epoch-1)
        break;
    end
end
info.name = 'ANLS active set givens';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch-1);
info.epoch = epoch;
if opt.print
    info %#ok
end
end