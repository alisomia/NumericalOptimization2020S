% we test the synthesis data, V1---V5
opt.metric = "euc";
opt.eps = eps;
opt.maxiter =300;
opt.admm_iter = 10;
opt.tol =0;
opt.rho=0.1;
opt.mu = 2;
opt.print = 1;
opt.print=0;

func = {@nmf_mu, @nmf_mumod, @nmf_muacc, @nmf_als, ...
    @nmf_hals, @nmf_halsacc, @nmf_apbb,  @nmf_alspgd, ...
    @nmf_anls_activeset, @nmf_anls_asgivens, @nmf_anls_blockpivot,...
    @nmf_admm_euc, @nmf_ao_admm_euc, @nmf_lm_admm};
tic
for i = 1:size(func,2)
    tic
[~,~,info] = feval(func{i}, yale64,5,opt);
toc
RESULT(i,1) = info.loss(10);
RESULT(i,2) = info.loss(50);
RESULT(i,3) = info.loss(100);
RESULT(i,4) = info.loss(info.epoch-1);
RESULT(i,5) = info.time;
end
toc