% random test;
r = 5;
opt.metric = "KL";
opt.eps = 1e-8;
opt.maxiter = 1000;
opt.print = 1;
[W,H,info_pgd] = nmf_admm_euc(V,r,opt);

