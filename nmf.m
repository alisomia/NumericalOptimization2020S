% random test;
if ~exist('V'); V =  rand(100,5)*rand(5,100); end %#ok
r = 5;
opt.metric = "KL";
opt.eps = 1e-8;
opt.maxiter = 1000;
opt.rho=1;
opt.print = 1;
[W,H,info] = nmf_ao_admm_kl(V,r,opt);

