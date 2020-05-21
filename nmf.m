% random test;
if ~exist('V'); V =  rand(100,5)*rand(5,100); end %#ok
r = 5;
opt.metric = "euc";
opt.eps = 1e-8;
opt.maxiter =100;
opt.admm_iter = 10;
opt.tol =1e-8;
opt.rho=0.1;
opt.mu = 2;
opt.print = 1;
[W,H,info] = nmf_ao_admm_euc(V,r,opt);
%[W,H,info1] = nmf_lm_admm(yale64,r_yale,opt);
%[W,H,info2] = nmf_ao_admm_euc(yale64, r_yale, opt);
draw(info1,info2);