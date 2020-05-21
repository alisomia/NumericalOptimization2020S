% random test;
if ~exist('V'); V =  rand(100,5)*rand(5,100); end %#ok
r = 5;
opt.metric = "euc";
opt.eps = 1e-8;
opt.maxiter =100;
opt.admm_iter = 20;
opt.tol =1e-8;
opt.rho=4;
opt.print = 1;
nmf_lm_admm(V,r,opt);
%[W,H,info1] = nmf_lm_admm(yale64,r_yale,opt);
%[W,H,info2] = nmf_ao_admm_euc(yale64, r_yale, opt);
draw(info1,info2);