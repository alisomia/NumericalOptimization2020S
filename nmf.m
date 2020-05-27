% random test;
if ~exist('V'); V =  rand(100,5)*rand(5,100); end %#ok
r = 5;
opt.metric = "euc";
opt.eps = eps;
opt.maxiter =300;
opt.admm_iter = 10;
opt.tol =1e-8;
opt.rho=0.1;
opt.mu = 2;
opt.print = 1;
[W1,H1,~] = nmf_ao_admm_euc(V(1:50,1:50),r,opt);
[W2,H2,~] = nmf_ao_admm_euc(V(51:100,51:100),r,opt);
W = [W1 ;W2]; H = [H1 H2];
 %[W,H,info] = nmf_lm_admm2(V,r,opt);
% [W,H,infoo] = nmf_ao_admm_euc(VV,r,opt);
%[W,H,info1] = nmf_lm_admm2(yale64/256,r_yale,opt);
%[W,H,info2] = nmf_mu(yale64/256, r_yale, opt);
%draw(info1,info2);

%opt = rmfield(opt,'init');