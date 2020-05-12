% random test;
r = 5;
opt.metric = "KL";
opt.eps = 1e-8;
opt.maxiter = 1000;
opt.print = 1;
[W,H,info_halsacc] = nmf_halsacc(V,r,opt);
[W,H,info_hals] = nmf_hals(V,r,opt);
draw(info_halsacc, info_hals);
