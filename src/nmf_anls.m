function [W,H, info] = nmf_anls(V, r, opt)
if isfield(opt,'method') && (opt.method=="bp")
    [W,H,info] = nmf_anls_blockpivot(V,r,opt);
elseif isfield(opt,'method') && (opt.method=="givens")
    [W,H,info] = nmf_anls_asgivens(V,r,opt);
else
    [W,H,info] = nmf_anls_activeset(V,r,opt);