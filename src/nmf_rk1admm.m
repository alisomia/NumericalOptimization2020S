function [U,V,info] = nmf_rk1admm(X,r,opt)
% Alternative Direction Multiplier method (ADMM) for NMF with metic
% euclidean
%  ref:
%  DongjinSong?,DavidA.Meyer?,MartinRenqiangMin?,  "FastNonnegativeMatrixFactorizationwith Rank-oneADMM
%
%
%
[m,n] = size(X);
U = zeros(m,r);
V = zeros(r,n);
t0 = cputime;
pprint = opt.print;
opt.print = 0;
for riter = 1:r
    [u,v,~] = nmf_admm_euc(X,1,opt);
    X = X - u*v;
    U(:,riter) = u;
    V(riter,:) = v;
end
info.name = 'Rank 1 ADMM';
info.time = cputime - t0;
info.fvalue = metric_euc(X,U,V);
if pprint
    info %#ok
end
end