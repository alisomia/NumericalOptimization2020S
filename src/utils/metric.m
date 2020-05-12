function s = metric(V,W,H,met)
if met=="euc"; s = metric_euc(V,W,H); else s = metric_kl(V,W,H); end