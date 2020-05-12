function s = metric_kl(V,W,H)
s = sum(sum(-V.*(log((W*H+eps)./(V+eps))+1)+W*H));
end