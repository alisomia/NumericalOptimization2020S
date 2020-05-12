function s = metric_euc(V,W,H)
    %s = sum(sum((V-W*H).^2));
    s = norm(V-W*H,'fro');
end