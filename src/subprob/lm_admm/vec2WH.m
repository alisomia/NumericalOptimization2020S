function [WW,HH] = vec2WH(vec,m,n,r)
    WW = vec(1:m*r);
    HH  = vec(m*r+1:m*r+r*n);
    WW = reshape(WW,[m,r]);
    HH = reshape(HH,[r,n]);
end