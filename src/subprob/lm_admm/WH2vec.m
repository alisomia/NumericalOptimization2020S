function [vec] = WH2vec(WW, HH,m,n,r)
vec = [reshape(WW,[m*r,1]); reshape(HH,[r*n,1])];
end