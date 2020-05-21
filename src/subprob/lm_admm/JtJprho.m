function [resvec] = JtJprho(tildevec, W, H, rho,m,n,r)
    [tildeW, tildeH] = vec2WH(tildevec,m,n,r);
    resW = tildeW * (H*H') + W * (tildeH*H') + rho*tildeW;
    resH = (W'*tildeW)*H + (W'*W) * tildeH+ rho*tildeH;
    resvec = WH2vec(resW, resH,m,n,r);
end