function projnorm = pgradtest_euc(V,W,H)
        gradW = W*(H*H') - V*H'; 
        gradH = (W'*W)*H - W'*V; 
        projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)],'fro');
end