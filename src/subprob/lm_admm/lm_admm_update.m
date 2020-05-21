function [W,H, Waux, Haux,alphaW, alphaH] = lm_admm_update(W, H, Waux, Haux, d, W0, H0, alphaW, alphaH, nu, rho, admm_iter, m, n, r, eps, mu)

 Jtd_W =  d*H0';
    Jtd_H = W0'*d;
for i = 1:admm_iter
    % update aux var
    rhsW = Jtd_W + rho*(W - W0 )+ alphaW;
    rhsH = Jtd_H + rho*(H- H0) + alphaH;
    vec = pcg(@(tildevec) JtJprho(tildevec, W0, H0, rho,m,n,r), WH2vec(rhsW, rhsH,m,n,r), 1e-6, 1000 ,[] , [], WH2vec(Waux, Haux,m,n,r));
    %vec = minres(@(tildevec) JtJprho(tildevec, W0, H0, rho,m,n,r), WH2vec(rhsW, rhsH,m,n,r), 1e-6, 10000);% ,[] , [], WH2vec(Waux, Haux,m,n,r));
    %vec = gmres(@(tildevec) JtJprho(tildevec, W0, H0, rho,m,n,r),WH2vec(rhsW, rhsH,m,n,r), 1000, [], 1000);
    [Wadd, Hadd] = vec2WH(vec,m,n,r);
    Waux = W0 + Wadd;
    Haux = H0 + Hadd;
    %Waux(Waux<eps) = eps;
    %Haux(Haux<eps) = eps;
    
    % update primal var
    W = (nu*W0 + rho*Waux - alphaW)./(nu+rho);
    H = (nu*H0 + rho*Haux - alphaH)./(nu+rho);
    
    % project the primal var
    W(W<eps) = eps;
    H(H<eps) = eps;
    
    
    tolW = norm(W-Waux,'fro')/norm(W,'fro');
    tolH = norm(H - Haux,'fro')/norm(H,'fro');
    if tolW<1e-5&&tolH<1e-5; break; end
    % update multiplier
    alphaW = alphaW + mu*rho*(W - Waux);
    alphaH = alphaH + mu*rho*(H - Haux);

% 
%     % update aux var
%     rhsW = -Jtd_W + rho*(W0 - W+ alphaW);
%     rhsH = -Jtd_H + rho*(H0- H + alphaH);
%     vec = pcg(@(tildevec) JtJprho(tildevec, W0, H0, rho + mu,m,n,r), WH2vec(rhsW, rhsH,m,n,r), 1e-6, 1000 ,[] , [], WH2vec(Waux, Haux,m,n,r));
%     %vec = minres(@(tildevec) JtJprho(tildevec, W0, H0, rho,m,n,r), WH2vec(rhsW, rhsH,m,n,r), 1e-6, 10000);% ,[] , [], WH2vec(Waux, Haux,m,n,r));
%     %vec = gmres(@(tildevec) JtJprho(tildevec, W0, H0, rho,m,n,r),WH2vec(rhsW, rhsH,m,n,r), 1000, [], 1000);
%     [Wadd, Hadd] = vec2WH(vec,m,n,r);
%     Waux = W0 - Wadd;
%     Haux = H0 - Hadd;
%     %Waux(Waux<eps) = eps;
%     %Haux(Haux<eps) = eps;
%     
%     % update primal var
%     %W = (nu*W0 + rho*Waux - alphaW)./(nu+rho);
%     %H = (nu*H0 + rho*Haux - alphaH)./(nu+rho);
%     W = Waux + alphaW;
%     H = Haux + alphaH;
%     % project the primal var
%     W(W<eps) = eps;
%     H(H<eps) = eps;
%     
%     % update multiplier
%     alphaW = alphaW - mu*(W - Waux);
%     alphaH = alphaH - mu*(H - Haux);
% 






end
end