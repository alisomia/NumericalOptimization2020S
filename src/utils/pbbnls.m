function [X,gradX] = pbbnls(AtB, AtA, X, gradX, rho, gamma, M, lambda_min, lambda_max, eps, iter)
% Projected BB method of Nonlinear Least Square
% 
%   min_{X >  0} ||B-AX||_F^2
%

% initialization
Q = @(X) -sum(sum(AtB.*X)) + 0.5*sum(sum(X.*(AtA*X)));
lambda= 1/max(max(gradX));
Qlist = -inf(1,M);
Qlist(1) = Q(X);
for it = 1:iter
    if(norm(gradX(gradX<0 | X>0), 'fro') < rho)
        break;
    end
    
    alpha = 1;
    Xplus = X - lambda*gradX;
    Xplus(Xplus<eps) = eps;
    Qplus = Q(Xplus);
    D = Xplus -X;
    for i = 1:20
        if(Qplus < max(Qlist) + gamma * alpha * sum(sum(gradX.* D)))
            break
        end
        alpha = alpha/4;
        Xplus = X - alpha*lambda*gradX;
        Xplus(Xplus<eps) = eps;
        Qplus = Q(Xplus);
    end
    
    % compute BB stepsize
    s = Xplus - X;
    gradXplus = AtA*X - AtB;
    y = gradXplus - gradX;
    sinnery = sum(sum(s.*y));
    if sinnery< eps
        lambda = lambda_max;
    else
        sinners = sum(sum(s.*s));
        lambda = min(lambda_max, max(lambda_min, sinners/sinnery));
    end
    gradX = gradXplus;
    X= Xplus;
    Qlist(mod(it,M)+1) = Qplus;
end
end
        
        
