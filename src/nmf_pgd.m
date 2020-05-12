function [W,H,info] = nmf_pgd(V,r, opt)
% Projected Gradient Descent (PGD) for non-negative matrix factorization (NMF).
%
%
% INPUT
% --------------------------------------------------------
% V [m*n double]: the original nonnegative matrix, need to factorize
% r [int]: the rank of NMF
% opt [struct]: additional options
%   - metric [string]: "euc" or "KL"
%   - eps [double]: precision, default = 1e-8
%   - maxiter [int] : max iteration, default = 1000
% OUTPUT
% ----------------------------------------------------------
% W [m*r double], V [r*n double]: the factor
% info [struct]
%   - fvalue [double]: fvalue
%   - epoch [int]: iteration time
%   - loss [array]: the loss at each epoch
%
% Author: Ting Lin @ PKU

[m,n] = size(V);
%if ~isfield(opt,'eps'); eps = 1e-6; else; eps = opt.eps; end
%if ~isfield(opt,'metric');metric = "KL"; else; metric = opt.metric; end
if ~isfield(opt,'alpha'); alpha = 1; else; alpha = opt.alpha; end
%if ~isfield(opt,'delta'); delta = 0.1; else; delta = opt.delta; end
if ~isfield(opt,'maxiter'); maxiter = 1000; else; maxiter = opt.maxiter; end
if ~isfield(opt,'tol_grad_ratio'); tol_grad_ratio = 1e-4; else; tol_grad_ratio = opt.tol_grad_ratio; end
W = rand(m,r); H = rand(r,n);
epoch = 0;
loss = zeros(1, maxiter);
t0   = cputime;
        gradW = W*(H*H') - V*H'; 
        gradH = (W'*W)*H - W'*V;           
        init_grad = norm([gradW; gradH'],'fro');    
        H = nlssubprob(V,W,H,0.001,1000);    
        obj = 0.5*(norm(V-W*H,'fro')^2);
        
for epoch = 1:maxiter
 gradW = W*(H*H') - V*H';
            gradH = (W'*W)*H - W'*V;

            projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)]);  
            if projnorm < tol_grad_ratio*init_grad
                fprintf('final grad norm %f\n', projnorm);
            else
                Wn = max(W - alpha*gradW,0);    
                Hn = max(H - alpha*gradH,0);    
                newobj = 0.5*(norm(V-Wn*Hn,'fro')^2);

                if newobj-obj > 0.01*(sum(sum(gradW.*(Wn-W)))+ sum(sum(gradH.*(Hn-H))))
                    % decrease stepsize    
                    while 1
                        alpha = alpha/10;
                        Wn = max(W - alpha*gradW,0);    
                        Hn = max(H - alpha*gradH,0);    
                        newobj = 0.5*(norm(V-Wn*Hn,'fro')^2);

                        if newobj - obj <= 0.01*(sum(sum(gradW.*(Wn-W)))+ sum(sum(gradH.*(Hn-H))))
                            W = Wn; H = Hn;
                            obj = newobj;
                        break;

                        end
                    end
                else 
                    % increase stepsize
                    while 1
                        Wp = Wn; 
                        Hp = Hn; 
                        objp = newobj;
                        alpha = alpha*10;
                        Wn = max(W - alpha*gradW,0);    
                        Hn = max(H - alpha*gradH,0);    
                        newobj = 0.5*(norm(V-Wn*Hn,'fro')^2);

                        %if (newobj - obj > 0.01*(sum(sum(gradW.*(Wn-W)))+ ...
                        %    sum(sum(gradH.*(Hn-H))))) | (Wn==Wp & Hn==Hp)
                        if (newobj - obj > 0.01*(sum(sum(gradW.*(Wn-W)))+ sum(sum(gradH.*(Hn-H))))) ...
                                || (isequal(Wn, Wp) && isequal(Hn, Hp))               
                            W = Wp; 
                            H = Hp;
                            obj = objp; 
                            alpha = alpha/10;
                            break;
                        end
                    end
                end 
            end
            loss(epoch) = metric_euc(V,W,H);
end
info.name = 'PGD';
info.time = cputime - t0;
info.loss = loss;
info.fvalue = loss(epoch);
info.epoch = epoch;
if opt.print
    info %#ok
end
end