function [h,dual] = admm_ls_update(y,w,h,dual, k, admm_iter,~,  eps)
    g = w'*w;
    rho = trace(g)/k;
    cho = chol(g+  rho*eye(size(g,1)));
    wty = w'*y;
    for i = 1:admm_iter
        h_aux = cho\(cho'\(wty+rho*(h+dual)));
        h_prev = h;
        h = h_aux - dual;
        h(h<eps) = eps;
        dual = dual + h - h_aux;
%         if terminate(h, h_prev, h_aux, dual, tol)
%             break;
%         end
    end
end
