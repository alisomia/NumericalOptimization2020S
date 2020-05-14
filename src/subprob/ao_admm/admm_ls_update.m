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
%         if terminate(h, h_prev, h_aux, dual, 1e-2)
%             break;
%         end
    end
end

function b = terminate(mat, mat_prev, aux, dual, tol)
r = norm(mat - aux,'fro')/norm(mat,'fro');
s= norm(mat - mat_prev, 'fro')/norm(dual);
b = (r<tol && s<tol);
end
