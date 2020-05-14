function [h, dual_h, v_aux, dual_v] = admm_kl_update(v, v_aux, dual_v, w, h, dual_h, k, admm_iter, ~, eps)
    g = w'*w;
    rho = trace(g)/k;
     cho = chol(g+  rho*eye(size(g,1)));
     for i = 1:admm_iter
         h_aux = cho\(cho'\((w'*(v_aux+dual_v)+rho*(h+dual_h))));
        h = h_aux - dual_h;
        h(h<eps) = eps;
        
        % update v_aux
        v_bar = w*h_aux - dual_v;
        v_aux = 0.5*((v_bar - 1) + sqrt((v_bar-1).^2+4*v));
        
        % dual variable updating
        dual_h = dual_h + h - h_aux;
        dual_v = dual_v + v_aux - w*h_aux;
     end
end