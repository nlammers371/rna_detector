function dydt = ncr_solver_param(t,y,r_vec,Q,rate_params)
  
    % rate vector
    v_vec = zeros(size(Q,2),1);
    % apply mass action
    for i = 1:numel(v_vec)
        v_vec(i) = r_vec(i)*prod(y(Q(:,i)<0));
    end

    % update vector
    dydt = Q*v_vec;   
end