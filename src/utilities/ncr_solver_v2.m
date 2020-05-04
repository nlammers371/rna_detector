function dydt = ncr_solver(t,y,p_vec,Q)
  
    % rate vector
    v_vec = zeros(size(Q,2),1);
    % apply mass action
    for i = 1:numel(v_vec)
        v_vec(i) = p_vec(i)*prod(y(Q(:,i)<0));
    end

    % update vector
    dydt = Q*v_vec;   
    
%     if t > 100 && dydt(3) > 1
%         error('af')
%     end
end