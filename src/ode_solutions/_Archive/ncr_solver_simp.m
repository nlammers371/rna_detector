function dydt = ncr_solver_simp(t,y,Kd,kcat,ns_rat,Q)
  
    % bound A
    bound_A = y(1) * y(3)/(Kd + y(3) + y(5));
    % rate vector
    v_vec = zeros(size(Q,2),1);
    v_vec(1) = y(2)*kcat*(bound_A + kcat); % AI
    v_vec(2) = 0;%y(6)*(ns_rat*kcat*bound_A + kcat*y(4)); %QI
    v_vec(3) = y(7)*kcat*(bound_A + y(4)); %QI
    
%     if t > 2
%         error('asfa')
%     end
            
    % update vector
    dydt = Q*v_vec;
    
end