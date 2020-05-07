function f_trend_interp = ode_objfunction_v2(t_max,rate_params,rate_fun,stoich_mat,y0_cell,fluo_indices,time_exp_array)

    % generate full rate vector
    rate_params = num2cell(rate_params);
    rate_vec = rate_fun(rate_params{:});
    
    % enforce proper shape
    time_exp_array = reshape(time_exp_array,[],numel(y0_cell));
    
    % initialize array to store data
    f_trend_interp = NaN(size(time_exp_array));
    
    % iterate through y0 data and solve ODE for each case
    for i = 1:numel(y0_cell)
        y0_vec = y0_cell{i};
        
        % solve ODE numerically
        [t_sol,y_sol]=ode15s(@(t,y) ncr_solver_v2(t,y,rate_vec,stoich_mat), [0 t_max], y0_vec);

        % get fluorescence trend
        f_trend_raw = sum(y_sol(:,fluo_indices),2);
        f_trend_interp(:,i) = interp1(t_sol,f_trend_raw,time_exp_array(:,i));
    end
end
