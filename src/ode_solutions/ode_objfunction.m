function f_trend_interp = ode_objfunction(t_max,rate_params,rate_fun,stoich_mat,y0_vec,fluo_indices,time_exp)
% generate full rate vector
rate_params = num2cell(rate_params);
rate_vec = rate_fun(rate_params{:});
% solve ODE numerically
[t_sol,y_sol]=ode15s(@(t,y) ncr_solver(t,y,rate_vec,stoich_mat), [0 t_max], y0_vec);

% get fluorescence trend
f_trend_raw = sum(y_sol(:,fluo_indices),2);
f_trend_interp = interp1(t_sol,f_trend_raw,time_exp);
% delta_f = sum(f_trend_interp-f_exp);
end
