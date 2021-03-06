clear 
close all

% basic path info
addpath('../utilities')
DataPath = '../../out/ode_studies_v2/';
mkdir(DataPath)

% specify project to load
project = 'ncr_basic';
% load 
load([DataPath project '.mat'])

% make figure path
FigPath = ['../../fig/ode_studies_v2/' project '/' ];
mkdir(FigPath)


% extract reaction rate quantities that need to be specified/fit
reaction_rate_index = unique(rate_vec);

% decompose rates into fundamental units
reaction_units = [];
for i = 1:numel(reaction_rate_index)
    r_sym = reaction_rate_index(i);
    t_char = char(r_sym);
    char_list = strsplit(t_char,'*');
    for c = 1:numel(char_list)
        eval(['syms ' char_list{c}])
        reaction_units = [reaction_units eval(char_list{c})];
    end
end
reaction_unit_index = unique(reaction_units);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% set key parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
t_max = 3600; % duration of time to solve for
A10 = 1e-1;
S0 = 200;
%%%% specify initial reactant concentrations 
% vectors for substitutions
init_sym_array =     {'A1' 'S' 'IA2' 'C13' 'G1' 'G2'};
init_val_vec_cage =    [0    S0 20   20.8   20  0.8]; % cage only
init_val_vec_full =    [A10  S0 20   20.8   20  0.8]; % primary only
init_val_vec_primary = [A10  S0 0    20.8   20  0.8]; % primary only
init_val_vec_neg =     [0    S0 0    20.8   20  0.8]; % no activator

% initialize concentration vector. Everything not explicitly assigned taken
% to be zero
y0_vec_full = zeros(size(full_reactant_list));
y0_vec_primary = zeros(size(full_reactant_list));
y0_vec_cage = zeros(size(full_reactant_list));
y0_vec_neg = zeros(size(full_reactant_list));
for i = 1:numel(init_sym_array)
    species = init_sym_array{i};
    spec_indices = find(strcmp(full_reactant_list,species));    
    y0_vec_full(spec_indices) = init_val_vec_full(i);
    y0_vec_primary(spec_indices) = init_val_vec_primary(i);
    y0_vec_cage(spec_indices) = init_val_vec_cage(i);
    y0_vec_neg(spec_indices) = init_val_vec_neg(i);
end



%%%%%%%%%%%%%%%%%
%%% Specify which rate parameters to fit and specify values for "given"
%%% parameters

% specify subset of parameters to fit
fit_indices_round1 = [1:5 7];
fit_indices_round2 = 6;
n_inference_runs = 20;
sim_noise = 1e-2*S0;

% define "true" rates to fit
b = 1e-6; % Cas13 1/SNR
k = 0.025; % association rate (s^-1 nM^-1)
kc = 200; % Cas13 catalytic rate
rcg = 1e-2; % off rate for Cas13:gRNA
rga = 1e-2; % off rate for gRNA:Activator
ria = 1e-2; % off rate for CleavedCage:Activator
rns = 1e3; % off rate for all nonspecific interactions

% get set of true param values
true_param_vec = eval(reaction_unit_index);

% specify degree of uncertainty for each param value
sigma_vec = [1e3, 1e1, 1e1, 1e1, 1e1, 1e1, 1e1];

% rate vector to substitute into
rate_vec_val = rate_vec;
% generate valued long rate vec
for r = 1:numel(reaction_unit_index)
    rr = reaction_unit_index(r);        
    rate_vec_val = subs(rate_vec_val,rr,eval(rr));    
end
rate_vec_val = double(rate_vec_val);


%%%%%%%%%%%%%%% Generate "experimental data" using true param values %%%%%%%%%%%%%%%%%%
close all
f_indices = find(strcmp(full_reactant_list,{'F'}));

[t_ncr_full,y_ncr_full] = ode15s(@(t,y) ncr_solver_v2(t,y,rate_vec_val,Q),[0 t_max],y0_vec_full);
f_ncr_full = sum(y_ncr_full(:,f_indices),2);

[t_ncr_primary,y_ncr_primary] = ode15s(@(t,y) ncr_solver_v2(t,y,rate_vec_val,Q),[0 t_max],y0_vec_primary);
f_ncr_primary = sum(y_ncr_primary(:,f_indices),2);

[t_ncr_cage,y_ncr_cage] = ode15s(@(t,y) ncr_solver_v2(t,y,rate_vec_val,Q),[0 t_max],y0_vec_cage);
f_ncr_cage = sum(y_ncr_cage(:,f_indices),2);

[t_ncr_neg,y_ncr_neg] = ode15s(@(t,y) ncr_solver_v2(t,y,rate_vec_val,Q),[0 t_max],y0_vec_neg);
f_ncr_neg = sum(y_ncr_neg(:,f_indices),2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Prepare for fits %%%%%%%%%%%%%%%%%%%%%%%%%%

% substitute given values into rate vector
fixed_indices_round1 = find(~ismember(1:numel(reaction_unit_index),fit_indices_round1));
rate_vec_fit_round1 = rate_vec;
for r = fixed_indices_round1
    rr = reaction_unit_index(r);        
    rate_vec_fit_round1 = subs(rate_vec_fit_round1,rr,eval(rr));    
end

% convert rate vector into a matlabFunction
rate_vec_fun_round1 = matlabFunction(rate_vec_fit_round1);

% now set bounds for variables to be inferred
ub_vec = NaN(size(fit_indices_round1));
lb_vec = NaN(size(fit_indices_round1));
for r = 1:numel(reaction_unit_index)
    rate_init = eval(reaction_unit_index(r));        
    ub_vec(r) = 2*sigma_vec(r)*rate_init;
    lb_vec(r) = rate_init/sigma_vec(r)/2;    
end

% specify conditions to include in fit
y0_cell = {y0_vec_primary, y0_vec_neg};%{y0_vec_full, y0_vec_primary, y0_vec_cage, y0_vec_neg};
fluo_cell_raw = {f_ncr_primary, f_ncr_neg};%{f_ncr_full, f_ncr_primary, f_ncr_cage, f_ncr_neg};
time_cell_raw = {t_ncr_primary, t_ncr_neg};%{t_ncr_full, t_ncr_primary, t_ncr_cage, t_ncr_neg};

% define fitting function
fit_fun = @(rate_params,time_exp_array) ode_objfunction_v2(t_max,rate_params,rate_vec_fun_round1,...
    Q,y0_cell,f_indices,time_exp_array);

% generate interpolated vectors for fitting
time_exp = 0:0.1:t_max;
time_exp_array = repmat(time_exp',1,numel(y0_cell));

fluo_exp_array = NaN(numel(time_exp),numel(y0_cell));
for f = 1:numel(fluo_cell_raw)
    fluo_exp_array(:,f) = interp1(time_cell_raw{f},fluo_cell_raw{f},time_exp) + normrnd(0,sim_noise,1,numel(time_exp));
end


% (1) perform optimization on simple non-NCR reaction

% set optimization options
options = optimoptions('lsqcurvefit','MaxIterations',200);

% initialize structures to store results
rate_init_array_r1 = NaN(n_inference_runs,length(fit_indices_round1));
rate_fit_array_r1 = NaN(n_inference_runs,length(fit_indices_round1));
fluo_fit_array_r1 = NaN(size(fluo_exp_array,1),size(fluo_exp_array,2),n_inference_runs);
res_vec_r1 = NaN(1,n_inference_runs);
parfor n = 1:n_inference_runs
    % generate randomized starting values
    iter = 1;    
    rate_init_temp = NaN(1,numel(fit_indices_round1));
    for r = fit_indices_round1
        pass_flag = false;
        init_val = true_param_vec(r);        
        while ~pass_flag
            mFactor = normrnd(0,log10(sigma_vec(r)));
            prop_val = init_val*10^mFactor;
            pass_flag = prop_val <= ub_vec(r) && prop_val >= lb_vec(r);
        end
        rate_init_temp(iter) = prop_val;
        iter = iter + 1;
    end
    % perform fit
    tic
    [rate_vec_out, resnorm]= lsqcurvefit(fit_fun,rate_init_temp,time_exp_array,...
        fluo_exp_array,lb_vec(fit_indices_round1),ub_vec(fit_indices_round1),options);
    toc
    % store results    
    fluo_fit_array_r1(:,:,n) = fit_fun(rate_vec_out,time_exp_array);
    rate_fit_array_r1(n,:) = rate_vec_out;
    rate_init_array_r1(n,:) = rate_init_temp;
    res_vec_r1(n) = resnorm;
end

% Identify iteration with smallest residual and use those parameters to 
%%  fit additional unknowns
[min_res,mi_index] = min(res_vec_r1);
rate_vec_out = rate_fit_array_r1(mi_index,:);
% substitute given values into rate vector
fixed_indices_r2 = find(~ismember(1:numel(reaction_unit_index),fit_indices_round2));
rate_vec_fit_r2 = rate_vec;
iter = 1;
for r = fixed_indices_r2
    rr = reaction_unit_index(r);        
    rate_vec_fit_r2 = subs(rate_vec_fit_r2,rr,rate_vec_out(iter));    
    iter = iter + 1;
end

% convert rate vector into a matlabFunction
rate_vec_fun_r2 = matlabFunction(rate_vec_fit_r2);
% Round 2: Fig NCR curves

% specify conditions to include in fit
y0_cell_r2 = {y0_vec_cage, y0_vec_full};%{y0_vec_full, y0_vec_primary, y0_vec_cage, y0_vec_neg};
fluo_cell_raw_r2 = {f_ncr_cage, f_ncr_full};%{f_ncr_full, f_ncr_primary, f_ncr_cage, f_ncr_neg};
time_cell_raw_r2 = {t_ncr_cage, t_ncr_full};%{t_ncr_full, t_ncr_primary, t_ncr_cage, t_ncr_neg};

% define fitting function
fit_fun_r2 = @(rate_params,time_exp_array) ode_objfunction_v2(t_max,rate_params,rate_vec_fun_r2,...
    Q,y0_cell_r2,f_indices,time_exp_array);

% generate interpolated vectors for fitting
time_exp_array_r2 = repmat(time_exp',1,numel(y0_cell_r2));

fluo_exp_array_r2 = NaN(numel(time_exp),numel(y0_cell_r2));
for f = 1:numel(fluo_cell_raw)
    fluo_exp_array_r2(:,f) = interp1(time_cell_raw_r2{f},fluo_cell_raw_r2{f},time_exp) + normrnd(0,sim_noise,1,numel(time_exp));
end


%% initialize structures to store results
rate_init_array_r2 = NaN(n_inference_runs,length(fit_indices_round2));
rate_fit_array_r2 = NaN(n_inference_runs,length(fit_indices_round2));
fluo_fit_array_r2 = NaN(size(fluo_exp_array_r2,1),size(fluo_exp_array_r2,2),n_inference_runs);
res_vec_r2 = NaN(1,n_inference_runs);
parfor n = 1:n_inference_runs
    % generate randomized starting values
    iter = 1;    
    rate_init_temp = NaN(1,numel(fit_indices_round2));
    for r = fit_indices_round2
        pass_flag = false;
        init_val = true_param_vec(r);        
        while ~pass_flag
            mFactor = normrnd(0,log10(sigma_vec(r)));
            prop_val = init_val*10^mFactor;
            pass_flag = prop_val <= ub_vec(r) && prop_val >= lb_vec(r);
        end
        rate_init_temp(iter) = prop_val;
        iter = iter + 1;
    end
    % perform fit
    tic
    [rate_vec_out, resnorm]= lsqcurvefit(fit_fun_r2,rate_init_temp,time_exp_array,...
        fluo_exp_array_r2,lb_vec(fit_indices_round2),ub_vec(fit_indices_round2),options);
    toc
    % store results    
    fluo_fit_array_r2(:,:,n) = fit_fun_r2(rate_vec_out,time_exp_array);
    rate_fit_array_r2(n,:) = rate_vec_out;
    rate_init_array_r2(n,:) = rate_init_temp;
    res_vec_r2(n) = resnorm;
end

%% Evaluate fit
[~,mi_index2] = min(res_vec_r2);
best_fit_vec = NaN(1,numel(true_param_vec));
best_fit_vec(fit_indices_round1) = rate_fit_array_r1(mi_index,:);
best_fit_vec(fit_indices_round2) = rate_fit_array_r2(mi_index2);

best_fit_vec_ste = NaN(1,numel(true_param_vec));
best_fit_vec_ste(fit_indices_round1) = std(rate_fit_array_r1);
best_fit_vec_ste(fit_indices_round2) = std(rate_fit_array_r2);

best_fit_fluo_r1 = fit_fun(rate_fit_array_r1(mi_index,:),time_exp_array);
best_fit_fluo_r2 = fit_fun_r2(rate_fit_array_r2(mi_index2),time_exp_array);

%=% make fit figure
fit_fig1 = figure;
cmap1 = brewermap(9,'Set2');
hold on
p = [];
a = [];
for f = 1:size(fluo_exp_array,2)
    a = [a scatter(time_exp(1:10:end),fluo_exp_array(1:10:end,f),5,'MarkerFaceColor',cmap1(f,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.7)];
%     a = [a plot(time_exp,fluo_exp_array(:,f),'Color',cmap1(f,:),'LineWidth',1.5)];
    p = [p plot(time_exp,best_fit_fluo_r1(:,f),'Color','black','LineWidth',1.5)];
end
legend([a p(1)],'Primary Only','No Activator','Model Fits','Location','southeast')
ylim([0 210])
xlim([0 t_max])
xlabel('time (seconds)')
ylabel('signal (au)')
set(gca,'Fontsize',14)

saveas(fit_fig1,[FigPath 'NCR_fit_plots_r1.png'])

fit_fig2 = figure;
cmap1 = brewermap(9,'Set2');
hold on
p = [];
a = [];
for f = 1:size(fluo_exp_array_r2,2)
    a = [a scatter(time_exp(1:10:end),fluo_exp_array_r2(1:10:end,f),5,'MarkerFaceColor',cmap1(f+2,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.7)];
%     a = [a plot(time_exp,fluo_exp_array_r2(:,f),'Color',cmap1(f+2,:),'LineWidth',1.5)];
    p = [p plot(time_exp,best_fit_fluo_r2(:,f),'Color','black','LineWidth',1.5)];
end
legend([a p(1)],'Cage Only','NCR','Model Fits','Location','southeast')
ylim([0 210])
xlim([0 t_max])
xlabel('time (seconds)')
ylabel('signal (au)')
set(gca,'Fontsize',14)

saveas(fit_fig2,[FigPath 'NCR_fit_plots_r2.png'])

fit_fig_full = figure;
cmap1 = brewermap(9,'Set2');
hold on
p = [];
for f = 1:size(fluo_exp_array,2)
    scatter(time_exp(1:10:end),fluo_exp_array(1:10:end,f),5,'MarkerFaceColor',cmap1(f,:)/1.2,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.1)
    p = [p plot(time_exp,best_fit_fluo_r1(:,f),'Color',cmap1(f,:),'LineWidth',1.5)];
end
for f = 1:size(fluo_exp_array_r2,2)
    scatter(time_exp(1:10:end),fluo_exp_array_r2(1:10:end,f),5,'MarkerFaceColor',cmap1(f+2,:)/1.2,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.1)
    p = [p plot(time_exp,best_fit_fluo_r2(:,f),'Color',cmap1(f+2,:),'LineWidth',1.5)];
end
legend(p,'Primary Only','No Activator','Cage Only','NCR','Location','southeast')
ylim([0 210])
xlim([0 t_max])
xlabel('time (seconds)')
ylabel('signal (au)')
set(gca,'Fontsize',14)

saveas(fit_fig_full,[FigPath 'NCR_fit_plots.png'])
%%
param_fig = figure;
hold on
% e = errorbar(true_param_vec,best_fit_vec, best_fit_vec_ste,'o', 'Color','black');
% e.CapSize = 0;
scatter(true_param_vec,best_fit_vec,'MarkerFaceColor',cmap1(3,:),'MarkerEdgeColor','black')
% scatter(best_fit_vec_med,true_param_vec,'MarkerFaceColor',cmap1(2,:),'MarkerEdgeColor','black')
set(gca,'YScale','log')
set(gca,'XScale','log')
grid on
ylabel('inferred parameter value')
xlabel('true parameter value')
set(gca,'Fontsize',12)
ylim([1e-7 1e5])
xlim([1e-7 1e5])
saveas(param_fig,[FigPath 'NCR_param_plots.png'])