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
t_max = 3600; % duration of time to solve for (need to change this to calculate dynamically)
t_res = 0.1;
A10 = 1e-1; % concentration of target RNA
IA20 = 20; % caged secondary activator
S0 = 200; % reporter substrate
RNP1 = 20; % Cas13:gRNA1
RNP2 = 0.8; % Cas13:gRNA2

%%%% specify initial reactant concentrations 
% vectors for substitutions
init_sym_array =      {'A1'  'S' 'IA2' 'C13'       'G1'   'G2'};
init_val_vec_cage =    [0    S0  IA20   RNP1+RNP2   RNP1  RNP2]; % cage only
init_val_vec_full =    [A10  S0  IA20   RNP1+RNP2   RNP1  RNP2]; % primary only
init_val_vec_primary = [A10  S0  0      RNP1+RNP2   RNP1  RNP2]; % primary only
init_val_vec_neg =     [0    S0  0      RNP1+RNP2   RNP1  RNP2]; % no activator

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
y0_cell = {y0_vec_full y0_vec_primary y0_vec_cage y0_vec_neg};
name_cell = {'NCR', 'Primary Only', 'Cage Only', 'No Activator'};


%%%%%%%%%%%%%%%%%
%%% Specify which rate parameters to fit and specify values for "given"
%%% parameters

% specify subset of parameters to fit
fit_indices = 1:7;
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
time_cell_raw = cell(1,numel(y0_cell));
fluo_cell_raw = cell(1,numel(y0_cell));
raw_full_output_cell = cell(1,numel(y0_cell));
for y0 = 1:numel(y0_cell)
    [t_ncr,y_ncr] = ode15s(@(t,y) ncr_solver_v2(t,y,rate_vec_val,Q),[0 t_max],y0_cell{y0});
    time_cell_raw{y0} = t_ncr;
    fluo_cell_raw{y0} = sum(y_ncr(:,f_indices),2);
    raw_full_output_cell{y0} = y_ncr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Prepare for fits %%%%%%%%%%%%%%%%%%%%%%%%%%

% substitute given values into rate vector
fixed_indices = find(~ismember(1:numel(reaction_unit_index),fit_indices));
rate_vec_fit_round1 = rate_vec;
for r = fixed_indices
    rr = reaction_unit_index(r);        
    rate_vec_fit_round1 = subs(rate_vec_fit_round1,rr,eval(rr));    
end

% convert rate vector into a matlabFunction
rate_vec_fun = matlabFunction(rate_vec_fit_round1);

% now set bounds for variables to be inferred
ub_vec = NaN(size(fit_indices));
lb_vec = NaN(size(fit_indices));
for r = 1:numel(reaction_unit_index)
    rate_init = eval(reaction_unit_index(r));        
    ub_vec(r) = 2*sigma_vec(r)*rate_init;
    lb_vec(r) = rate_init/sigma_vec(r)/2;    
end

% define fitting function
fit_fun = @(rate_params,time_exp_array) ode_objfunction_v2(t_max,rate_params,rate_vec_fun,...
    Q,y0_cell,f_indices,time_exp_array);

% generate interpolated vectors for fitting
time_exp = 0:0.1:t_max;
time_exp_array = repmat(time_exp',1,numel(y0_cell));

fluo_exp_array = NaN(numel(time_exp),numel(y0_cell));
for f = 1:numel(fluo_cell_raw)
    fluo_exp_array(:,f) = interp1(time_cell_raw{f},fluo_cell_raw{f},time_exp) + normrnd(0,sim_noise,1,numel(time_exp));
end


%% perform MCMC optimization 

% define MCMC parameters
temp_init = 300; % impacts rate of acceptance for jumps (will be tuned)
sample_mem = 10; % sets # of previous samples to use for calc acceptance rate
acceptance_target = 0.3; % target rate at which proposals are accepted
sample_batch_size = 3; % number of parameters to sample for proposal step
jump_sigma = .1; % defined as fraction of param value
n_steps = 1e2; % number of sampling steps
n_chains = 1; % number of chains to run in parallel (kept at 1 for now)
noise = 5; % experimental error--sets scale for fit likelihood calculations

% initialize arrays to store results
param_array = NaN(n_steps,numel(reaction_unit_index),n_chains); % parameters
logL_array = NaN(n_steps,n_chains); % likelihood score
temperature_array = NaN(n_steps,n_chains); % keep track of temp over time
acc_rate_array = NaN(n_steps,n_chains); % keep track of acceptance rate

chain = 1;

parfor chain = 1:n_chains
    % select starting values
    param_array(1,fixed_indices,chain) = true_param_vec(fixed_indices);
    for r = fit_indices
        pass_flag = false;
        init_val = true_param_vec(r);        
        while ~pass_flag
            mFactor = normrnd(0,log10(sigma_vec(r)));
            prop_val = init_val*10^mFactor;
            pass_flag = prop_val <= ub_vec(r) && prop_val >= lb_vec(r);
        end
        param_array(1,r) = prop_val;
    end

    % initialize
    params_current = param_array(1,:);
    prediction_current = fit_fun(params_current,time_exp_array);
    logL_current = -mean(0.5*((fluo_exp_array(:) - prediction_current(:))/noise).^2);
    acc_rate_array(1,chain) = true;
    temperature_array(1) = true;

    % iterate through sampling steps
    for step = 2:n_steps

        tic    
        %%%%%%%%%%%%%%%%%%%%%%%%
        % propose new move
        %%%%%%%%%%%%%%%%%%%%%%%%
        % select params to sample
        prop_param_indices = randsample(fit_indices,sample_batch_size,false);
        old_param_vals = params_current(prop_param_indices);

        pass_flag = false;
        while ~pass_flag
            % generate new prop vals
            m_factors = normrnd(0,jump_sigma,1,sample_batch_size);
            prop_param_vals = old_param_vals.*10.^m_factors;
            % check that vals are inside bounds
            pass_flag = all(prop_param_vals<=ub_vec(prop_param_indices)&prop_param_vals>=lb_vec(prop_param_indices));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%
        % perform MH step 
        %%%%%%%%%%%%%%%%%%%%%%%%
        params_prop = params_current;
        params_prop(prop_param_indices) = prop_param_vals;
        prediction_prop = fit_fun(params_prop,time_exp_array);
        logL_prop = -mean(0.5*((fluo_exp_array(:) - prediction_prop(:))/noise).^2);

        % calculate likelhood diff
        deltaL = (logL_prop-logL_current) / .1;    
        move_flag = exp(deltaL) > rand();

        if move_flag        
            % load prediction into memory
            prediction_current = prediction_prop;
            logL_current= logL_prop;
            params_current = params_prop;
        end

        param_array(step,:) = params_current;
        logL_array(step,chain) = logL_current;
        temperature_array(step,chain) = temp_init;
        acc_rate_array(step,chain) = move_flag;
        toc
    end
end   














% % set optimization options
% options = optimoptions('lsqcurvefit','MaxIterations',200);
% 
% % initialize structures to store results
% rate_init_array_r1 = NaN(n_inference_runs,length(fit_indices_round1));
% rate_fit_array_r1 = NaN(n_inference_runs,length(fit_indices_round1));
% fluo_fit_array_r1 = NaN(size(fluo_exp_array,1),size(fluo_exp_array,2),n_inference_runs);
% res_vec_r1 = NaN(1,n_inference_runs);
% parfor n = 1:n_inference_runs
%     % generate randomized starting values
%     iter = 1;    
%     rate_init_temp = NaN(1,numel(fit_indices_round1));
%     for r = fit_indices_round1
%         pass_flag = false;
%         init_val = true_param_vec(r);        
%         while ~pass_flag
%             mFactor = normrnd(0,log10(sigma_vec(r)));
%             prop_val = init_val*10^mFactor;
%             pass_flag = prop_val <= ub_vec(r) && prop_val >= lb_vec(r);
%         end
%         rate_init_temp(iter) = prop_val;
%         iter = iter + 1;
%     end
%     % perform fit
%     tic
%     [rate_vec_out, resnorm]= lsqcurvefit(fit_fun,rate_init_temp,time_exp_array,...
%         fluo_exp_array,lb_vec(fit_indices_round1),ub_vec(fit_indices_round1),options);
%     toc
%     % store results    
%     fluo_fit_array_r1(:,:,n) = fit_fun(rate_vec_out,time_exp_array);
%     rate_fit_array_r1(n,:) = rate_vec_out;
%     rate_init_array_r1(n,:) = rate_init_temp;
%     res_vec_r1(n) = resnorm;
% end
% 
% % Identify iteration with smallest residual and use those parameters to 
% %%  fit additional unknowns
% [min_res,mi_index] = min(res_vec_r1);
% rate_vec_out = rate_fit_array_r1(mi_index,:);
% % substitute given values into rate vector
% fixed_indices_r2 = find(~ismember(1:numel(reaction_unit_index),fit_indices_round2));
% rate_vec_fit_r2 = rate_vec;
% iter = 1;
% for r = fixed_indices_r2
%     rr = reaction_unit_index(r);        
%     rate_vec_fit_r2 = subs(rate_vec_fit_r2,rr,rate_vec_out(iter));    
%     iter = iter + 1;
% end
% 
% % convert rate vector into a matlabFunction
% rate_vec_fun_r2 = matlabFunction(rate_vec_fit_r2);
% % Round 2: Fig NCR curves
% 
% % specify conditions to include in fit
% y0_cell_r2 = {y0_vec_cage, y0_vec_full};%{y0_vec_full, y0_vec_primary, y0_vec_cage, y0_vec_neg};
% fluo_cell_raw_r2 = {f_ncr_cage, f_ncr_full};%{f_ncr_full, f_ncr_primary, f_ncr_cage, f_ncr_neg};
% time_cell_raw_r2 = {t_ncr_cage, t_ncr};%{t_ncr_full, t_ncr_primary, t_ncr_cage, t_ncr_neg};
% 
% % define fitting function
% fit_fun_r2 = @(rate_params,time_exp_array) ode_objfunction_v2(t_max,rate_params,rate_vec_fun_r2,...
%     Q,y0_cell_r2,f_indices,time_exp_array);
% 
% % generate interpolated vectors for fitting
% time_exp_array_r2 = repmat(time_exp',1,numel(y0_cell_r2));
% 
% fluo_exp_array_r2 = NaN(numel(time_exp),numel(y0_cell_r2));
% for f = 1:numel(fluo_cell_raw)
%     fluo_exp_array_r2(:,f) = interp1(time_cell_raw_r2{f},fluo_cell_raw_r2{f},time_exp) + normrnd(0,sim_noise,1,numel(time_exp));
% end
% 
% 
% %% initialize structures to store results
% rate_init_array_r2 = NaN(n_inference_runs,length(fit_indices_round2));
% rate_fit_array_r2 = NaN(n_inference_runs,length(fit_indices_round2));
% fluo_fit_array_r2 = NaN(size(fluo_exp_array_r2,1),size(fluo_exp_array_r2,2),n_inference_runs);
% res_vec_r2 = NaN(1,n_inference_runs);
% parfor n = 1:n_inference_runs
%     % generate randomized starting values
%     iter = 1;    
%     rate_init_temp = NaN(1,numel(fit_indices_round2));
%     for r = fit_indices_round2
%         pass_flag = false;
%         init_val = true_param_vec(r);        
%         while ~pass_flag
%             mFactor = normrnd(0,log10(sigma_vec(r)));
%             prop_val = init_val*10^mFactor;
%             pass_flag = prop_val <= ub_vec(r) && prop_val >= lb_vec(r);
%         end
%         rate_init_temp(iter) = prop_val;
%         iter = iter + 1;
%     end
%     % perform fit
%     tic
%     [rate_vec_out, resnorm]= lsqcurvefit(fit_fun_r2,rate_init_temp,time_exp_array,...
%         fluo_exp_array_r2,lb_vec(fit_indices_round2),ub_vec(fit_indices_round2),options);
%     toc
%     % store results    
%     fluo_fit_array_r2(:,:,n) = fit_fun_r2(rate_vec_out,time_exp_array);
%     rate_fit_array_r2(n,:) = rate_vec_out;
%     rate_init_array_r2(n,:) = rate_init_temp;
%     res_vec_r2(n) = resnorm;
% end
% 
% %% Evaluate fit
% [~,mi_index2] = min(res_vec_r2);
% best_fit_vec = NaN(1,numel(true_param_vec));
% best_fit_vec(fit_indices_round1) = rate_fit_array_r1(mi_index,:);
% best_fit_vec(fit_indices_round2) = rate_fit_array_r2(mi_index2);
% 
% best_fit_fluo_r1 = fit_fun(rate_fit_array_r1(mi_index,:),time_exp_array);
% best_fit_fluo_r2 = fit_fun_r2(rate_fit_array_r2(mi_index2),time_exp_array);
% 
% % make fit figure
% fit_fig = figure;
% cmap1 = brewermap(9,'Set2');
% hold on
% p = [];
% for f = 1:size(fluo_exp_array,2)
%     scatter(time_exp(1:10:end),fluo_exp_array(1:10:end,f),5,'MarkerFaceColor',cmap1(f,:)/1.2,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.1)
%     p = [p plot(time_exp,best_fit_fluo_r1(:,f),'Color',cmap1(f,:),'LineWidth',1.5)];
% end
% for f = 1:size(fluo_exp_array_r2,2)
%     scatter(time_exp(1:10:end),fluo_exp_array_r2(1:10:end,f),5,'MarkerFaceColor',cmap1(f+2,:)/1.2,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.1)
%     p = [p plot(time_exp,best_fit_fluo_r2(:,f),'Color',cmap1(f+2,:),'LineWidth',1.5)];
% end
% legend(p,'Primary Only','No Activator','Cage Only','NCR','Location','southeast')
% ylim([0 210])
% xlim([0 t_max])
% xlabel('time (seconds)')
% ylabel('signal (au)')
% set(gca,'Fontsize',14)
% 
% saveas(fit_fig,[FigPath 'NCR_fit_plots.png'])
% %%
% param_fig = figure;
% scatter(best_fit_vec,true_param_vec,'MarkerFaceColor',cmap1(3,:),'MarkerEdgeColor','black')
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% grid on
% xlabel('inferred parameter value')
% ylabel('true parameter value')
% set(gca,'Fontsize',12)
% ylim([1e-7 1e5])
% xlim([1e-7 1e5])
% saveas(param_fig,[FigPath 'NCR_param_plots.png'])