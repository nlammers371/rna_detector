clear 
close all

% basic path info
addpath(genpath('../utilities'))
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
reaction_parameter_index = unique(reaction_units);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% set key parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
t_max = 2000; % duration of time to solve for (need to change this to calculate dynamically)
t_res = 0.1;
A10 = 0.5; % concentration of target RNA
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
fit_indices = 1:numel(reaction_parameter_index);
n_inference_runs = 20;
sim_noise = 1e-2*S0;

% define "true" rates to fit
b = 1e-6; % Cas13 1/SNR
k = 0.025; % association rate (s^-1 nM^-1)
kc = 100; % Cas13 catalytic rate
kd_cga = 1e-8;
rcg = 5e-3; % off rate for Cas13:gRNA
rcga = 1e-5; % off rate for gRNA:Activator
rga = 1e-2;
ria = 1e-2; % off rate for CleavedCage:Activator
rns = 1e3; % off rate for all nonspecific interactions

% get set of true param values
true_param_vec = eval(reaction_parameter_index);

% specify degree of uncertainty for each param value
param_sigma_vec = ones(1,9)*1e3;

% rate vector to substitute into
rate_vec_val = rate_vec;
% generate valued long rate vec
for r = 1:numel(reaction_parameter_index)
    rr = reaction_parameter_index(r);        
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
fixed_indices = find(~ismember(1:numel(reaction_parameter_index),fit_indices));
rate_vec_fit_round1 = rate_vec;
for r = fixed_indices
    rr = reaction_parameter_index(r);        
    rate_vec_fit_round1 = subs(rate_vec_fit_round1,rr,eval(rr));    
end

% convert rate vector into a matlabFunction
rate_vec_fun = matlabFunction(rate_vec_fit_round1);

% now set bounds for variables to be inferred
ub_vec = NaN(size(fit_indices));
lb_vec = NaN(size(fit_indices));
for r = 1:numel(reaction_parameter_index)
    rate_init = eval(reaction_parameter_index(r));        
    ub_vec(r) = 2*param_sigma_vec(r)*rate_init;
    lb_vec(r) = rate_init/param_sigma_vec(r)/2;    
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

% define helper functions for sampling
mu_fun = @(mu_in,s_factor) log(mu_in.^2 ./ sqrt(s_factor.^2.*mu_in.^2 + mu_in.^2));
sigma_fun = @(mu_in,s_factor) sqrt(log(s_factor.^2.*mu_in.^2./mu_in.^2 + 1));

% define MCMC parameters
temp_init = 1; % impacts rate of acceptance for jumps (will be tuned)
sample_mem = 10; % sets # of previous samples to use for calc acceptance rate
acceptance_target = 0.3; % target rate at which proposals are accepted
sample_batch_size = 3; % number of parameters to sample for proposal step
jump_sigma = .05; % defined as fraction of param value
n_steps = 1e3; % number of sampling steps
n_chains = 1; % number of chains to run in parallel (kept at 1 for now)
noise = 6; % experimental error--sets scale for fit likelihood calculations

% initialize arrays to store results
param_array = NaN(n_steps,numel(reaction_parameter_index),n_chains); % parameters
logL_array = NaN(n_steps,n_chains); % likelihood score
temperature_array = NaN(n_steps,n_chains); % keep track of temp over time
deltaL_array = NaN(n_steps,n_chains);
acc_rate_array = NaN(n_steps,n_chains); % keep track of acceptance rate

f = waitbar(0,'generating mcmc samples...');

for chain = 1:n_chains
    % select starting values
    param_array(1,fixed_indices,chain) = true_param_vec(fixed_indices);
    pass_flag = false;
    while ~pass_flag
        mu_vec = mu_fun(true_param_vec,param_sigma_vec);
        sigma_vec = sigma_fun(true_param_vec,param_sigma_vec);
        prop_param_vals = lognrnd(mu_vec,sigma_vec);
        pass_flag = all(prop_param_vals<=ub_vec&prop_param_vals>=lb_vec);   
    end
    param_array(1,:) = prop_param_vals;
    
    % initialize
    params_current = param_array(1,:);
    prediction_current = fit_fun(params_current,time_exp_array);
    logL_current = -sum(0.5*((fluo_exp_array(:) - prediction_current(:))/noise).^2);
    acc_rate_array(1,chain) = true;
    temperature_array(1) = temp_init;

    % iterate through sampling steps
    for step = 2:n_steps
        waitbar(step/n_steps,f)
%         tic                   
        %%%%%%%%%%%%%%%%%%%%%%%%
        % propose new move
        %%%%%%%%%%%%%%%%%%%%%%%%
        % select params to sample
        prop_param_indices = randsample(fit_indices,sample_batch_size,false);
        old_param_vals = params_current(prop_param_indices);

        pass_flag = false;
        iter = 1;
        while ~pass_flag
            % generate new prop vals
            m_factors = normrnd(0,jump_sigma,1,sample_batch_size);
            prop_param_vals = old_param_vals.*10.^m_factors;
            % check that vals are inside bounds
            pass_flag = all(prop_param_vals<=ub_vec(prop_param_indices)&prop_param_vals>=lb_vec(prop_param_indices));
%             m_factors = normrnd(0,jump_sigma,1,sample_batch_size);
%             mu_vec = mu_fun(old_param_vals,jump_sigma);
%             sigma_vec = sigma_fun(old_param_vals,jump_sigma);
%             prop_param_vals = lognrnd(mu_vec,sigma_vec);%old_param_vals.*10.^m_factors;
%             % check that vals are inside bounds
%             pass_flag = all(prop_param_vals<=ub_vec(prop_param_indices)&prop_param_vals>=lb_vec(prop_param_indices));
%             iter = iter + 1;
%             if iter > 100
%                 error('asfa')
%             end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%
        % perform MH step 
        %%%%%%%%%%%%%%%%%%%%%%%%
        params_prop = params_current;
        params_prop(prop_param_indices) = prop_param_vals;
        prediction_prop = fit_fun(params_prop,time_exp_array);
        logL_prop = -sum(0.5*((fluo_exp_array(:) - prediction_prop(:))/noise).^2);

        % calculate likelhood diff
        deltaL = (logL_prop-logL_current) / .1;    
        deltaL_array(step,chain) = deltaL;
        move_flag = exp(deltaL/temp_init) > rand();

        if move_flag        
            % load prediction into memory
            prediction_current = prediction_prop;
            logL_current = logL_prop;
            params_current = params_prop;
        end

        param_array(step,:) = params_current;
        logL_array(step,chain) = logL_current;
        temperature_array(step,chain) = temp_init;
        acc_rate_array(step,chain) = move_flag;
%         toc
    end
end   














