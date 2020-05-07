% Script to use numerical ODE solver and nonlinear least square fitting to
% estimate reaction parameters from experimental data
clear 
close all

% basic path info
addpath('../utilities')
DataPath = '../../out/ode_studies_v2/';
mkdir(DataPath)

% specify project to load
project = 'ncr_basic';

% make figure path
FigPath = ['../../fig/ode_studies_v2/' project '/' ];
mkdir(FigPath)

% load 
load([DataPath project '.mat'])

% extract reaction parameters that need to be specified or fit
reaction_rate_index = unique(rate_vec);

% decompose rates into fundamental units
reaction_units = [];
for i = 1:numel(reaction_rate_index)
    r_sym = reaction_rate_index(i);
    t_char = char(r_sym);
    char_list = strsplit(t_char,'*');
    for c = 1:numel(char_list)
%         char_list{c}
        eval(['syms ' char_list{c}])
        reaction_units = [reaction_units eval(char_list{c})];
    end
end
reaction_parameter_index = unique(reaction_units); % list of unique parameters
reaction_parameter_index


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%specify conditions to fit (ideally this info would be loaded with data)

exp_data = readtable('./data/042520 20nM wobble+mismatch raw data_full_by_cond_no extra_rows.csv');

num_replicates = 3;

t_max = max(exp_data.time_s); % duration of time to solve for
A10 = 0.02; % initial target RNA concentration (nM)
S0 = 200; % reporter concentration
%S0 = 3500; % reporter concentration
RNP1 = 20; 
IA2 = 20;
RNP2 = 0.05;
fluorescence_multiplier = 8000;


%%%% specify initial reactant concentrations 
% vectors for substitutions
init_sym_array =     {'A1' 'S' 'C13' 'G1' 'G2' 'IA2'};
init_val_vec_primary =   [A10  S0  (RNP1+RNP2) (2*RNP1) (2*RNP2) 0]; % primary only
init_val_vec_neg =    [0  S0  (RNP1+RNP2) (2*RNP1) (2*RNP2) 0]; % no activators
init_val_vec_cage =   [0  S0  (RNP1+RNP2) (2*RNP1) (2*RNP2) IA2]; % Cage only
%init_val_vec_ncr =    [A10  S0  (RNP1+RNP2) (2*RNP1) (2*RNP2) IA2]; % NCR

% initialize concentration vector. Everything not explicitly assigned taken
% to be zero
y0_vec_primary = zeros(size(full_reactant_list));
y0_vec_neg = zeros(size(full_reactant_list));
y0_vec_cage = zeros(size(full_reactant_list));
%y0_vec_ncr = zeros(size(full_reactant_list));
for i = 1:numel(init_sym_array)
    species = init_sym_array{i};
    spec_indices = find(strcmp(full_reactant_list,species));      
    y0_vec_primary(spec_indices) = init_val_vec_primary(i);    
    y0_vec_neg(spec_indices) = init_val_vec_neg(i);
    y0_vec_cage(spec_indices) = init_val_vec_cage(i);
    %y0_vec_ncr(spec_indices) = init_val_vec_ncr(i);
end

y0_cell = {y0_vec_primary, y0_vec_neg, y0_vec_cage};%, y0_vec_ncr}; % all conditions in this cell will be included in fit routine

%%%%%%%%%%%%%%%%%
%%% Specify which rate parameters to fit and specify values for "given"
%%% parameters
% reaction_parameter_index
% full_reactant_list
% return
% specify subset of parameters to fit
fit_indices = [1,2,3,4,5,6,7,8,9];
n_inference_runs = 20;
sim_noise = 1e-2*S0; % gaussian noise to add to simulated data

% define initial guess for parameters to fit
b = 6.3698e-06; % Cas13 1/SNR
k =  0.023604; % association rate (s^-1 nM^-1)
kc = 638.66; % Cas13 catalytic rate 
kd_cga =  1.4714e-08; %dissassociation constant for cas13:guide and activatior 
rcg =  0.0050912; % off rate for Cas13:gRNA      
rcga = 6.9116e-05; %on rate for cas13:guide and activatior      
rga =  0.062162; % off rate for gRNA:Activator
rns = 1313.1; % off rate for all nonspecific interactions

ria = 1e3; %off rate for cage and activator

% get vector of guessed param values
true_param_vec = eval(reaction_parameter_index);
% kd_cga_indices = find(strcmp(reaction_parameter_index,{'ka_cga'}));
% true_param_vec(kd_cga_indices) = 0;

% specify degree of uncertainty for each param value (need to play with
% this)
% param order [ b, k, kc, kd_cga, rcg, rcga, rga, ria, rns]
sigma_vec = [1e3, 1e1, 1e1, 1e3, 1e1, 1e3, 1e1, 1e3, 1e1];

%%%%%%%%%%% Generate "experimental data" using true param values %%%%%%%%%%

%{
% rate vector to substitute into
rate_vec_val = rate_vec;

% generate valued long rate vec
for r = 1:numel(reaction_parameter_index)
    rr = reaction_parameter_index(r);        
    rate_vec_val = subs(rate_vec_val,rr,eval(rr));    
end
rate_vec_val = double(rate_vec_val);

close all
f_indices = find(strcmp(full_reactant_list,{'F'})); % get indices of fluorescent reporter

% specify conditions to include in fit

fluo_cell_raw = cell(1,length(y0_cell));
time_cell_raw = cell(1,length(y0_cell));

for f = 1:length(y0_cell)
    [t_fit,y_fit] = ode15s(@(t,y) ncr_solver_v2(t,y,rate_vec_val,Q),[0 t_max],y0_cell{f});
    fluo_cell_raw{f} = sum(y_fit(:,f_indices),2);
    time_cell_raw{f} = t_fit;
end
%}
%%%%%%%%%%%%%%%%%%%%%% Process exp data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_indices = find(strcmp(full_reactant_list,{'F'})); % get indices of fluorescent reporter

% Make locations for data
fluo_cell_raw = cell(1,length(y0_cell)*num_replicates);
time_cell_raw = cell(1,length(y0_cell)*num_replicates);

%adding primary only data
fluo_cell_raw{1}=exp_data.x1_20nM_0_05nM_0_02_20____-exp_data.x1_20nM_0_05nM_0_02_20____(1);
fluo_cell_raw{2}=exp_data.x1_20nM_0_05nM_0_02_20_____1-exp_data.x1_20nM_0_05nM_0_02_20_____1(1);
fluo_cell_raw{3}=exp_data.x1_20nM_0_05nM_0_02_20_____2-exp_data.x1_20nM_0_05nM_0_02_20_____2(1);

%adding neg data
fluo_cell_raw{4}=exp_data.x1_20nM_0_05nM_0_0_20____-exp_data.x1_20nM_0_05nM_0_0_20____(1);
fluo_cell_raw{5}=exp_data.x1_20nM_0_05nM_0_0_20_____1-exp_data.x1_20nM_0_05nM_0_0_20_____1(1);
fluo_cell_raw{6}=exp_data.x1_20nM_0_05nM_0_0_20_____2-exp_data.x1_20nM_0_05nM_0_0_20_____2(1);

%adding cage data
fluo_cell_raw{7}=exp_data.x1_20nM_0_05nM_0_0_20_NCR061-exp_data.x1_20nM_0_05nM_0_0_20_NCR061(1);
fluo_cell_raw{8}=exp_data.x1_20nM_0_05nM_0_0_20_NCR061_1-exp_data.x1_20nM_0_05nM_0_0_20_NCR061_1(1);
fluo_cell_raw{9}=exp_data.x1_20nM_0_05nM_0_0_20_NCR061_2-exp_data.x1_20nM_0_05nM_0_0_20_NCR061_2(1);

%adding ncr data
% fluo_cell_raw{10}=exp_data.x1_20nM_0_05nM_0_02_20_NCR061-exp_data.x1_20nM_0_05nM_0_02_20_NCR061(1);
% fluo_cell_raw{11}=exp_data.x1_20nM_0_05nM_0_02_20_NCR061_1-exp_data.x1_20nM_0_05nM_0_02_20_NCR061_1(1);
% fluo_cell_raw{12}=exp_data.x1_20nM_0_05nM_0_02_20_NCR061_2-exp_data.x1_20nM_0_05nM_0_02_20_NCR061_2(1);

%adding times
time_cell_raw{1}=exp_data.time_s;
time_cell_raw{2}=exp_data.time_s;
time_cell_raw{3}=exp_data.time_s;

time_cell_raw{4}=exp_data.time_s;
time_cell_raw{5}=exp_data.time_s;
time_cell_raw{6}=exp_data.time_s;

time_cell_raw{7}=exp_data.time_s;
time_cell_raw{8}=exp_data.time_s;
time_cell_raw{9}=exp_data.time_s;

% time_cell_raw{10}=exp_data.time_s;
% time_cell_raw{11}=exp_data.time_s;
% time_cell_raw{12}=exp_data.time_s;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Prepare for fits %%%%%%%%%%%%%%%%%%%%%%%%%%

% substitute given values into rate vector
fixed_indices = find(~ismember(1:numel(reaction_parameter_index),fit_indices));
rate_vec_fit = rate_vec;
for r = fixed_indices
    rr = reaction_parameter_index(r);        
    rate_vec_fit = subs(rate_vec_fit,rr,eval(rr));    
end

% convert rate vector into a matlabFunction
rate_vec_fun = matlabFunction(rate_vec_fit);

% now set bounds for variables to be inferred
ub_vec = NaN(size(fit_indices));
lb_vec = NaN(size(fit_indices));
for r = 1:numel(reaction_parameter_index)
    rate_init = eval(reaction_parameter_index(r));        
    ub_vec(r) = 2*sigma_vec(r)*rate_init;
    lb_vec(r) = rate_init/sigma_vec(r)/2;    
end


% define fitting function
fit_fun = @(rate_params,time_exp_array) (fluorescence_multiplier/200)*ode_objfunction_v2(t_max,rate_params,rate_vec_fun,...
    Q,repelem(y0_cell, num_replicates),f_indices,time_exp_array);

% generate interpolated vectors for fitting
time_exp = 0:0.1:t_max; % need to make time res specification dynamic
time_exp_array = repmat(time_exp',1,numel(y0_cell)*num_replicates);

fluo_exp_array = NaN(numel(time_exp),numel(y0_cell)*num_replicates);
for f = 1:numel(fluo_cell_raw)
    fluo_exp_array(:,f) = interp1(time_cell_raw{f},fluo_cell_raw{f},time_exp);
end

%%%%%%%%%%%%%%%%
% perform optimization on simple non-NCR reaction

% set optimization options
options = optimoptions('lsqcurvefit','MaxIterations',200);

% initialize structures to store results
param_init_array = NaN(n_inference_runs,length(fit_indices)); % store initial values used
param_fit_array = NaN(n_inference_runs,length(fit_indices)); % store parameters
fluo_fit_array = NaN(size(fluo_exp_array,1),size(fluo_exp_array,2),n_inference_runs); % store fluo curves
res_vec = NaN(1,n_inference_runs); % store fit residuals
parfor n = 1:n_inference_runs
    % generate randomized starting values
    iter = 1;    
    rate_init_temp = NaN(1,numel(fit_indices));
    for r = fit_indices
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
        fluo_exp_array,lb_vec(fit_indices),ub_vec(fit_indices),options);
    toc
    % store results    
    fluo_fit_array(:,:,n) = fit_fun(rate_vec_out,time_exp_array);
    param_fit_array(n,:) = rate_vec_out;
    param_init_array(n,:) = rate_init_temp;
    res_vec(n) = resnorm;
end


[~,best_fit_index] = min(res_vec);

params = param_fit_array(best_fit_index, :)
param_order = reaction_parameter_index

figure;
hold on
p = plot(time_exp,fluo_exp_array);
q = plot(time_exp,fluo_fit_array(:,:,best_fit_index),'Color','black','LineWidth',1.5);
legend([p(1) p(2) q(1)],'primary only','no activator', 'model fits','Location','northwest')

figure(2);
p = plot(time_exp,fluo_exp_array-fluo_fit_array(:,:,best_fit_index));
legend([p(1) p(2)],'primary only','no activator','Location','northwest');