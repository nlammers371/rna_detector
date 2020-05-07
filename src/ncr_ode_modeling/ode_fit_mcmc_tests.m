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
t_max = 3600; % duration of time to solve for (need to change this to calculate dynamically)
t_res = 0.1;
A10 = 0.2; % concentration of target RNA
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
fit_indices = [1:7 8 9];%numel(reaction_parameter_index);
n_inference_runs = 20;
sim_noise = 1e-2*S0;

% define "true" rates to fit
b = 1e-6; % Cas13 1/SNR
k = 0.025; % association rate (s^-1 nM^-1)
kc = 600; % Cas13 catalytic rate
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

% generate data structure for mcmc algorithm to use
data = struct;
data.ydata = fluo_exp_array;
data.xdata = time_exp_array;

% high-level simulation paramters 
method      = 'dram'; % adaptation method, 'mh', 'dr', 'am', or 'dram'
nsimu       = 1000;   % number of simulations
adaptint    = 100;    % how often to adapt the proposal

% set options 
options.method      = method;        % adaptation method (mh,am,dr,dram)
options.nsimu       = nsimu;         % n:o of simulations
options.qcov        = eye(11)*0.001; % proposal covariance
options.adaptint    = adaptint; % adaptation interval
options.printint    = 200; % how often to show info on acceptance ratios
options.verbosity   = 0;  % how much to show output in Matlab window
options.waitbar     = 1;  % show garphical waitbar
options.updatesigma = 1;  % update error variance
options.stats       = 1;  % save extra statistics in results



% specify fit funtions
modelfun = @(xdata,theta) ode_objfunction_v2(t_max,theta,rate_vec_fun,...
                Q,y0_cell,f_indices,xdata);
            
ssfun = @(theta,data) sum(sum((data.ydata-modelfun(data.xdata,theta)).^2));

% specify parameters and set initial values
params = {};
iter = 1;
for p = fit_indices
    if iter == 1
        params = {{char(reaction_parameter_index(p)), rand()*true_param_vec(p), lb_vec(p), ub_vec(p)}};
    else
        params = {params{:} {char(reaction_parameter_index(p)), rand()*true_param_vec(p), lb_vec(p), ub_vec(p)}};    
    end
    iter = iter + 1;
end    

model.ssfun  = ssfun;
model.sigma2 = 1e-6; 

options.nsimu = 1000;
options.updatesigma = 0;
% options.qcov = tcov; % covariance from the initial fit

[results,chain,s2chain] = mcmcrun(model,data,params,options);

%%
close all

figure(2); clf
mcmcplot(chain,[],results,'chainpanel');

mcmc_sol = results.mean;
mcmc_curves = modelfun(time_exp_array,mcmc_sol);

figure;
hold on
plot(fluo_exp_array)
plot(mcmc_curves)

figure; clf
mcmcplot(chain,[],results,'chainpanel');

figure; clf
mcmcplot(chain,[],results,'pairs');
