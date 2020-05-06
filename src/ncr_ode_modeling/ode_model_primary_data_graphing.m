% Script to use numerical ODE solver and nonlinear least square fitting to
% estimate reaction parameters from experimental data
clear 
close all

% basic path info
addpath('../utilities')
DataPath = '../../out/ncr_ode_modeling/';
mkdir(DataPath)

% specify project to load
project = 'primary_only_v3';

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


%%%% specify initial reactant concentrations 
% vectors for substitutions

init_sym_array =     {'A1' 'S' 'C13' 'G1'};
points = 10;
% As = linspace(0, 1, points);
As = [0, 0.0000002, 0.000002, 0.00002, 0.0002, 0.002, 0.02, 0.2, 2, 20];
Ss = repelem(S0, points);
RNP1s = repelem(RNP1, points);
y0_vec = squeeze(zeros([size(full_reactant_list),points]));

spec_indices = find(strcmp(full_reactant_list,'A1'));
y0_vec(spec_indices, :) = As;
spec_indices = find(strcmp(full_reactant_list,'S'));
y0_vec(spec_indices,:) = Ss;
spec_indices = find(strcmp(full_reactant_list,'C13'));
y0_vec(spec_indices,:) = RNP1s;
spec_indices = find(strcmp(full_reactant_list,'G1'));
y0_vec(spec_indices,:) = RNP1s;

y0_cell=num2cell(y0_vec,1)'

%{
init_val_vec_primary =   [A10  S0  RNP1 RNP1]; % primary only
init_val_vec_neg =    [0  S0  RNP1 RNP1]; % no activator

% initialize concentration vector. Everything not explicitly assigned taken
% to be zero
init_sym_array =     {'A1' 'S' 'C13' 'G1'};
y0_vec_primary = zeros(size(full_reactant_list));
y0_vec_neg = zeros(size(full_reactant_list));
for i = 1:numel(init_sym_array)
    species = init_sym_array{i};
    spec_indices = find(strcmp(full_reactant_list,species));      
    y0_vec_primary(spec_indices) = init_val_vec_primary(i);    
    y0_vec_neg(spec_indices) = init_val_vec_neg(i);
end

y0_cell = {y0_vec_primary, y0_vec_neg}; % all conditions in this cell will be included in fit routine
y0_cell
%}

%%%%%%%%%%%%%%%%%
%%% Specify which rate parameters to fit and specify values for "given"
%%% parameters
% return
% specify subset of parameters to fit
fit_indices = [1,2,3,4,5,6,7,8];
n_inference_runs = 20;
sim_noise = 1e-2*S0; % gaussian noise to add to simulated data

% define initial guess for parameters to fit
b = 7.4162e-06; % Cas13 1/SNR
k =  0.012553; % association rate (s^-1 nM^-1)
kc = 229.98; % Cas13 catalytic rate 
kd_cga =  1.2136e-12; %dissassociation constant for cas13:guide and activatior 
rcg =  0.19425; % off rate for Cas13:gRNA      
rcga = 0.0019832; %on rate for cas13:guide and activatior      
rga =  0.0061094; % off rate for gRNA:Activator
rns = 137.67; % off rate for all nonspecific interactions

% get vector of guessed param values
true_param_vec = eval(reaction_parameter_index);

% specify degree of uncertainty for each param value (need to play with
% this)
sigma_vec = [1e3, 1e1, 1e1, 1e1, 1e1, 1e5, 1e1, 1e1];

%%%%%%%%%%% Generate "experimental data" using true param values %%%%%%%%%%


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
hold on
fluo_cell_raw = cell(1,length(y0_cell));
time_cell_raw = cell(1,length(y0_cell));

for f = 1:length(y0_cell)
    y0_cell{f}
    [t_fit,y_fit] = ode15s(@(t,y) ncr_solver_v2(t,y,rate_vec_val,Q),[0 t_max],y0_cell{f});
    plot(t_fit,(3500/200)*sum(y_fit(:,f_indices),2));
    fluo_cell_raw{f} = sum(y_fit(:,f_indices),2);
    time_cell_raw{f} = t_fit;
end
fluo_cell_raw


% p = scatter(time_cell_raw,fluo_cell_raw);
% q = plot(time_exp,fluo_fit_array(:,:,best_fit_index),'Color','black','LineWidth',1.5);
% legend([p(1) p(2) q(1)],'primary only','no activator', 'model fits','Location','northwest')
% 
% figure(2);
% p = plot(time_exp,fluo_exp_array-fluo_fit_array(:,:,best_fit_index));
% legend([p(1) p(2)],'primary only','no activator','Location','northwest');