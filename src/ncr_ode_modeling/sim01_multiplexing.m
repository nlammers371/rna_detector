% simulate impact of multiplexing
clear 
close all

% basic path info
addpath(genpath('../utilities'))
DataPath = '../../out/ncr_ode_modeling/';
mkdir(DataPath)

% specify project to load
project = 'ncr_v2';
% load 
load([DataPath project '.mat'])

% make figure path
FigPath = ['../../fig/ncr_ode_modeling/' project '/' ];
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
% specify subset of parameters to fit


% define "true" rates to fit
b = 1e-6; % Cas13 1/SNR
k = 0.025; % association rate (s^-1 nM^-1)
kc = 900; % Cas13 catalytic rate
kd_cga = 1e-8;
rcg = 5e-3; % off rate for Cas13:gRNA
rga = 1e-2;
ria = 1e-2; % off rate for CleavedCage:Activator
rns = 1e3; % off rate for all nonspecific interactions
keff = 1e-4;
kact = k*keff / (k - keff);

%
t_max = 3600; % duration of time to solve for (need to change this to calculate dynamically)
t_res = 0.1;
A10 = 0.2; % concentration of target RNA
IA20 = 20; % caged secondary activator
S0 = 200; % reporter substrate
RNP1 = 20; % Cas13:gRNA1
RNP2 = 0.8; % Cas13:gRNA2

% calculate fraction of cas/guide initially bound
CG10 = (rcg + 2*k*RNP1 - sqrt(rcg)* sqrt(rcg + 4*k*RNP1))/(2*k);
C10 = RNP1 - CG10;
G10 = RNP1 - CG10;

CG20 = (rcg + 2*k*RNP2 - sqrt(rcg)* sqrt(rcg + 4*k*RNP2))/(2*k);
C20 = RNP2 - CG20;
G20 = RNP2 - CG20;

%%%% specify initial reactant concentrations 
% vectors for substitutions
init_sym_array =      {'A1'  'S' 'IA2' 'C13'       'G1'  'G2'    'C13:G1'   'C13:G2'};
init_val_vec_cage =    [0    S0  IA20   C10+C20    G10   G20        CG10      CG20]; % cage only
init_val_vec_full =    [A10  S0  IA20   C10+C20    G10   G20        CG10      CG20]; % primary only
init_val_vec_primary = [A10  S0  0      C10+C20    G10   G20        CG10      CG20]; % primary only
init_val_vec_neg =     [0    S0  0      C10+C20    G10   G20        CG10      CG20]; % no activator

% noise term
sim_noise = 5e-3*S0;

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
non_ncr_indices = [2 4];
name_cell = {'NCR', 'Primary Only', 'Cage Only', 'No Activator'};


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

% generate interpolated vectors for fitting
time_exp = 0:0.1:t_max;
time_exp_array = repmat(time_exp',1,numel(y0_cell));

fluo_exp_array = NaN(numel(time_exp),numel(y0_cell));
for f = 1:numel(fluo_cell_raw)
    fluo_exp_array(:,f) = interp1(time_cell_raw{f},fluo_cell_raw{f},time_exp) + normrnd(0,sim_noise,1,numel(time_exp));
end

figure;
plot(time_exp/60,fluo_exp_array)