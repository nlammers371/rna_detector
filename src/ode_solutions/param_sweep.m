clear 
close all

% basic path info
addpath('../utilities')
DataPath = '../../out/ode_studies/';
mkdir(DataPath)

% specify project to load
project = 'naive_ncr';
% load 
load([DataPath project '_setup.mat'])

% make figure path
FigPath = ['../../fig/ode_studies/' project '/' ];
mkdir(FigPath)

%% Conduct parameter sweep
% ALPHA: ratio of initial Cas13 and RNA concentrations
% R: RNA concentration

% set basic param values
n_iter = 1e6;
t_max = 1e20;
atom_cell = {'C13' 'A13' 'A13I' 'S' 'F'};
f_ind = find(strcmp(atom_cell,'F'));
S0 = 50;
AI0 = 50;
R0_bounds = [-8,0];
alpha_bounds = [1,10];
rate_vec_sim = [1 1 1 200 200*2e-4];
y0_base = zeros(1,size(Q_mat,1));
y0_base(1:4) = [0 0 AI0 S0];

% set rate vec
p_vec_sim = NaN(size(rate_vec_first_order));

% generate valued long rate vec
for r = 1:numel(rate_vec_sim)
    rr = rate_vec(r);
    ri = find(rr==rate_vec_first_order);
    p_vec_sim(ri) = double(subs(rr,rr,rate_vec_sim(r))); 
end

% initialize vectors to store results
sim_struct = struct;
% iterate through conditions
parfor iter = 1:n_iter
    % set initial conditions
    y0_pos_iter = y0_base;
    y0_neg_iter = y0_base;        
    % randomly draw R0 and C0 
    R0_iter = 10^(rand()*(R0_bounds(2) - R0_bounds(1)) + R0_bounds(1));
    alpha_iter = 10^(rand()*(alpha_bounds(2) - alpha_bounds(1)) + alpha_bounds(1));
    C0_iter = alpha_iter*R0_iter;
    y0_pos_iter(1:2) = [C0_iter R0_iter];     
    y0_neg_iter(1:2) = [C0_iter 0];
    % solve ODEs
    Opt = odeset('Events', @myEvent);
    [t_pos,y_pos] = ode15s(@(t,y) ncr_solver(t,y,p_vec_sim,Q_mat),[0 t_max],y0_pos_iter,Opt);
    [t_neg,y_neg] = ode15s(@(t,y) ncr_solver(t,y,p_vec_sim,Q_mat),[0 t_max],y0_neg_iter,Opt);

    % store results
    sim_struct(iter).R0 = R0_iter;
    sim_struct(iter).C0 = C0_iter;
    sim_struct(iter).t_pos = t_pos;
    sim_struct(iter).t_neg = t_neg;
    sim_struct(iter).y_pos = y_pos;
    sim_struct(iter).y_neg = y_neg;
    sim_struct(iter).t2_pos = max(t_pos);
    sim_struct(iter).t2_neg = max(t_neg);        
end

save([DataPath 'param_sweep.mat'],'sim_struct')
obs_time_vec = [sim_struct.t2_pos];
time_ratio_vec = [sim_struct.t2_neg]./[sim_struct.t2_pos];
r0_vec = [sim_struct.R0];

%%


performance_fig = figure;
scatter(obs_time_vec(r0_vec<1e-4),time_ratio_vec(r0_vec<1e-4),50,log10(r0_vec(r0_vec<1e-4)),'filled','MarkerFaceAlpha',0.2)
set(gca,'Xscale','log','Yscale','log')
ylim([1 10])
colorbar
