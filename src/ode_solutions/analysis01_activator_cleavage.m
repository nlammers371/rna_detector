% script to determine if activator cleavage by Cas13 can replicate key
% features in NCR experimental data
clear
close all
addpath('../utilities')

DataPath = '../../out/ode_studies/';
mkdir(DataPath)

% specify project to load
project = 'ncr_general';
sub_project = 'activator_cleavage';
% load 
load([DataPath project '_setup.mat'])

% make figure path
FigPath = ['../../fig/ode_studies/' project '/' sub_project '/' ];
mkdir(FigPath)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Specify model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter for key pathway interactions only. Include Activator cleavage
% filters for catalytic and specific interactions 
cat_spec_filter = catalytic_activity_flag_vec | specific_activity_flag_vec; 
% filter to remove catalytic reactions that involve guide cleavage
cat_g_filter = ~contains(reaction_string_cell,{'G1Z','G2Z'});
% filter for all reactions for which guid RNA is an input or product
cg_indices = find(strcmp(full_atom_cell,'G1')|strcmp(full_atom_cell,'G2')|...
    strcmp(full_atom_cell,'C1')|strcmp(full_atom_cell,'C2')|...
    strcmp(full_atom_cell,'T1Z')|strcmp(full_atom_cell,'T1'));
free_g_filter = max(stoichiometry_matrix(cg_indices,:)~=0)==0;

final_filter = free_g_filter & cat_g_filter & cat_spec_filter;

% apply filters and generate truncated versions of relevant arrays

% remove undesired reactions
stoich_mat_truncated = stoichiometry_matrix(:,final_filter);
reaction_cell_truncated = reaction_string_cell(final_filter);
rate_vec_zero_order_truncated = rate_vec_zero_order_sym(final_filter);

% remove product entities not involved in final set of reactions
obs_flag = all(stoich_mat_truncated==0,2);
reactant_cell_truncated = full_reactant_cell(~obs_flag);
stoich_mat_truncated = stoich_mat_truncated(~obs_flag,:);

% specifiy reaction rate values and input concentrations
reaction_rate_index = unique(rate_vec_zero_order_truncated);
% decompose rates into fundamental units
reaction_units = [];
for i = 1:numel(reaction_rate_index)
    r_sym = reaction_rate_index(i);
    t_char = char(r_sym);
    char_list = strsplit(t_char,'*');
    for c = 1:numel(char_list)
        reaction_units = [reaction_units eval(char_list{c})];
    end
end
reaction_unit_index = unique(reaction_units);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% set key parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_max = 3600; % duration of time to solve for

%%%% specify initial reactant concentrations 
% vectors for substitutions
init_sym_array = {'T2','S' 'T2I' 'G1:C1' 'G2:C2'};
init_val_vec_ncr =    [0    200 20 50 0.8]; % for simulated negative NCR raction
init_val_vec_target = [3e-2 200  0 0  0.8]; % for target "sub critical" reaction

% initialize concentration vector. Everything not explicitly assigned taken
% to be zero
y0_vec_ncr = zeros(size(reactant_cell_truncated));
y0_vec_target = zeros(size(reactant_cell_truncated));
for i = 1:numel(init_sym_array)
    species = init_sym_array{i};
    spec_indices = find(strcmp(reactant_cell_truncated,species));
    y0_vec_ncr(spec_indices) = init_val_vec_ncr(i);
    y0_vec_target(spec_indices) = init_val_vec_target(i);
end
y0_cell = {y0_vec_ncr};
%%%%%%%%%%%%%%%%%
%%% Specify which rate parameters to fit and specify values for "given"
%%% parameters

% specify subset of parameters to fit
fit_indices = [1:8];
% define basic on/off and cleavaeg rates
kon = 0.25;
koff_s = 1;
kcat_cas = 200*2e-7; % basal cas cleavage rate
kcat_ratio = 1/2e-7;
koff_ns = 1.5e3; 
% define substrate cleavage susceptibilities (these are the key parameters
% here)
sS = 1; % reporter
sT2 = 1e-4; % activator 2
sT2I = 1e-1; % caged activator


% rate vector to substitute into
rate_vec_target = rate_vec_zero_order_truncated;
% generate valued long rate vec
for r = 1:numel(reaction_unit_index)
    rr = reaction_unit_index(r);        
    rate_vec_target = subs(rate_vec_target,rr,eval(rr));    
end
rate_vec_target = double(rate_vec_target);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% solve ODEs for target numerically %%%%%%%%%%%%%%%%%%
% 
close all
f_indices = find(contains(reactant_cell_truncated,'SZ'));
[t_target,y_target] = ode15s(@(t,y) ncr_solver(t,y,rate_vec_target,stoich_mat_truncated),[0 t_max],y0_vec_target);
f_target = sum(y_target(:,f_indices),2);
[t_ncr_null,y_ncr_null] = ode15s(@(t,y) ncr_solver(t,y,rate_vec_target,stoich_mat_truncated),[0 t_max],y0_vec_ncr);
f_ncr_null = sum(y_ncr_null(:,f_indices),2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% attempt to fit rate parameters %%%%%%%%%%%%%%%%%%

% substitute "given" values into rate vector
fixed_indices = find(~ismember(1:numel(reaction_unit_index),fit_indices));
rate_vec_fit = rate_vec_zero_order_truncated;
for r = fixed_indices
    rr = reaction_unit_index(r);        
    rate_vec_fit = subs(rate_vec_fit,rr,eval(rr));    
end
% convert this to a matlabFunction
rate_vec_fun = matlabFunction(rate_vec_fit);

% now set bounds for variables to be inferred
ub_vec = NaN(size(fit_indices));
lb_vec = NaN(size(fit_indices));
rate_params_init = rand(size(fit_indices));
iter = 1;
for r = fit_indices
    rate_init = eval(reaction_unit_index(r));        
    ub_vec(iter) = 1e4*rate_init;
    lb_vec(iter) = 1e-4*rate_init;
    rate_params_init(iter) = eval(reaction_unit_index(r))*rand()*1e1+lb_vec(iter);
    iter = iter + 1;
end

% define fitting function
fit_fun = @(rate_params,time_exp_array) ode_objfunction(t_max,rate_params,rate_vec_fun,...
    stoich_mat_truncated,y0_cell,f_indices,time_exp_array);

% set optimization options
options = optimoptions('lsqcurvefit','MaxIterations',200);

% generate interpolated vectors for fitting
time_exp = 0:0.1:t_max;
time_exp_array = [time_exp'];
fluo_exp_array(:,1) = interp1(t_target,f_target,time_exp);
% fluo_exp_array(:,2) = fluo_exp_array(:,1);
% fluo_exp_array(:,2) = interp1(t_ncr_null,f_ncr_null,time_exp);

% perform optimization
tic
rate_vec_out=lsqcurvefit(fit_fun,rate_params_init,time_exp_array,fluo_exp_array,lb_vec,ub_vec,options);
toc

initial_sol = fit_fun(rate_params_init,time_exp_array);
final_sol = fit_fun(rate_vec_out,time_exp_array);
% %%
% fluo_factor = 12000;
% close all
% 
% % under updated background assumption
% match_fig_3 = figure;
% cmap2 = brewermap([],'set2');
% hold on
% p1 = plot(t_pos_high,y_pos_high(:,f_indices)/S0*fluo_factor,'-','Color','black','LineWidth',1.5);
% p2 = plot(t_pos_low,y_pos_low(:,f_indices)/S0*fluo_factor,'--','Color','black','LineWidth',1.5);
% 
% ylim([0 fluo_factor*1.05])
% xlim([0 2000])
% 
% xlabel('time (seconds)')
% ylabel('fluorescent signal (au)')
% legend([p1 p2],'0.01 nM activator (NCR)','0.001 nM activator (NCR)','Location','southeast')
% 
% set(gca,'FontSize',14)
% grid on
% saveas(match_fig_3,[FigPath 'cas13_ncr_pos_only_sl' num2str(SI_ratio) '.png'])
% 
% 
% match_fig_4 = figure;
% hold on
% p1 = plot(t_pos_high,y_pos_high(:,f_indices)/S0*fluo_factor,'-','Color','black','LineWidth',1.5);
% p3 = plot(t_neg,y_neg(:,f_indices)/S0*fluo_factor,'-','Color','red','LineWidth',1.5);
% p2 = plot(t_pos_low,y_pos_low(:,f_indices)/S0*fluo_factor,'--','Color','black','LineWidth',1.5);
% 
% ylim([0 fluo_factor*1.05])
% xlim([0 2000])
% 
% xlabel('time (seconds)')
% ylabel('fluorescent signal (au)')
% legend([p1 p2 p3],'0.01 nM activator','0.001 nM activator','Cas/G only (5e6 SNR)','Location','southeast')
% 
% set(gca,'FontSize',14)
% grid on
% saveas(match_fig_4,[FigPath 'cas13_ncr_pos_only_neg_sl' num2str(SI_ratio) '.png'])
% 
% 
% match_fig_5 = figure;
% hold on
% p1 = plot(t_pos_high,y_pos_high(:,f_indices)/S0*fluo_factor,'-','Color','black','LineWidth',1.5);
% % p3 = plot(t_neg,y_neg(:,f_ind)/S0*fluo_factor,'-','Color',[.5 .5 .5],'LineWidth',1.5);
% p2 = plot(t_pos_low,y_pos_low(:,f_indices)/S0*fluo_factor,'--','Color','black','LineWidth',1.5);
% p4 = plot(t_neg_low,y_neg_low(:,f_indices)/S0*fluo_factor,'--','Color','red','LineWidth',1.5);
% 
% ylim([0 fluo_factor*1.05])
% xlim([0 t_max])
% 
% xlabel('time (seconds)')
% ylabel('fluorescent signal (au)')
% legend([p1 p2 p4],'0.01 nM activator','0.001 nM activator','Cas/G only (5e8 SNR)','Location','southeast')
% 
% set(gca,'FontSize',14)
% grid on
% saveas(match_fig_5,[FigPath 'cas13_ncr_pos_only_neg_low_sl' num2str(SI_ratio) '.png'])