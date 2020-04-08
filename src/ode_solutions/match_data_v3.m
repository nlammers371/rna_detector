clear
close all
addpath('../utilities')

DataPath = '../../out/ode_studies/';
mkdir(DataPath)

% specify project to load
project = 'ncr_general';
% load 
load([DataPath project '_setup.mat'])

% make figure path
FigPath = ['../../fig/ode_studies/' project '/' ];
mkdir(FigPath)

t_max = 2e3;

% Filter for key pathway interactions only. Include Activator cleavage

% filters for catalytic and specific interactions 
cat_spec_filter = catalytic_activity_flag_vec | specific_activity_flag_vec; 
% filter to remove catalytic reactions that involve guide cleavage
cat_g_filter = ~contains(reaction_string_cell,{'G1Z','G2Z'});
% filter for all reactions for which guid RNA is an input or product
cg_indices = find(strcmp(full_atom_cell,'G1')|strcmp(full_atom_cell,'G2')|...
    strcmp(full_atom_cell,'C1')|strcmp(full_atom_cell,'C2'));
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

%% initial reactant concentrations
S0 = 200;
C1_G10 = 50;
C2_G20 = .8;
T10 = 0;
T20 = 0;
T2I0 = 20;
% vectors for substitutions
init_sym_array = {'T1','T2','S' 'T2I' 'G1:C1' 'G2:C2'};
init_val_vec = [T10 T20 S0 T2I0 C1_G10 C2_G20];
% initialize concentration vector. Everything not explicitly assigned taken
% to be zero
y0_vec = zeros(size(reactant_cell_truncated));
for i = 1:numel(init_sym_array)
    species = init_sym_array{i};
    spec_indices = find(strcmp(reactant_cell_truncated,species));
    y0_vec(spec_indices) = init_val_vec(i);
end

% define basic on/off and cleavaeg rates
kon = 0.25;
koff_s = 1;
kcat_cas = 200*2e-5; % basal cas cleavage rate
kcat_ratio = 1/2e-5;
koff_ns = 1.5e3; 
% define substrate cleavage susceptibilities
sS = 1; % reporter
sT1 = 1; % activator 1
sT2 = 1; % activator 2
sT2I = 1; % caged activator

% rate vector to substitute into
rate_vec_val = rate_vec_zero_order_truncated;

% generate valued long rate vec
for r = 1:numel(reaction_unit_index)
    rr = reaction_unit_index(r);    
    
    rate_vec_val = subs(rate_vec_val,rr,eval(rr));
    
end
rate_vec_val = double(rate_vec_val);
%% solve ODEs numerically
close all
f_indices = find(contains(reactant_cell_truncated,'SZ'));
ti_indices = find(contains(reactant_cell_truncated,'T2I'));
t2_indices = find(contains(reactant_cell_truncated,'T2')...
    &~contains(reactant_cell_truncated,'T2I')...
    &~contains(reactant_cell_truncated,'T2Z'));
t2z_indices = find(contains(reactant_cell_truncated,'T2Z'));

[t_pos,y_pos] = ode15s(@(t,y) ncr_solver(t,y,rate_vec_val,stoich_mat_truncated),[0 t_max],y0_vec);


figure;
plot(t_pos,sum(y_pos(:,ti_indices),2))
title('caged activator 2')

figure;
plot(t_pos,sum(y_pos(:,t2_indices),2))
title('free activator 2')

figure;
plot(t_pos,sum(y_pos(:,t2z_indices),2))
title('cleaved activator 2')

figure;
plot(t_pos,sum(y_pos(:,f_indices),2))
title('reporter')

%% attempt to fit rate parameters

% specify subset of parameters to fit
fit_indices = [1 2 3 4 7:9];
fixed_indices = find(~ismember(1:numel(reaction_unit_index),fit_indices));
rate_vec_fit = rate_vec_zero_order_truncated;
for r = fixed_indices
    rr = reaction_unit_index(r);        
    rate_vec_fit = subs(rate_vec_fit,rr,eval(rr));    
end

rate_vec_fun = matlabFunction(rate_vec_fit);

% generate fake experimental data
time_exp = t_pos;
f_vec = sum(y_pos(:,f_indices),2);
fluo_exp = f_vec;%interp1(t_pos,f_vec,time_exp);%sum(y_pos(:,f_indices),2);%interp1([0 t_max],[0,S0*.4],time_exp);
rate_params_init = rand(size(fit_indices));

fit_fun = @(rate_params,time_exp) ode_objfunction(t_max,rate_params,rate_vec_fun,...
    stoich_mat_truncated,y0_vec,f_indices,time_exp);
%%
options = optimoptions('lsqcurvefit','MaxIterations',50);
kon = 0.25;
koff_s = 1;
kcat_cas = 200*2e-5; % basal cas cleavage rate
kcat_ratio = 1/2e-5;
koff_ns = 1.5e3; 
% define substrate cleavage susceptibilities
sS = 1; % reporter
sT1 = 1; % activator 1
sT2 = 1; % activator 2
sT2I = 1; % caged activator

% set bounds for inference parameters
ub = [1    1e7 10  5e3 10 10 10];
lb = [2e-4 1e3 0.1 100 0 0 0];


tic
rate_vec_out=lsqcurvefit(fit_fun,[.1 1e5 1 1e3 1 1 1],time_exp,fluo_exp,lb,ub,options);
toc
%%
initial_sol = fit_fun(rate_params_init,time_exp);
final_sol = fit_fun(rate_vec_out,time_exp);
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