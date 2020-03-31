% Script to illustrate the sensitivity issue with naive NCR

% Even "negative" samples with low level background activity will trigger
% the cascade, and the timing with which the signals for positve and
% negative samples is quite close, even given relatively large initial
% differences

clear
close all
addpath('../utilities')

FigPath = '../../fig/ode_studies/simple_model/';
mkdir(FigPath)
DataPath = '../out/ode_studies/';
mkdir(DataPath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reactant key (first order)
%   species 1: A13I (caged activator)
%   species 2: A13 (free activator and/or virus)
%   species 3: C13 (free Cas13)
%   species 4: S (dark reporter)
%   species 5: F (Cleaved (fluorescent) reporter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize rates
rate_vec = sym({'kon', 'koff_ns', 'koff_s', 'kcat_high', 'kcat_low'});
syms kon koff_ns koff_s kcat_high kcat_low

% generate list of raction atoms 
atom_cell = {'C13' 'A13' 'A13I' 'S' 'F'};
atom_list = sym(atom_cell);
n_atoms = numel(atom_list);

% compound reaction components
compound_cell = {'C13_A13','C13_A13I','C13_S','C13_A13__A13I','C13_A13__S'};
off_rate_vec = [koff_s koff_ns koff_ns koff_ns koff_ns];
compound_list = sym(compound_cell);
n_compounds = numel(compound_cell); 

% define compounds with catalytic activity
cat_compounds =  {'C13_A13I','C13_S','C13_A13__A13I','C13_A13__S'};
cat_rates = [kcat_low kcat_low kcat_high kcat_high];
cat_products = {'A13','F','A13','F'};

% initialize stoichiometry matrix 
full_reactant_cell = [atom_cell compound_cell];
full_reactant_list = [atom_list compound_list];
n_reactants = n_atoms + n_compounds; % total number of ractions
n_reactions = 2*n_compounds + numel(cat_compounds);
Q_mat = zeros(n_reactants,n_reactions);
         
rate_vec_first_order = sym('NA',[1,n_reactions]);
rate_cell_fun = {};
cat_rate_indices = [];
cat_cp_indices = [];
rate_ind = 1;
% first do association and disassociation reactions
for i = 1:numel(compound_cell)
    % extract compound
    cp_str = compound_cell{i};
    cp_ind = i + n_atoms;
    % extract reactant components
    cp3_flag = contains(cp_str,'__');
    if cp3_flag
        reactants = strsplit(cp_str,'__');
    else
        reactants = strsplit(cp_str,'_');
    end
    r_ind1 = find(strcmp(full_reactant_cell,reactants(1)));
    r_ind2 = find(strcmp(full_reactant_cell,reactants(2)));
    

    % add association reaction and update Q
    rate_vec_sym(rate_ind) = full_reactant_list(r_ind1)*full_reactant_list(r_ind2)*kon;
    rate_vec_fun{rate_ind} = @(y,r) y(r_ind1)*y(r_ind2)*kon;
    rate_vec_first_order(rate_ind) = kon;
    Q_mat([r_ind1 r_ind2],rate_ind) = -1;
    Q_mat(cp_ind,rate_ind) = 1;

    
    %%%%%%%%%%% increment 
    rate_ind = rate_ind + 1;

    % add disassociation reaction
    off_rate = off_rate_vec(i);
    
    rate_vec_sym(rate_ind) = full_reactant_list(cp_ind)*off_rate;
    rate_vec_fun{rate_ind} = @(y) y(cp_ind) * off_rate;
    rate_vec_first_order(rate_ind) = off_rate;

    Q_mat([r_ind1 r_ind2],rate_ind) = 1;
    Q_mat(cp_ind,rate_ind) = -1;

    % increment 
    rate_ind = rate_ind + 1;   
    
    % check for catalytic activity 
    cat_ind = find(strcmp(cat_compounds,cp_str));
    if ~isempty(cat_ind)

        % update rate vec
        rate_vec_sym(rate_ind) = full_reactant_list(cp_ind)*cat_rates(cat_ind);
        rate_vec_fun{rate_ind} = @(y) y(cp_ind) * cat_rates(cat_ind);
        rate_vec_first_order(rate_ind) = cat_rates(cat_ind);
        % update Q
        product_str = cat_products{cat_ind};
        product_ind = find(strcmp(full_reactant_cell,product_str)); % catalyzed product
        
        Q_mat([r_ind1 product_ind],rate_ind) = 1;
        Q_mat(cp_ind,rate_ind) = -1;

        % increment
        cat_rate_indices = [cat_rate_indices rate_ind];
        cat_cp_indices = [cat_cp_indices cp_ind];
        rate_ind = rate_ind + 1;            
    end   
end


%% (1) make initial plots to illustrate the promise and challenge
close all

% basic parameters
t_max = 360;
atom_cell = {'C13' 'A13' 'A13I' 'S' 'F'};
f_ind = find(strcmp(atom_cell,'F'));
S0 = 200;
AI0 = 200;
C0 = 200;
R0 = 1e-4;
% positive and negative samples
y0_pos_vec = zeros(1,size(Q_mat,1));
y0_pos_vec(1:4) = [C0 R0 AI0 S0];
y0_neg_vec = y0_pos_vec;
y0_neg_vec(2) = 0;
y0_pos_vec_sh = y0_pos_vec;
y0_pos_vec_sh(3) = 0;

% idealized and non idealized rate vectors
rate_vec_ideal = [1 1 1 0.3 0];
rate_vec_actual = [1 1 1 0.1 0.1*2e-4];
p_vec_ideal = NaN(size(rate_vec_first_order));
p_vec_actual = NaN(size(rate_vec_first_order));

% generate valued long rate vec
for r = 1:numel(rate_vec_ideal)
    rr = rate_vec(r);
    ri = find(rr==rate_vec_first_order);
%     rep_indices = rate_vec_first_order==rep_rate;
    p_vec_ideal(ri) = double(subs(rr,rr,rate_vec_ideal(r))); 
    p_vec_actual(ri) = double(subs(rr,rr,rate_vec_actual(r))); 
end

% solve odes numerically
[t_pos_ideal ,y_pos_ideal] = ode15s(@(t,y) ncr_solver(t,y,p_vec_ideal,Q_mat),[0 t_max],y0_pos_vec);
[t_pos_ideal_sh ,y_pos_ideal_sh] = ode15s(@(t,y) ncr_solver(t,y,p_vec_ideal,Q_mat),[0 t_max],y0_pos_vec_sh);
[t_neg_ideal ,y_neg_ideal] = ode15s(@(t,y) ncr_solver(t,y,p_vec_ideal,Q_mat),[0 t_max],y0_neg_vec);

[t_pos_actual ,y_pos_actual] = ode15s(@(t,y) ncr_solver(t,y,p_vec_actual,Q_mat),[0 t_max],y0_pos_vec);
[t_neg_actual,y_neg_actual] = ode15s(@(t,y) ncr_solver(t,y,p_vec_actual,Q_mat),[0 t_max],y0_neg_vec);

detection_threshold = repelem(50,numel(t_pos_ideal));
% "Sherlock" only plot
ideal_fig_sh = figure;
cmap1 = brewermap(9,'Set2');
hold on
plot(t_pos_ideal_sh/60,y_pos_ideal_sh(:,f_ind)/S0*100,'-','Color',cmap1(1,:),'LineWidth',1.5);
plot(t_pos_ideal/60,detection_threshold,'--','Color','black')
ylim([0 110])
xlim([0 t_max/60])
xlabel('time (minutes)')
ylabel('fluorescent signal (% of max)')
set(gca,'FontSize',14)
saveas(ideal_fig_sh,[FigPath 'illustrative_sher_ideal.png'])

% add in NCR (idealized)
ideal_fig_sh = figure;
hold on
p1 = plot(t_pos_ideal/60,y_pos_ideal(:,f_ind)/S0*100,'-','Color',cmap1(2,:),'LineWidth',1.5);
p2 = plot(t_pos_ideal_sh/60,y_pos_ideal_sh(:,f_ind)/S0*100,'-','Color',cmap1(1,:),'LineWidth',1.5);
plot(t_pos_ideal/60,detection_threshold,'--','Color','black')
ylim([0 110])
xlim([0 t_max/60])
xlabel('time (minutes)')
ylabel('fluorescent signal (% of max)')
legend([p1 p2],'with NCR','without NCR','Location','southeast')
set(gca,'FontSize',14)
saveas(ideal_fig_sh,[FigPath 'illustrative_ideal.png'])
% set(gca,'Yscale','log')
% grid on
% saveas(ideal_fig_sh,[FigPath 'illustrative_ideal_ylog.png'])

% NCR (actual, positive)
actual_fig_pos = figure;
hold on
p1 = plot(t_pos_actual/60,y_pos_actual(:,f_ind)/S0*100,'-','Color',cmap1(2,:),'LineWidth',1.5);
plot(t_pos_ideal/60,detection_threshold,'--','Color','black')
ylim([0 110])
xlim([0 t_max/60])
xlabel('time (minutes)')
ylabel('fluorescent signal (% of max)')
% legend([p1 p2],'with NCR','without NCR','Location','southeast')
set(gca,'FontSize',14)
saveas(actual_fig_pos,[FigPath 'illustrative_actual.png'])
% set(gca,'Yscale','log')
% grid on
% saveas(actual_fig_pos,[FigPath 'illustrative_actual_ylog.png'])
%
actual_fig_pos_neg = figure;
hold on
p1 = plot(t_pos_actual/60,y_pos_actual(:,f_ind)/S0*100,'-','Color',cmap1(2,:),'LineWidth',1.5);
p2 = plot(t_neg_actual/60,y_neg_actual(:,f_ind)/S0*100,'--','Color',cmap1(3,:),'LineWidth',1.5);
plot(t_pos_ideal/60,detection_threshold,'--','Color','black')
ylim([0 110])
xlim([0 t_max/60])
xlabel('time (minutes)')
ylabel('fluorescent signal (% of max)')
legend([p1 p2],'positive sample','negative sample','Location','southeast')
set(gca,'FontSize',14)
saveas(actual_fig_pos,[FigPath 'illustrative_actual_pos_neg.png'])

%%* Can we estimate contribution from Cas13 bound to virus
close all
RC = 60*rate_vec_actual(4)*R0*(y_pos_actual(:,6) ./ (y_pos_actual(:,2) + y_pos_actual(:,6)));
CC = 60*rate_vec_actual(5)*y_pos_actual(:,1);
AC = 60*rate_vec_actual(4)*(y_pos_actual(:,6) - R0*(y_pos_actual(:,6) ./ (y_pos_actual(:,2) + y_pos_actual(:,6))));
AC(t_pos_actual<=8.8e-6) = 0;
cat_cap_fig = figure;
hold on
area(t_pos_actual/60,CC,'FaceColor',cmap1(3,:),'FaceAlpha',0.75)
area(t_pos_actual/60,RC,'FaceColor',cmap1(2,:),'FaceAlpha',0.75)
area(t_pos_actual/60,AC,'FaceColor',cmap1(5,:),'FaceAlpha',1)
legend('Cas13 (free)','Cas13 (bound to target)','Cas13 (uncaged)','Location','northeast')
xlabel('time (minutes)')
ylabel('catalytic potential (nM^{-1}min^{-1})')
xlim([0 t_max/60])
ylim([0 1.3e3])
grid on
set(gca,'Fontsize',14)
saveas(cat_cap_fig,[FigPath 'catalytic_capacities.png'])
set(gca,'XScale','log','YScale','log')
saveas(cat_cap_fig,[FigPath 'catalytic_capacities_log.png'])

%% Conduct parameter sweep
% ALPHA: ratio of initial Cas13 and RNA concentrations
% R: RNA concentration

% set basic param values
n_iter = 1e2;
t_max = 1e20;
atom_cell = {'C13' 'A13' 'A13I' 'S' 'F'};
f_ind = find(strcmp(atom_cell,'F'));
S0 = 200;
AI0 = 200;
R0_bounds = [-8,0];
alpha_bounds = [1,10];
rate_vec_sim = [1 1 1 0.1 0.1*2e-4];
y0_base = zeros(1,size(Q_mat,1));
y0_base(1:4) = [0 0 AI0 S0];

% set rate vec
p_vec_sim = NaN(size(rate_vec_first_order));

% generate valued long rate vec
for r = 1:numel(rate_vec_ideal)
    rr = rate_vec(r);
    ri = find(rr==rate_vec_first_order);
    p_vec_sim(ri) = double(subs(rr,rr,rate_vec_sim(r))); 
end

% initialize vectors to store results
sim_struct = struct;
% iterate through conditions
for iter = 1:n_iter
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

%%
