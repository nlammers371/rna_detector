% Script to simulate early exponential phase of basic NCR using langevin
% framework to account for noise in reaction process

clear
close all

addpath('../utilities')

% make paths 
FigPath = '../../fig/detection_limits/';
mkdir(FigPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define static reaction paramters
Kd = 3e4; % fudging this high to account for competition from other reactions. It does not alter the results, just timescales
S0 = 50; 
k_cat = 200; % catalytic rate for activated Cas13
n_mol = 6.022e23; % mole 
reaction_vol = 1e-5; % reaction volume


%% (1) generate plots illustrating critical window concept
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define variables
C0 = .1; % concentration of Cas13 in nM
A0 = 1e-8; % concentration of activator moleules in nM
b_cat = 1e-8; % ratio of background Cas13 activity to activated Cas13;
% calculate key derivative quantities for simulation
NC = C0 / 1e9 * n_mol * reaction_vol; % number of Cas13 molecules
NA = A0 / 1e9 * n_mol * reaction_vol; % number of target activator molecules
k_cat_eff = k_cat * S0 / (Kd + S0); % effective growth rate
NA_eff = NA + b_cat*NC; % render C0 in terms of effective numbers of activator

% calculate time step and reaction time
t_sim = log(NC/NA/100)/k_cat_eff; % time to simulate
dt =  1 / k_cat_eff / 1e3;
% n_steps = t_sim / dt;
n_sim = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate forward in time
t_vec = 0:dt:t_sim;
a_array = NaN(numel(t_vec),n_sim);
a_array(1,:) = NA_eff;%normrnd(NA_eff,sqrt(NA_eff),1,n_sim);

% pre-draw gaussian noise terms
gauss_array = normrnd(0,1,numel(t_vec)-1,n_sim);


for t = 2:numel(t_vec)                
    a_curr_vec = a_array(t-1,:);        
    % update
    a_new = a_curr_vec + a_curr_vec*dt*k_cat_eff + gauss_array(t-1,:).*sqrt(a_curr_vec*dt*k_cat_eff);
    a_new(a_new<0) = 0;
    a_array(t,:) = a_new;    
end

t_ref = 1;
[~, sec_ind] = min(abs(t_vec-t_ref));
% render a_array in nM
a_array_nM = a_array / n_mol * 1e9 /reaction_vol;

% find max and min values and indices
[P01, p1_ind] = max(a_array(sec_ind,:));
[P02 ,p2_ind]= min(a_array(sec_ind,:));





%% make critical window figure
c_window_fig = figure;
cmap = brewermap(9,'Set2');
colormap(cmap);
hold on
p1 = plot(t_vec,a_array_nM,'LineWidth',1,'Color',[cmap(8,:) .4]);
p2 = plot(t_vec,nanmean(a_array_nM,2),'Color','black','LineWidth',1.5);
plot(t_vec,a_array_nM(:,p1_ind),'LineWidth',1.5,'Color',cmap(2,:));
plot(t_vec,a_array_nM(:,p2_ind),'LineWidth',1.5,'Color',cmap(3,:));
xlim([0 10])
grid on
xlabel('time (seconds)')
ylabel('free activator (nM)')
set(gca,'YScale','log')
legend([p1(1) p2],'single experiments','mean response','Location','northwest')
set(gca,'Fontsize',12)
saveas(c_window_fig,[FigPath 'critical_window_pt1.png'])

%% make normalized critical window figure
c_window_fig = figure;
cmap = brewermap(9,'Set2');
colormap(cmap);
hold on
plot(t_vec,a_array_nM./nanmean(a_array_nM,2),'LineWidth',1,'Color',[cmap(8,:) .4]);
% p2 = plot(t_vec,ones(size(t_vec)),'Color','black','LineWidth',1.5);
p1 = plot(t_vec,a_array_nM(:,p1_ind)./nanmean(a_array_nM,2),'LineWidth',1.5,'Color',cmap(2,:));
p2 = plot(t_vec,a_array_nM(:,p2_ind)./nanmean(a_array_nM,2),'LineWidth',1.5,'Color',cmap(3,:));
xlim([0 25])
grid on
xlabel('time (seconds)')
ylabel('free activator (de-trended)')
legend([p1 p2],'replicate 1','replicate 2','Location','northwest')
set(gca,'Fontsize',12)
saveas(c_window_fig,[FigPath 'critical_window_pt1_norm.png'])


% make sigmoid plot
NS = S0 / 1e9 * n_mol * reaction_vol;
mean_p = mean(a_array(sec_ind,:));
t_max = 3/k_cat_eff * log((NS-mean_p)/mean_p);
t_long = linspace(0,t_max,1e4);

%% first solve for starting concentration of each
sig_trend1 = NS ./ (1 + (NS - P01)/P01 * exp(-k_cat_eff*t_long)) / n_mol * 1e9 /reaction_vol;
sig_trend2 = NS ./ (1 + (NS - P02)/P02 * exp(-k_cat_eff*t_long)) / n_mol * 1e9 /reaction_vol;

sigmoid_plt = figure;
hold on
plot(t_long,sig_trend1,'Color',cmap(2,:),'LineWidth',1.5);
plot(t_long,sig_trend2,'Color',cmap(3,:),'LineWidth',1.5);
xlim([0 80])
xlabel('time (seconds)')
ylabel('free activator (nM)')
set(gca,'Fontsize',12)

h = breakxaxis([20 45]);
% legend(h,'experimental replicate 1','experimental replicate 2','Location','northwest')

saveas(sigmoid_plt,[FigPath 'critical_window_pt2.png'])


% pd_cov = sqrt(2/NA_eff);   
% figure;    
% plot(t_vec,nanstd(a_array,[],2)./nanmean(a_array,2));    
% hold on
% plot(t_vec,repelem(pd_cov,numel(t_vec)));

%% (2) illustrate basic testing scheme
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define variables
C0 = 0.1; % concentration of Cas13 in nM
A0 = 1e-8; % concentration of activator moleules in nM
b_cat = 10^-6.5; % ratio of background Cas13 activity to activated Cas13;
% calculate key derivative quantities for simulation
NC = C0 / 1e9 * n_mol * reaction_vol; % number of Cas13 molecules
NA = A0 / 1e9 * n_mol * reaction_vol; % number of target activator molecules
k_cat_eff = k_cat * S0 / (Kd + S0); % effective growth rate
NA_eff_pos = NA + b_cat*NC; % render C0 in terms of effective numbers of activator
NA_eff_neg = b_cat*NC;
% calculate time step and reaction time
t_sim = log(NC/NA/100)/k_cat_eff; % time to simulate
dt =  1 / k_cat_eff / 1e3;
n_steps = t_sim / dt;
n_sim = 2e3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate forward in time
t_vec = 0:dt:t_sim;
a_array_pos = NaN(numel(t_vec),n_sim);
a_array_neg = NaN(numel(t_vec),n_sim);
a_array_pos(1,:) = NA_eff_pos;%normrnd(NA_eff_pos,sqrt(NA_eff_pos),1,n_sim);
a_array_neg(1,:) = NA_eff_neg;%normrnd(NA_eff_neg,sqrt(NA_eff_pos),1,n_sim);

% pre-draw gaussian noise terms
% gauss_array_pos = normrnd(0,1,numel(t_vec)-1,n_sim);
% gauss_array_neg = normrnd(0,1,numel(t_vec)-1,n_sim);

%
tic
for t = 2:numel(t_vec)                
    a_curr_pos = a_array_pos(t-1,:);        
    a_curr_neg = a_array_neg(t-1,:); 
    % draw noise terms
    gauss_vec_pos = normrnd(0,1,1,n_sim);
    gauss_vec_neg = normrnd(0,1,1,n_sim);

    % update
    a_new_pos = a_curr_pos + a_curr_pos*dt*k_cat_eff + gauss_vec_pos.*sqrt(a_curr_pos*dt*k_cat_eff);
    a_new_pos(a_new_pos<0) = 0;
    a_array_pos(t,:) = a_new_pos; 
    
    a_new_neg = a_curr_neg + a_curr_neg*dt*k_cat_eff + gauss_vec_neg.*sqrt(a_curr_neg*dt*k_cat_eff);
    a_new_neg(a_new_neg<0) = 0;
    a_array_neg(t,:) = a_new_neg; 
end
toc

%
t_ref_2 = 10;
[~, sec_ind] = min(abs(t_vec-t_ref_2));
% render a_array in nM
a_array_pos_nM = a_array_pos / n_mol * 1e9 /reaction_vol;
a_array_neg_nM = a_array_neg / n_mol * 1e9 /reaction_vol;


% make pos vs neg plots
close all
plot_ids = randsample(1:n_sim,200,false);

pos_vs_neg_fig = figure;
cmap = brewermap(9,'Set2');
colormap(cmap);
hold on
plot(t_vec,a_array_pos_nM(:,plot_ids),'LineWidth',.5,'Color',[cmap(2,:) .2]);
plot(t_vec,a_array_neg_nM(:,plot_ids),'LineWidth',.5,'Color',[cmap(3,:) .2]);
p1 = plot(t_vec,nanmean(a_array_pos_nM,2),'LineWidth',1.5,'Color',cmap(2,:));
p2 = plot(t_vec,nanmean(a_array_neg_nM,2),'LineWidth',1.5,'Color',cmap(3,:));

xlim([0 t_ref_2*1.1])
grid on
xlabel('time (seconds)')
ylabel('free activator (nM)')
% set(gca,'YScale','log')
legend([p1 p2],'replicates (positive sample)','replicates (negative sample)','Location','northwest')
set(gca,'Fontsize',12)
saveas(pos_vs_neg_fig,[FigPath 'pos_vs_neg_pt1.png'])

% plot distributions over time 
neg_vec = a_array_neg_nM(sec_ind,:);
pos_vec = a_array_pos_nM(sec_ind,:);
n_bins = 50;
pos_bins = linspace(prctile(pos_vec,.2),prctile(pos_vec,99.8),n_bins);
neg_bins = pos_bins + nanmean(neg_vec) - nanmean(pos_vec);
pos_vs_neg_hist = figure;
hold on
h2 = histogram(neg_vec,neg_bins,'FaceColor',cmap(3,:),'EdgeAlpha',.1,'Normalization','probability');
h1 = histogram(pos_vec,pos_bins,'FaceColor',cmap(2,:),'EdgeAlpha',.1,'Normalization','probability');
xlabel('free activator (nM)')
ylabel('share')
set(gca,'Fontsize',12)
legend([h1 h2], 'replicates (positive sample)','replicates (negative sample)','Location','northeast')
saveas(pos_vs_neg_hist,[FigPath 'pos_vs_neg_hist.png'])

%% (3) plot limit of detection as a function of C and b_c

A0_fun = @(bc,NC) 4*(bc*NC)^.5;
% specify input vectors
NC_vec = [.1 1 10 100] * n_mol / 1e9 * reaction_vol;
bc_vec = logspace(-10,0);
% initialize output array
A_min_array = NaN(numel(bc_vec), numel(NC_vec));
for i = 1:numel(NC_vec)
    for j = 1:numel(bc_vec)
        A_min_array(j,i) = A0_fun(bc_vec(j),NC_vec(i));
    end
end

A_min_array_fM =  A_min_array / n_mol * 1e18 /reaction_vol;


detection_lim = figure;
hold on
cmap = brewermap(numel(NC_vec)+1,'Reds');
p = [];
for i = 1:numel(NC_vec)
    plot(bc_vec,A_min_array_fM(:,i),'Color',cmap(i+1,:),'LineWidth',1.5);
end
set(gca,'XScale','log')
set(gca,'YScale','log')
grid on
xlabel('relative Cas13 background activity (b_c)')
ylabel('minimum detectable RNA level (aM)')
legend(p,'0.1 nM RNP','1 nM RNP','10 nM RNP','100 nM RNP','Location','northwest')
set(gca,'Fontsize',12)
saveas(detection_lim,[FigPath 'detection_limit_plot.png'])