% Script to simulate early exponential phase of basic NCR using langevin
% framework to account for noise in reaction process

clear
close all

addpath('../utilities')

% make paths 
FigPath = '../../fig/detection_limits_v2/';
mkdir(FigPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define static reaction paramters
Kd = 3e4; % fudging this high to account for competition from other reactions. It does not alter the results, just timescales
S0 = 50; 
k_cat = 200; % catalytic rate for activated Cas13
n_mol = 6.022e23; % mole 
reaction_vol = 1e-5; % reaction volume


% (1) generate plots illustrating critical window concept
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define variables
C0 = 1; % concentration of Cas13 in nM
A0 = 1e-4; % concentration of activator moleules in nM
b_cat = 1e-8; % ratio of background Cas13 activity to activated Cas13;
% calculate key derivative quantities for simulation
NC = C0 / 1e9 * n_mol * reaction_vol; % number of Cas13 molecules
NA = A0 / 1e9 * n_mol * reaction_vol; % number of target activator molecules
k_cat_eff = k_cat * S0 / (Kd + S0); % effective growth rate
NA_eff = NA + b_cat*NC; % render C0 in terms of effective numbers of activator

% calculate time step and reaction time
% t_sim = log(NC/NA/10)/k_cat_eff; % time to simulate
t_sim = 100;
dt =  1 / k_cat_eff / 1e3;
% n_steps = t_sim / dt;
n_sim = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate forward in time
t_vec = 0:dt:t_sim;
a_array = NaN(numel(t_vec),n_sim);
sherlock_array = NaN(numel(t_vec),n_sim);
a_array(1,:) = NA_eff;%normrnd(NA_eff,sqrt(NA_eff),1,n_sim);
sherlock_array(1,:) = NA_eff;%normrnd(NA_eff,sqrt(NA_eff),1,n_sim);

% pre-draw gaussian noise terms
gauss_array1 = normrnd(0,1,numel(t_vec)-1,n_sim);
gauss_array2 = normrnd(0,1,numel(t_vec)-1,n_sim);


for t = 2:numel(t_vec)                
    a_curr_vec = a_array(t-1,:);        
    % update
    a_new = a_curr_vec + a_curr_vec*dt*k_cat_eff + gauss_array1(t-1,:).*sqrt(a_curr_vec*dt*k_cat_eff);
    a_new(a_new<0) = 0;
    a_array(t,:) = a_new;
    
    % sherlock
    s_curr_vec = sherlock_array(t-1,:);        
    % update
    s_new = s_curr_vec + NA_eff*dt*k_cat_eff + gauss_array2(t-1,:).*sqrt(NA_eff*dt*k_cat_eff);
    s_new(s_new<0) = 0;
    sherlock_array(t,:) = s_new;       
end

% render a_array in nM
a_array_nM = a_array / n_mol * 1e9 /reaction_vol;
sherlock_array_nM = sherlock_array / n_mol * 1e9 /reaction_vol;

%% find max and min values and indices
close all

ncr_mean = nanmean(a_array_nM,2);
ncr_std = nanstd(a_array_nM,[],2);

sherlock_mean = nanmean(sherlock_array_nM,2);
sherlock_std = nanstd(sherlock_array_nM,[],2);

snr_fig = figure;
hold on
plot(t_vec(2:end),ncr_mean(2:end),'LineWidth',1.5)
plot(t_vec(2:end),sherlock_mean(2:end),'LineWidth',1.5)
xlabel('time (seconds)')
ylabel('average signal')
% set(gca,'YScale','log')
legend('NCR-like process','Sherlock-like process','Location','northwest')
set(gca,'Fontsize',12)
set(gca,'YScale','log')
xlim([1 t_sim])
grid on
saveas(snr_fig,[FigPath 'snr_plot.png'])


snr_fig = figure;
hold on
plot(t_vec(2:end),ncr_mean(2:end)./ncr_std(2:end),'LineWidth',1.5)
plot(t_vec(2:end),sherlock_mean(2:end)./sherlock_std(2:end),'LineWidth',1.5)
xlabel('time (seconds)')
ylabel('SNR (\mu / \sigma)')
% set(gca,'YScale','log')
% set(gca,'YScale','log')
legend('NCR-like process','Sherlock-like process','Location','northwest')
% set(gca,'Fontsize',12)
xlim([1 t_sim])
ylim([650 4500])
grid on
saveas(snr_fig,[FigPath 'snr_plot_detailed.png'])
