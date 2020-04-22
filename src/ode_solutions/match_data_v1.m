clear
close all
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

t_max = 1e4;
f_ind = 5;

% Generate expected signal curves using full ODE model
S0 = 200;
AI0 = 0;
C0 = 250;
R0_high = 1e-2;
R0_low = 1e-6;


kon = 0.25;
koff_s = 1;
kcat_high = 200; % this is much lower
%%%% expectation under 1/5000 SNR assumption
kcat_low_1 = kcat_high*2e-4;
koff_ns_1 = 1e4; % NOTE: bumped this up by an order of magnitude
km1 = (kcat_high + koff_ns_1)/kon;
%%%% adjusted
kcat_low_2 = kcat_high*2e-7;
koff_ns_2 = 1.5e3; % NOTE: bumped this up by an order of magnitude

%%%% set initial conditions

y0_pos_vec_high = zeros(1,size(Q_mat,1));
y0_pos_vec_low = zeros(1,size(Q_mat,1));
y0_pos_vec_high(1:4) = [C0 R0_high AI0 S0];
y0_pos_vec_low(1:4) = [C0 R0_low AI0 S0];

y0_neg_vec = y0_pos_vec_high;
y0_neg_vec(2) = 0;

% rate vector
rate_vec_1 = [kon koff_ns_1 koff_s kcat_high kcat_low_1];
rate_vec_2 = [kon koff_ns_2 koff_s kcat_high kcat_low_2];

% generate rate vec
p_vec_1 = NaN(size(rate_vec_first_order));
p_vec_2 = NaN(size(rate_vec_first_order));
% generate valued long rate vec
for r = 1:numel(rate_vec_1)
    rr = rate_vec(r);
    ri = find(rr==rate_vec_first_order);    
    p_vec_1(ri) = double(subs(rr,rr,rate_vec_1(r))); 
    p_vec_2(ri) = double(subs(rr,rr,rate_vec_2(r))); 
end

% solve ODEs numerically
[t_pos_high_1,y_pos_high_1] = ode15s(@(t,y) ncr_solver(t,y,p_vec_1,Q_mat),[0 t_max],y0_pos_vec_high);
[t_pos_low_1,y_pos_low_1] = ode15s(@(t,y) ncr_solver(t,y,p_vec_1,Q_mat),[0 t_max],y0_pos_vec_low);

[t_pos_high_2,y_pos_high_2] = ode15s(@(t,y) ncr_solver(t,y,p_vec_2,Q_mat),[0 t_max],y0_pos_vec_high);
[t_pos_low_2,y_pos_low_2] = ode15s(@(t,y) ncr_solver(t,y,p_vec_2,Q_mat),[0 t_max],y0_pos_vec_low);

[t_neg_1,y_neg_1] = ode15s(@(t,y) ncr_solver(t,y,p_vec_1,Q_mat),[0 t_max],y0_neg_vec);
[t_neg_2,y_neg_2] = ode15s(@(t,y) ncr_solver(t,y,p_vec_2,Q_mat),[0 t_max],y0_neg_vec);



% add factor to mimic experimental results
fluo_factor = 12000;
close all

% under original background assumption
match_fig_1 = figure;
cmap2 = brewermap([],'Set2');
hold on
p1 = plot(t_pos_high_1,y_pos_high_1(:,f_ind)/S0*fluo_factor,'-','Color',cmap2(3,:),'LineWidth',1.5);
p3 = plot(t_neg_1,y_neg_1(:,f_ind)/S0*fluo_factor,'-','Color',cmap2(7,:),'LineWidth',1.5);
p2 = plot(t_pos_low_1,y_pos_low_1(:,f_ind)/S0*fluo_factor,'--','Color',cmap2(3,:),'LineWidth',1.5);

ylim([0 fluo_factor*1.05])
xlim([0 2000])

xlabel('time (seconds)')
ylabel('fluorescent signal (au)')
legend([p1 p2 p3],'0.01 nM activator','0.001 nM activator','Cas/G only (5e3 SNR)','Location','northwest')

set(gca,'FontSize',14)
grid on
saveas(match_fig_1,[FigPath 'cas13_snr_orig.png'])

% under updated background assumption
match_fig_2 = figure;
hold on
p1 = plot(t_pos_high_2,y_pos_high_2(:,f_ind)/S0*fluo_factor,'-','Color',cmap2(3,:),'LineWidth',1.5);
p3 = plot(t_neg_2,y_neg_2(:,f_ind)/S0*fluo_factor,'-','Color',cmap2(7,:),'LineWidth',1.5);
p2 = plot(t_pos_low_2,y_pos_low_2(:,f_ind)/S0*fluo_factor,'--','Color',cmap2(3,:),'LineWidth',1.5);

ylim([0 fluo_factor*1.05])
xlim([0 2000])

xlabel('time (seconds)')
ylabel('fluorescent signal (au)')
legend([p1 p2 p3],'0.01 nM activator','0.001 nM activator','Cas/G only (5e6 SNR)','Location','northwest')

set(gca,'FontSize',14)
grid on
saveas(match_fig_2,[FigPath 'cas13_snr_updated.png'])


