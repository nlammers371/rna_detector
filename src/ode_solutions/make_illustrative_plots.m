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

% basic parameters
t_max = 1e4;
f_ind = find(strcmp(atom_string_cell,'F'));
% specify initial concentrations
S0 = 200;
AI0 = 20;
C0 = 20;
R0 = 1e-4;
R0_high = 1e-1;

% specify rates
kon = 1;
koff_s = 1;
koff_ns = 2e3;
kcat_high = 2e2;
kcat_low = 4e-2;

% set initial conditions positive and negative samples

% NCR samples (low R0)
y0_pos_vec = zeros(1,size(Q_mat,1));
y0_pos_vec(1:4) = [C0 R0 AI0 S0];
y0_neg_vec = y0_pos_vec;
y0_neg_vec(2) = 0;

% sherlock (no NCR)
y0_pos_vec_sh = y0_pos_vec;
y0_pos_vec_sh(3) = 0;

% samples (with high R0)
y0_pos_vec_high = y0_pos_vec;
y0_pos_vec_high(2) = R0_high;
y0_pos_vec_high_sh = y0_pos_vec_high;
y0_pos_vec_high_sh(3) = 0;
y0_neg_vec_high_sh = y0_pos_vec_high_sh;
y0_neg_vec_high_sh(2) = 0;

% idealized and non idealized rate vectors
rate_vec_ideal = [kon koff_ns koff_s kcat_high 0];
rate_vec_actual = [kon koff_ns koff_s kcat_high kcat_low];
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

[t_pos_high,y_pos_high] = ode15s(@(t,y) ncr_solver(t,y,p_vec_actual,Q_mat),[0 t_max],y0_pos_vec_high);
[t_pos_high_sh,y_pos_high_sh] = ode15s(@(t,y) ncr_solver(t,y,p_vec_actual,Q_mat),[0 t_max],y0_pos_vec_high_sh);
[t_neg_high_sh,y_neg_high_sh] = ode15s(@(t,y) ncr_solver(t,y,p_vec_actual,Q_mat),[0 t_max],y0_neg_vec_high_sh);

%% (1) plots for idealized case
detection_threshold = repelem(50,numel(t_pos_ideal));
fluo_pos_ideal_sh = y_pos_ideal_sh(:,f_ind)/S0*100;
fluo_pos_ideal = y_pos_ideal(:,f_ind)/S0*100;

% determine x limits
t1_ideal = round(10*2*t_pos_ideal(find(fluo_pos_ideal>1e-4,1)))/10;
t2_ideal = round(10*t_pos_ideal(find(fluo_pos_ideal>100-1e-3,1)))/10;

% "Sherlock" only plot
close all
ideal_fig_sh = figure;
cmap1 = brewermap(9,'Set2');
hold on
plot(t_pos_ideal_sh,fluo_pos_ideal_sh,'Color',cmap1(1,:),'LineWidth',1.5);
plot(t_pos_ideal,detection_threshold,'--','Color','black')
ylim([0 110])
xlim([t1_ideal t2_ideal])
xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')
set(gca,'FontSize',14)
% set(gca,'XScale','log')
grid on
saveas(ideal_fig_sh,[FigPath 'illustrative_sher_ideal.png'])

% add in NCR (idealized)
ideal_fig_pos = figure;
hold on
p1 = plot(t_pos_ideal,y_pos_ideal(:,f_ind)/S0*100,'-','Color',cmap1(2,:),'LineWidth',1.5);
p2 = plot(t_pos_ideal_sh,y_pos_ideal_sh(:,f_ind)/S0*100,'-','Color',cmap1(1,:),'LineWidth',1.5);
plot(t_pos_ideal,detection_threshold,'--','Color','black')
ylim([0 110])
xlim([t1_ideal t2_ideal])
xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')
legend([p1 p2],'with NCR','without NCR','Location','southwest')
set(gca,'FontSize',14)
% set(gca,'XScale','log')
grid on
saveas(ideal_fig_pos,[FigPath 'illustrative_ideal.png'])

ideal_fig_pos_log = figure;
hold on
p1 = plot(t_pos_ideal,y_pos_ideal(:,f_ind)/S0*100,'-','Color',cmap1(2,:),'LineWidth',1.5);
p2 = plot(t_pos_ideal_sh,y_pos_ideal_sh(:,f_ind)/S0*100,'-','Color',cmap1(1,:),'LineWidth',1.5);
plot(t_pos_ideal,detection_threshold,'--','Color','black')

ylim([1e-5 110])
xlim([t1_ideal t2_ideal])

xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')
legend([p1 p2],'with NCR','without NCR','Location','southwest')

set(gca,'FontSize',14)
% set(gca,'XScale','log')
set(gca,'YScale','log')

grid on
saveas(ideal_fig_pos_log,[FigPath 'illustrative_ideal_log.png'])

%% (2) plots for situation when background Cas13 activity is included
close all

fluo_pos_actual = y_pos_actual(:,f_ind)/S0*100;
% determine x limits
t1_act = 1e-4*t_pos_actual(find(fluo_pos_ideal>1e-4,1));
t2_act = 1.1*t_pos_actual(find(fluo_pos_ideal>100-1,1));

% positive control only
actual_fig_pos = figure;
hold on
p1 = plot(t_pos_actual,y_pos_actual(:,f_ind)/S0*100,'-','Color',cmap1(2,:),'LineWidth',1.5);
plot(t_pos_ideal,detection_threshold,'--','Color','black')

ylim([0 110])
xlim([t1_act t2_act])

xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')

% set(gca,'XScale','log')
set(gca,'FontSize',14)
grid on

saveas(actual_fig_pos,[FigPath 'illustrative_actual.png'])


% Add in negative control
actual_fig_pos_neg = figure;
hold on
p1 = plot(t_pos_actual,y_pos_actual(:,f_ind)/S0*100,'-','Color',cmap1(2,:),'LineWidth',1.5);
p2 = plot(t_neg_actual,y_neg_actual(:,f_ind)/S0*100,'--','Color',cmap1(3,:),'LineWidth',1.5);
plot(t_pos_ideal,detection_threshold,'--','Color','black')

ylim([0 110])
xlim([t1_act t2_act])

xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')
legend([p1 p2],'positive sample','negative sample','Location','southeast')

set(gca,'FontSize',14)
% set(gca,'XScale','log')
grid on
saveas(actual_fig_pos_neg,[FigPath 'illustrative_actual_pos_neg.png'])

actual_fig_pos_neg_log = figure;
hold on
p1 = plot(t_pos_actual,y_pos_actual(:,f_ind)/S0*100,'-','Color',cmap1(2,:),'LineWidth',1.5);
p2 = plot(t_neg_actual,y_neg_actual(:,f_ind)/S0*100,'--','Color',cmap1(3,:),'LineWidth',1.5);
plot(t_pos_ideal,detection_threshold,'--','Color','black')

ylim([1e-4 110])
xlim([t1_act t2_act])

xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')
legend([p1 p2],'positive sample','negative sample','Location','southeast')

set(gca,'FontSize',14)
set(gca,'YScale','log')

grid on
saveas(actual_fig_pos_neg_log,[FigPath 'illustrative_actual_pos_neg_log.png'])

%% Show that R0 high enough to discern using NCR is high enough to detect using sherlock
close all

% Add in negative control
high_r_fig = figure;
hold on
p1 = plot(t_pos_high,y_pos_high(:,f_ind)/S0*100,'-','Color',cmap1(2,:),'LineWidth',1.5);
p2 = plot(t_neg_actual,y_neg_actual(:,f_ind)/S0*100,'--','Color',cmap1(3,:),'LineWidth',1.5);
plot(t_pos_ideal,detection_threshold,'--','Color','black')

ylim([0 110])
xlim([t1_act t2_act])

xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')
legend([p1 p2],'positive sample','negative sample','Location','southeast')

set(gca,'FontSize',14)
% set(gca,'XScale','log')
grid on
saveas(high_r_fig,[FigPath 'illustrative_high_r0.png'])

% add sherlock
high_r_sh_fig = figure;
hold on
% NCR
p1 = plot(t_pos_high,y_pos_high(:,f_ind)/S0*100,'-','Color',cmap1(2,:),'LineWidth',1.5);
p2 = plot(t_neg_actual,y_neg_actual(:,f_ind)/S0*100,'--','Color',cmap1(3,:),'LineWidth',1.5);
% Sherlock
p3 = plot(t_pos_high_sh,y_pos_high_sh(:,f_ind)/S0*100,'-*','Color',cmap1(2,:),'LineWidth',1.5);
p4 = plot(t_neg_high_sh,y_neg_high_sh(:,f_ind)/S0*100,'-o','Color',cmap1(3,:),'LineWidth',1.5);

plot(t_pos_ideal,detection_threshold,'--','Color','black')

ylim([0 110])
xlim([t1_act 5e3])

xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')
legend([p1 p2 p3 p4],'positive sample (NCR)','negative sample (NCR)',...
    'positive sample (sherlock)','negative sample (sherlock)','Location','southeast')

set(gca,'FontSize',14)
% set(gca,'XScale','log')
grid on
saveas(high_r_sh_fig,[FigPath 'illustrative_high_r0_sherlock.png'])


high_r_sh_fig = figure;
hold on
% NCR
p1 = plot(t_pos_high,y_pos_high(:,f_ind)/S0*100,'-','Color',cmap1(2,:),'LineWidth',1.5);
p2 = plot(t_neg_actual,y_neg_actual(:,f_ind)/S0*100,'--','Color',cmap1(3,:),'LineWidth',1.5);
% Sherlock
p3 = plot(t_pos_high_sh,y_pos_high_sh(:,f_ind)/S0*100,'-*','Color',cmap1(2,:),'LineWidth',1.5);
p4 = plot(t_neg_high_sh,y_neg_high_sh(:,f_ind)/S0*100,'-o','Color',cmap1(3,:),'LineWidth',1.5);

plot(t_pos_ideal,detection_threshold,'--','Color','black')

ylim([0 110])
xlim([1 t_max])

xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')
legend([p1 p2 p3 p4],'positive sample (NCR)','negative sample (NCR)',...
    'positive sample (sherlock)','negative sample (sherlock)','Location','southeast')

set(gca,'FontSize',14)
set(gca,'XScale','log')
grid on
saveas(high_r_sh_fig,[FigPath 'illustrative_high_r0_sherlock_log.png'])

%%

high_low_r_fig = figure;
hold on
p1 = plot(t_pos_high,y_pos_high(:,f_ind)/S0*100,'-','Color',cmap1(2,:),'LineWidth',1.5);
p2 = plot(t_pos_actual,y_pos_actual(:,f_ind)/S0*100,'-','Color',cmap1(5,:),'LineWidth',1.5);
p3 = plot(t_neg_actual,y_neg_actual(:,f_ind)/S0*100,'--','Color',cmap1(3,:),'LineWidth',1.5);
plot(t_pos_ideal,detection_threshold,'--','Color','black')

ylim([0 110])
xlim([t1_act t2_act])

xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')
legend([p1 p2 p3],'high RNA','low RNA','negative control','Location','southeast')

set(gca,'FontSize',14)
% set(gca,'XScale','log')
grid on
saveas(high_low_r_fig,[FigPath 'illustrative_high_low_r0.png'])

%% Estimate contribution from Cas13 bound to virus
close all
% low condition
f_bound_low = y_pos_actual(:,6) ./ (y_pos_actual(:,2) + y_pos_actual(:,6));
RC_low = rate_vec_actual(4)*R0*f_bound_low;
CC_low = rate_vec_actual(5)*y_pos_actual(:,1);
AC_low = rate_vec_actual(4)*(y_pos_actual(:,6)-R0*f_bound_low);
AC_low(AC_low<0) = 0;
% high condition
f_bound_high = y_pos_high(:,6) ./ (y_pos_high(:,2) + y_pos_high(:,6));
RC_high = rate_vec_actual(4)*R0_high*f_bound_high;
CC_high = rate_vec_actual(5)*y_pos_high(:,1);
AC_high = rate_vec_actual(4)*(y_pos_high(:,6)-R0*f_bound_high);
AC_high(AC_high<0) = 0;
% AC(t_pos_actual<=8.8e-6) = 0;


cat_cap_fig = figure;
hold on
a3 = area(t_pos_actual,AC_low,'FaceColor',cmap1(5,:),'FaceAlpha',0.8);
a1 = area(t_pos_actual,CC_low,'FaceColor',cmap1(3,:),'FaceAlpha',0.8);
a2 = area(t_pos_actual,RC_low,'FaceColor',cmap1(2,:),'FaceAlpha',0.8);

legend([a1 a2 a3],'Cas13 (free)','Cas13 (bound to target)','Cas13 (uncaged)','Location','northwest')
xlabel('time (seconds)')
ylabel('catalytic potential (nM^{-1}min^{-1})')
xlim([0 t_max])
ylim([0 1.1*max(AC_low)])
grid on
set(gca,'Fontsize',14)
saveas(cat_cap_fig,[FigPath 'catalytic_capacities.png'])
set(gca,'XScale','log','YScale','log')
saveas(cat_cap_fig,[FigPath 'catalytic_capacities_log.png'])


rc_frac_vec = RC_low ./ (RC_low + CC_low + AC_low) * 100;

rna_cat = figure;
hold on
p1 = area(t_pos_actual,rc_frac_vec,'FaceColor',cmap1(2,:),'EdgeAlpha',0,'FaceAlpha',0.8);

xlim([1e-5 t_max])
xlabel('time (seconds)')
ylabel('catalytic activity (percent of total)')
% legend([p1 p2],'positive sample','negative sample','Location','southeast')

set(gca,'FontSize',14)
set(gca,'XScale','log')
grid on

yyaxis right
a2 = area(t_pos_actual,y_pos_actual(:,f_ind)/S0*100,'FaceColor',cmap1(5,:),'FaceAlpha',0.5,'EdgeAlpha',0);
ylabel('fluorescent signal (percent of max)')
ylim([0 105])
ax = gca;
ax.YColor = 'black';

legend('Cas13 (bound to target)','reporter fluorescence','Location','southwest')
saveas(rna_cat,[FigPath 'target_cat_contribution.png'])

%%%%
rna_cat = figure;
hold on
p1 = area(t_pos_actual,rc_frac_vec,'FaceColor',cmap1(2,:),'EdgeAlpha',0,'FaceAlpha',0.8);
ylabel('catalytic activity (percent of total)')
xlim([1e-5 t2_act])
xlabel('time (seconds)')

% legend([p1 p2],'positive sample','negative sample','Location','southeast')

set(gca,'FontSize',14)
% set(gca,'XScale','log')
grid on

yyaxis right
a2 = area(t_pos_actual,y_pos_actual(:,f_ind)/S0*100,'FaceColor',cmap1(5,:),'FaceAlpha',0.5,'EdgeAlpha',0);
ylabel('fluorescent signal (percent of max)')
ylim([0 105])
ax = gca;
ax.YColor = 'black';

legend('Cas13 (bound to target)','reporter fluorescence','Location','southeast')
saveas(rna_cat,[FigPath 'target_cat_contribution_lin.png'])


%% Make stacked area plot
close  all

y_raw_low = [RC_low AC_low CC_low];
y_norm_low = y_raw_low ./ sum(y_raw_low,2);


% High condition
area_abs = figure;
h = area(t_pos_actual,y_raw_low,'EdgeAlpha',0,'FaceAlpha',0.8);
h(1).FaceColor = cmap1(2,:);
h(2).FaceColor = cmap1(5,:);
h(3).FaceColor = cmap1(3,:);
xlabel('time (seconds)')
ylabel('catalytic activity (nM sec^{-1})')
set(gca,'FontSize',12)
grid on
xlim([0 t2_act])
saveas(area_abs,[FigPath 'cat_area_plot_abs.png'])


area_norm = figure;
h = area(t_pos_actual,y_norm_low,'EdgeAlpha',0,'FaceAlpha',0.6);
h(1).FaceColor = cmap1(2,:);
h(2).FaceColor = cmap1(5,:);
h(3).FaceColor = cmap1(3,:);

xlim([0 t2_act])
ylim([0 1])

legend('Cas13 (bound to target)','Cas13 (bound to activator)','Cas13 (free)','Location','northwest')
xlabel('time (seconds)')
ylabel('catalytic activity (share of total)')
set(gca,'FontSize',12)
set(gca,'XScale','log')
grid on

saveas(area_norm,[FigPath 'cat_area_plot_norm.png'])


%%%
y_raw_high = [RC_high AC_high CC_high];
y_norm_high = y_raw_high ./ sum(y_raw_high,2);


% High condition
area_abs = figure;
h = area(t_pos_high,y_raw_high,'EdgeAlpha',0,'FaceAlpha',0.8);
h(1).FaceColor = cmap1(2,:);
h(2).FaceColor = cmap1(5,:);
h(3).FaceColor = cmap1(3,:);
xlabel('time (seconds)')
ylabel('catalytic activity (nM sec^{-1})')
set(gca,'FontSize',12)
grid on
xlim([0 t2_act])
% legend('Cas13 (bound to target)','Cas13 (bound to activator)','Cas13 (free)')
saveas(area_abs,[FigPath 'cat_area_plot_abs_high.png'])


area_norm = figure;
h = area(t_pos_high,y_norm_high,'EdgeAlpha',0,'FaceAlpha',0.6);
h(1).FaceColor = cmap1(2,:);
h(2).FaceColor = cmap1(5,:);
h(3).FaceColor = cmap1(3,:);

xlim([0 t2_act])
ylim([0 1])

xlabel('time (seconds)')
ylabel('catalytic activity (share of total)')
set(gca,'FontSize',12)
set(gca,'XScale','log')
grid on

legend('Cas13 (bound to target)','Cas13 (bound to activator)','Cas13 (free)','Location','northwest')
saveas(area_norm,[FigPath 'cat_area_plot_norm_high.png'])


%% overlay reporter trend with total cat activity
last_fig = figure;

hold on
h = area(t_pos_high,sum(y_raw_high,2),'EdgeAlpha',0,'FaceAlpha',0.5);
ylabel('catalytic activity (nM sec^{-1})')

yyaxis right
p1 = plot(t_pos_actual,y_pos_actual(:,f_ind)/S0*100,'-','Color','black','LineWidth',1.5);
ylabel('fluorescent signal (percent of max)')
ylim([0 105])
legend('catalytic activity','reporter signal','Location','southeast')

xlabel('time (seconds)')
xlim([0 t2_act])
set(gca,'FontSize',12)
ax = gca;
ax.YColor = 'black';
grid on

saveas(last_fig,[FigPath 'cat_fluo_trend.png'])


%%

free_a = AI0 - y_pos_actual(:,3);
r = nanmedian(f_bound_low);
start_ind = find(t_pos_actual>1,1);
P0 = R0+C0/5000;%AI0-y_pos_actual(start_ind,3);

t_vec = t_pos_actual;%(start_ind:end);
p_test = AI0 ./ (1 + (AI0-P0)/P0 *exp(-r*t_vec));