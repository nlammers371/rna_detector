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

%% basic parameters
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
rate_vec_ideal = [1 1 1 200 0];
rate_vec_actual = [1 1 1 200 200*2e-4];
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

% Can we estimate contribution from Cas13 bound to virus
close all
RC = 60*rate_vec_actual(4)*R0*(y_pos_actual(:,6) ./ (y_pos_actual(:,2) + y_pos_actual(:,6)));
CC = 60*rate_vec_actual(5)*y_pos_actual(:,1);
AC = 60*rate_vec_actual(4)*(y_pos_actual(:,6) - R0*(y_pos_actual(:,6) ./ (y_pos_actual(:,2) + y_pos_actual(:,6))));
% AC(t_pos_actual<=8.8e-6) = 0;
cat_cap_fig = figure;
hold on
a3 = area(t_pos_actual/60,AC,'FaceColor',cmap1(5,:),'FaceAlpha',0.8);
a1 = area(t_pos_actual/60,CC,'FaceColor',cmap1(3,:),'FaceAlpha',0.8);
a2 = area(t_pos_actual/60,RC,'FaceColor',cmap1(2,:),'FaceAlpha',0.8);

legend([a1 a2 a3],'Cas13 (free)','Cas13 (bound to target)','Cas13 (uncaged)','Location','northwest')
xlabel('time (minutes)')
ylabel('catalytic potential (nM^{-1}min^{-1})')
xlim([0 t_max/60])
ylim([0 1.1*max(AC)])
grid on
set(gca,'Fontsize',14)
saveas(cat_cap_fig,[FigPath 'catalytic_capacities.png'])
set(gca,'XScale','log','YScale','log')
saveas(cat_cap_fig,[FigPath 'catalytic_capacities_log.png'])