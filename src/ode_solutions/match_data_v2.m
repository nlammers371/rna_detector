clear
close all
addpath('../utilities')

DataPath = '../../out/ode_studies/';
mkdir(DataPath)

% specify project to load
project = 'ncr_2activator';
% load 
load([DataPath project '_setup.mat'])

% make figure path
FigPath = ['../../fig/ode_studies/' project '/' ];
mkdir(FigPath)

t_max = 4000;
f_ind = find(strcmp(atom_string_cell,'F'));

% Generate expected signal curves using full ODE model
S0 = 200;
A2I0 = 20;
C10 = 50;
C20 = 0.8;
A10_high = 1e-2;
A10_low = 10^-5.5;
A20 = 0;


% define rates
kon = 0.25;
koff_s = 1;
kcat_high = 200; % this is much lower
kcat_low = kcat_high*2e-7;
hl_ratio = 1e-2;
koff_ns = 1.5e3; % NOTE: bumped this up by an order of magnitude
SI_ratio = 0.065;
%%%% set initial conditions

y0_pos_high = zeros(1,size(Q_mat,1));
y0_pos_high(1:6) = [C10 C20 A10_high A20 A2I0 S0];
y0_pos_low = y0_pos_high;
y0_pos_low(3) = A10_low;

y0_neg_vec = y0_pos_high;
y0_neg_vec(3) = 0;

% rate vector
rate_vec_ref = unique(rate_vec_first_order);

% generate rate vec
p_vec = NaN(size(rate_vec_first_order));
p_vec_low = NaN(size(rate_vec_first_order));
% generate valued long rate vec
for r = 1:numel(rate_vec_ref)
    rr = rate_vec_ref(r);
    ri = find(rr==rate_vec_first_order);    
    p_vec(ri) = eval(rr); 
    p_vec_low(ri) = p_vec(ri);
    if contains(char(rr),'kcat_low')
        p_vec_low(ri) = p_vec_low(ri)*hl_ratio;
    end
end

%% solve ODEs numerically
[t_pos_high,y_pos_high] = ode15s(@(t,y) ncr_solver(t,y,p_vec,Q_mat),[0 t_max],y0_pos_high);
[t_pos_low,y_pos_low] = ode15s(@(t,y) ncr_solver(t,y,p_vec,Q_mat),[0 t_max],y0_pos_low);
[t_neg,y_neg] = ode15s(@(t,y) ncr_solver(t,y,p_vec,Q_mat),[0 t_max],y0_neg_vec);
[t_neg_low,y_neg_low] = ode15s(@(t,y) ncr_solver(t,y,p_vec_low,Q_mat),[0 t_max],y0_neg_vec);

fluo_factor = 12000;
close all

% under updated background assumption
match_fig_3 = figure;
cmap2 = brewermap([],'set2');
hold on
p1 = plot(t_pos_high,y_pos_high(:,f_ind)/S0*fluo_factor,'-','Color','black','LineWidth',1.5);
p2 = plot(t_pos_low,y_pos_low(:,f_ind)/S0*fluo_factor,'--','Color','black','LineWidth',1.5);

ylim([0 fluo_factor*1.05])
xlim([0 2000])

xlabel('time (seconds)')
ylabel('fluorescent signal (au)')
legend([p1 p2],'10 pM activator (NCR)','1 fM activator (NCR)','Location','southeast')

set(gca,'FontSize',14)
grid on
saveas(match_fig_3,[FigPath 'cas13_ncr_pos_only_sl' num2str(SI_ratio) '.png'])


match_fig_4 = figure;
hold on
p1 = plot(t_pos_high,y_pos_high(:,f_ind)/S0*fluo_factor,'-','Color','black','LineWidth',1.5);
p3 = plot(t_neg,y_neg(:,f_ind)/S0*fluo_factor,'-','Color','red','LineWidth',1.5);
p2 = plot(t_pos_low,y_pos_low(:,f_ind)/S0*fluo_factor,'--','Color','black','LineWidth',1.5);

ylim([0 fluo_factor*1.05])
xlim([0 2000])

xlabel('time (seconds)')
ylabel('fluorescent signal (au)')
legend([p1 p2 p3],'10 pM activator','1 fM activator','Cas/G only (5e6 SNR)','Location','southeast')

set(gca,'FontSize',14)
grid on
saveas(match_fig_4,[FigPath 'cas13_ncr_pos_only_neg_sl' num2str(SI_ratio) '.png'])


match_fig_5 = figure;
hold on
p1 = plot(t_pos_high,y_pos_high(:,f_ind)/S0*fluo_factor,'-','Color','black','LineWidth',1.5);
% p3 = plot(t_neg,y_neg(:,f_ind)/S0*fluo_factor,'-','Color',[.5 .5 .5],'LineWidth',1.5);
p2 = plot(t_pos_low,y_pos_low(:,f_ind)/S0*fluo_factor,'--','Color','black','LineWidth',1.5);
p4 = plot(t_neg_low,y_neg_low(:,f_ind)/S0*fluo_factor,'--','Color','red','LineWidth',1.5);

ylim([0 fluo_factor*1.05])
xlim([0 t_max])

xlabel('time (seconds)')
ylabel('fluorescent signal (au)')
legend([p1 p2 p4],'10 pM activator','1 fM activator','Cas/G only (5e8 SNR)','Location','southeast')

set(gca,'FontSize',14)
grid on
saveas(match_fig_5,[FigPath 'cas13_ncr_pos_only_neg_low_sl' num2str(SI_ratio) '.png'])