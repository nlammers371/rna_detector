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
t_max = 1e4;
f_ind = find(strcmp(atom_string_cell,'F'));

% specify initial concentrations
S0 = 200;
AI0 = 20;
C0 = 20;
R0 = 1e-4;
C0_vec = logspace(-4,2);

factor =  1e4;

% specify rates
kon = 1;
koff_s = 1;
koff_ns = 2e3;
kcat_high = 2e2;
kcat_low = 4e-2;

kcat_low_vec = logspace(-8,0);

% set initial conditions positive and negative samples

% NCR samples (low R0)
y0_pos_vec = zeros(1,size(Q_mat,1));
y0_pos_vec(1:4) = [C0 R0 AI0 S0]*factor;
y0_neg_vec = y0_pos_vec;
y0_neg_vec(2) = 0;

% idealized and non idealized rate vectors
rate_vec_actual = [kon koff_ns koff_s kcat_high kcat_low];
p_vec = NaN(size(rate_vec_first_order));

% generate valued long rate vec
for r = 1:numel(rate_vec_ideal)
    rr = rate_vec(r);
    ri = find(rr==rate_vec_first_order);    
    p_vec(ri) = double(subs(rr,rr,rate_vec_actual(r))); 
end

% Cas13 SNR titration

sim_struct_snr = struct;
Opt = odeset('Events', @myEvent);
for c = 1:numel(kcat_low_vec)
    rate_vec_iter = rate_vec_actual;
    rate_vec_iter(end) = kcat_low;%kcat_low_vec(c);
    % generate rate vec
    p_vec_iter = NaN(size(rate_vec_first_order));

    y0_pos_iter = y0_pos_vec;
    y0_pos_iter(1) = C0_vec(c)*factor;%y0_pos_vec;
    y0_neg_iter = y0_neg_vec;
    y0_neg_iter(1) = C0_vec(c)*factor;%
    % generate valued long rate vec
    for r = 1:numel(rate_vec_ideal)
        rr = rate_vec(r);
        ri = find(rr==rate_vec_first_order);    
        p_vec_iter(ri) = double(subs(rr,rr,rate_vec_iter(r))); 
    end
    % solve ODEs numerically
    [t_pos,y_pos] = ode15s(@(t,y) ncr_solver(t,y,p_vec_iter,Q_mat),[0 t_max],y0_pos_iter,Opt);
    [t_neg,y_neg] = ode15s(@(t,y) ncr_solver(t,y,p_vec_iter,Q_mat),[0 t_max],y0_neg_iter,Opt);
    % record results
    sim_struct_snr(c).kcat_low = rate_vec_iter(end);
    sim_struct_snr(c).t_pos = t_pos;
    sim_struct_snr(c).t_neg = t_neg;
    sim_struct_snr(c).t2_pos = t_pos(end);
    sim_struct_snr(c).t2_neg = t_neg(end);
    sim_struct_snr(c).y_pos = y_pos;
    sim_struct_snr(c).y_neg = y_neg;    
end

%% look at results
t_obs_vec = [sim_struct_snr.t2_pos];
t_ratio_vec = [sim_struct_snr.t2_neg]./[sim_struct_snr.t2_pos];
%%
ex_plot_ind = find(t_ratio_vec >= 1.1,1,'last');

plot_vec = kcat_high./kcat_low_vec;

[~,curr_ind] = min(abs(kcat_low_vec-kcat_low));
close all


snr_figure = figure;
hold on
cmap1 = brewermap([],'Set3');

% plot detection time ratio
scatter(plot_vec,t_ratio_vec,'MarkerFaceColor',cmap1(4,:),'MarkerEdgeColor','black')
scatter(plot_vec(ex_plot_ind),t_ratio_vec(ex_plot_ind),50,'MarkerFaceColor',...
        cmap1(3,:),'MarkerEdgeColor','black')
scatter(plot_vec(curr_ind),t_ratio_vec(curr_ind),50,'MarkerFaceColor',...
    cmap1(5,:),'MarkerEdgeColor','black')

xlim([plot_vec(end) plot_vec(1)])

xlabel('SNR (kcat_{b}/kcat_{f})')
ylabel('relative detection time (t_{n}/t_{p})')
% set(gca,'Fontsize',12)
set(gca,'XScale','log')

grid on
% 
% y_limits = get(gca,'ylim');
% plot([kcat_high/kcat_low kcat_high/kcat_low],y_limits,'-','Color','black')

saveas(snr_figure,[FigPath 'Cas13_snr_titration_updated.png'])

%% make example plot
rate_vec_iter = rate_vec_actual;
rate_vec_iter(end) = kcat_low_vec(ex_plot_ind);
% generate rate vec
p_vec_iter = NaN(size(rate_vec_first_order));

% generate valued long rate vec
for r = 1:numel(rate_vec_ideal)
    rr = rate_vec(r);
    ri = find(rr==rate_vec_first_order);    
    p_vec_iter(ri) = double(subs(rr,rr,rate_vec_iter(r))); 
end
% solve ODEs numerically
[t_pos,y_pos] = ode15s(@(t,y) ncr_solver(t,y,p_vec_iter,Q_mat),[0 t_max],y0_pos_vec/factor);
[t_neg,y_neg] = ode15s(@(t,y) ncr_solver(t,y,p_vec_iter,Q_mat),[0 t_max],y0_neg_vec/factor);


close all
ex_fig = figure;
cmap2 = brewermap([],'Set2');
hold on
p1 = plot(t_pos,y_pos(:,f_ind)/S0*100,'-','Color',cmap2(2,:),'LineWidth',1.5);
p2 = plot(t_neg,y_neg(:,f_ind)/S0*100,'--','Color',cmap2(3,:),'LineWidth',1.5);
plot(t_pos,repelem(50,numel(t_pos)),'--','Color','black')

ylim([0 110])
xlim([0 20])

xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')
legend([p1 p2],'positive sample','negative control','Location','southeast')

set(gca,'FontSize',14)
% set(gca,'XScale','log')
grid on
saveas(ex_fig,[FigPath 'illustrative_improved_cas13.png'])

