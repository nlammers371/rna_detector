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
R0 = 1e-2;

bkg_cat_ratio_vec = logspace(-8,-4,10);%1e-8;
% specify ratesq
kon = 1;
koff_s = 1;
koff_ns = 1e3;
kcat_high = 5e1;
kcat_low = kcat_high*1e-8;

% set initial conditions positive and negative samples

% NCR samples (low R0)
y0_pos_vec = zeros(1,size(Q_mat,1));
y0_pos_vec(1:4) = [C0 R0 AI0 S0];
y0_neg_vec = y0_pos_vec;
y0_neg_vec(2) = 0;

% rate vector
rate_vec_actual = [kon koff_ns koff_s kcat_high kcat_low];


sim_struct = struct;

for k = 1:numel(bkg_cat_ratio_vec)
    rate_vec_iter = rate_vec_actual;
    rate_vec_iter(end) = rate_vec_iter(end-1)*bkg_cat_ratio_vec(k);
    % generate rate vec
    p_vec_iter = NaN(size(rate_vec_first_order));

    % generate valued long rate vec
    for r = 1:numel(rate_vec_actual)
        rr = rate_vec(r);
        ri = find(rr==rate_vec_first_order);    
        p_vec_iter(ri) = double(subs(rr,rr,rate_vec_iter(r))); 
    end
    % solve ODEs numerically    
    [t_neg,y_neg] = ode15s(@(t,y) ncr_solver(t,y,p_vec_iter,Q_mat),[0 t_max],y0_neg_vec);
    sim_struct(k).t_neg = t_neg;
    sim_struct(k).y_neg = y_neg;
end

% generate rate vec
p_vec_pos = NaN(size(rate_vec_first_order));

% generate valued long rate vec
for r = 1:numel(rate_vec_actual)
    rr = rate_vec(r);
    ri = find(rr==rate_vec_first_order);    
    p_vec_pos(ri) = double(subs(rr,rr,rate_vec_actual(r))); 
end
[t_pos,y_pos] = ode15s(@(t,y) ncr_solver(t,y,p_vec_pos,Q_mat),[0 t_max],y0_pos_vec);

%%
fluo_factor = 12000;
close all


bkg_fig = figure;
cmap1 = brewermap(10,'Blues');
cmap2 = brewermap([],'Set2');
hold on
for i = 1:numel(sim_struct)
    plot(sim_struct(i).t_neg,sim_struct(i).y_neg(:,f_ind)/S0*fluo_factor,'-','Color',cmap1(i,:),'LineWidth',1.5);
end
plot(t_pos,y_pos(:,f_ind)/S0*fluo_factor,'-','Color',cmap2(2,:),'LineWidth',1.5);


ylim([0 fluo_factor*1.05])
xlim([0 2000])

xlabel('time (seconds)')
ylabel('fluorescent signal (percent of max)')
% legend([p1 p2 p3],'0.01 nM activator','0.0001 nM activator','negative control','Location','northwest')

set(gca,'FontSize',14)
% set(gca,'XScale','log')
grid on
saveas(bkg_fig,[FigPath 'cas13_snr_titration.png'])

