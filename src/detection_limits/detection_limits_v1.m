clear
close all

% make figure path
FigPath = '../../fig/detection_limits/';
mkdir(FigPath);

% (1) First simulate drawing from gaussian for positive and negative 
n_mol = 6.022e23;
reaction_vol = 1e-6;
% set simulation parameters
kc = 200; % catalytic rate of activated Cas13
C0 = 50 * n_mol /1e9 *reaction_vol; % concentration of RNP (Cas13 + Guide)
A0 = 1 * n_mol /1e9 *reaction_vol /1e6; % concentration of target
bc = 10^-5; % 1/SNR fro Cas13

t_vec = logspace(-6,0,1e2); % times to simulate

% define simple functions for pos and neg control
pos_mean_fun = @(t) t*(C0-A0)*bc*kc + t*A0*kc;
pos_sigma_fun = @(t) sqrt(t*(C0-A0)*bc*kc + t*A0*kc);

neg_mean_fun = @(t) t*C0*bc*kc;
neg_sigma_fun = @(t) sqrt(t*C0*bc*kc);

% plot simulated distributions over time
cmap = brewermap(9,'Set2');
n_replicates = 1e6;
% hist_bins1 = linspace(1-5e-6,1+5e-6,1e3);


time = 1e-3;
% calculate parameters
p_mu = pos_mean_fun(time);
p_sigma = pos_sigma_fun(time) / p_mu;
n_mu = neg_mean_fun(time) / p_mu;
n_sigma = neg_sigma_fun(time) / p_mu;
p_mu = p_mu / p_mu;
% "draw" samples
pos_samples = normrnd(p_mu,p_sigma,1,n_replicates);
neg_samples = normrnd(n_mu,n_sigma,1,n_replicates);

% make plots
ps_fig = figure;
hold on
histogram(pos_samples,1e3,'Normalization','probability','EdgeAlpha',0,'FaceColor',cmap(9,:));
histogram(neg_samples,1e3,'Normalization','probability','EdgeAlpha',0,'FaceColor',cmap(8,:));

xlabel('catalytic events (normalized)')
ylabel('share')
legend('positive sample (1 fM)','negetive sample','Location','southwest') 
set(gca,'Fontsize',12)

saveas(ps_fig,[FigPath 'pos_vs_neg_dist_plot.png'])

diff_fig = figure;

histogram(pos_samples-neg_samples,1e3,'Normalization','probability','EdgeAlpha',0,'FaceColor',cmap(3,:));

xlabel('difference btw positive and negative samples')
ylabel('share')
set(gca,'Fontsize',12)

saveas(diff_fig,[FigPath 'delta_dist_plot.png'])

%%
time_series_folder = [FigPath 'time_series_plots/'];
mkdir(time_series_folder);
fn_rate_vec = NaN(1,numel(t_vec));
for t = 25:80%numel(t_vec)
    time = t_vec(t);
    % calculate parameters
    p_mu = pos_mean_fun(time);
    p_sigma = pos_sigma_fun(time) / p_mu;
    n_mu = neg_mean_fun(time) / p_mu;
    n_sigma = neg_sigma_fun(time) / p_mu;
    p_mu = p_mu / p_mu;
    % "draw" samples
    pos_samples = normrnd(p_mu,p_sigma,1,n_replicates);
    neg_samples = normrnd(n_mu,n_sigma,1,n_replicates);
    
    % make plots
    delta_dist = pos_samples-neg_samples;
    mx = prctile(abs(delta_dist),99.9);
    h_bins = linspace(-mx,mx,5e2);
    delta_fig = figure;%('Visible','off');
    hold on
    h = histogram(pos_samples-neg_samples,h_bins,'Normalization','probability','EdgeAlpha',0,'FaceColor',cmap(3,:));    
    histogram([delta_dist(delta_dist<=0) repelem(1e7,sum(delta_dist>0))],h_bins,'Normalization','probability','EdgeAlpha',0,'FaceColor',cmap(2,:));    
    xlabel('catalytic events (normalized)')
    ylabel('share')
    legend('all replicates','false negatives','Location','southwest')
    set(gca,'Fontsize',12)    
    ylim([0 10e-3])
    y_lim = get(gca,'ylim');
    x_lim = get(gca,'xlim');
    text(0.95*x_lim(1),0.95*y_lim(2),['t = ' num2str(t-24)])
    
    saveas(delta_fig,[time_series_folder 'delta_dist_plot_t' sprintf('%03d',t) '.png'])
    
%     if t == 1
%         
%         fn_fig = figure;
%         hold on
%         
%         h = histogram(delta_dist,1e3,'Normalization','probability','EdgeAlpha',0,'FaceColor',cmap(3,:));    
%         histogram([delta_dist(delta_dist<=0) repelem(1e7,sum(delta_dist>0))],h.BinEdges,'Normalization','probability','EdgeAlpha',0,'FaceColor',cmap(2,:));    
%         xlabel('catalytic events (normalized)')
%         ylabel('share')
%         set(gca,'Fontsize',14)
%         legend('all replicates','false negetives')
%         saveas(fn_fig,[FigPath 'delta_dist_plot_fn.png'])
%     end   
    fn_rate_vec(t) = mean(delta_dist<=0);
    close all
end

%% plot false negative rate over time
for t = 1:numel(t_vec)
    time = t_vec(t);
    % calculate parameters
    p_mu = pos_mean_fun(time);
    p_sigma = pos_sigma_fun(time) / p_mu;
    n_mu = neg_mean_fun(time) / p_mu;
    n_sigma = neg_sigma_fun(time) / p_mu;
    p_mu = p_mu / p_mu;
    % "draw" samples
    pos_samples = normrnd(p_mu,p_sigma,1,n_replicates);
    neg_samples = normrnd(n_mu,n_sigma,1,n_replicates);
    
    % make plots
    delta_dist = pos_samples-neg_samples;
    fn_rate_vec(t) = mean(delta_dist<=0);
end
%%
fn_rate = figure;
hold on
plot(t_vec,fn_rate_vec,'Color',cmap(2,:),'LineWidth',1.5);
plot(t_vec,repelem(0.025,numel(t_vec)),'--','Color','black')
xlabel('time (seconds)')
ylabel('false negative rate')
set(gca,'Fontsize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
grid on
saveas(fn_rate,[FigPath 'fn_rate_plot.png'])

%% Plot minimum detectable A0 as a function of time

% function giving minimum detectable A given other parameters
a_min_fun = @(t)2.*((-1)+bc).^(-2).*kc.^(-1).*t.^(-1).*(1+(-1).*bc+(((-1)+bc).^2.*... 
  (1+2.*bc.*C0.*kc.*t)).^(1/2));


a_min_vec = a_min_fun(t_vec) / n_mol / reaction_vol / 1e-15;

min_a_fig = figure;
hold on
plot(t_vec,a_min_vec,'Color',cmap(2,:),'LineWidth',1.5);
xlabel('time (seconds)')
ylabel('idealized detection limit (fM)')
set(gca,'Fontsize',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
grid on
saveas(min_a_fig,[FigPath 'min_a_plot.png'])


%% now incorporate self-limitng nature of NCR
a_min_fun_ncr = @(C0,bc) (2/15)+(1/15).*(4+(-90).*((-1)+bc).^(-1).*bc.*C0).*(8+30.*(((-1)+ ...
          bc).^(-4).*((-3).*((-1)+bc).^2.*bc.^2.*C0.^2+2025.*bc.^4.*C0.^4)) ...
          .^(1/2)+270.*((-1)+bc).^(-2).*bc.*C0.*(1+bc.*((-1)+5.*C0))).^( ...
          -1/3)+(1/15).*(8+30.*(((-1)+bc).^(-4).*((-3).*((-1)+bc).^2.* ...
          bc.^2.*C0.^2+2025.*bc.^4.*C0.^4)).^(1/2)+270.*((-1)+bc).^(-2).* ...
          bc.*C0.*(1+bc.*((-1)+5.*C0))).^(1/3);