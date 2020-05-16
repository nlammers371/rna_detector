clear
close all
addpath('../utilities')

% specify which project we're dealing with
DataName = '051420_NCRData_round2';
DataRoot = '../../out/data_analysis/';
load([DataRoot DataName '_processed_data.mat'])

ncr_project_indices = 2:6;
primary_project_indices = [1 7];

max_delta_array = NaN(7,numel(ncr_project_indices));
prim_delta_vec = NaN(7,1);
FigPath = '../../fig/data_analysis/';

for p = ncr_project_indices    
    
  project = master_struct(p).project;
  % load data
  x_lim = [0 20];
  
  % make path for figures
  

  ncr_plot_struct = master_struct(p).plot_struct;
  primary_plot_struct = master_struct(primary_project_indices(1)).plot_struct;
  guide2_vec = [ncr_plot_struct.guide2_val];
  act1_vec = [ncr_plot_struct.act1_val];
  ncr_indices = find(guide2_vec>0 &  act1_vec>0);
  cage_index = find(guide2_vec>0 &  act1_vec==0);
  
  % extract reference cage only data and negative control data
  raw_cage_data = ncr_plot_struct(cage_index).raw_data;
  cage_outliers = ncr_plot_struct(cage_index).outlier_flags;

  cage_mean = ncr_plot_struct(cage_index).activator_mean';
  cage_std = ncr_plot_struct(cage_index).activator_std';
  
  negative_data = primary_plot_struct(end).raw_data;
  negative_outliers = primary_plot_struct(end).outlier_flags;

  negative_mean = primary_plot_struct(end).activator_mean';
  negative_std = primary_plot_struct(end).activator_std';
  
  cage_ub = cage_mean + cage_std;
  cage_lb = cage_mean - cage_std;
  negative_ub = negative_mean + negative_std;
  negative_lb = negative_mean - negative_std;

  % generate color map
  cmap = brewermap([],'Set2');
  % plot symbol shapes
  sym_cell = {'o','s'};

  % ncr titration array
  ncr_titration_array = NaN(length(cage_mean),1+numel(ncr_indices));
  ncr_titration_array(:,1) = cage_mean;
  primary_titration_array = NaN(length(cage_mean),1+numel(ncr_indices));
  primary_titration_array(:,1) = negative_mean;   
  
  for a = ncr_indices
    close all
    % get activator val
    act1_val = ncr_plot_struct(a).act1_val;
    % extract data
    ncr_data = ncr_plot_struct(a).raw_data;
    ncr_outliers = ncr_plot_struct(a).outlier_flags;
    primary_data = primary_plot_struct(a).raw_data;
    primary_outliers = primary_plot_struct(a).outlier_flags;

    primary_mean = primary_plot_struct(a).activator_mean';
    primary_std = primary_plot_struct(a).activator_std';
    primary_titration_array(:,a+1) = primary_mean;
    
    ncr_mean = ncr_plot_struct(a).activator_mean';
    ncr_std = ncr_plot_struct(a).activator_std';
    ncr_titration_array(:,a+1) = ncr_mean;
    
    time = ncr_plot_struct(a).time';

    % generate smoothed diff curves
    ncr_mean_sm = imgaussfilt(ncr_mean,1);
    cage_mean_sm = imgaussfilt(cage_mean,1);
    ncr_delta = ncr_mean_sm-cage_mean_sm;
    dd = diff(diff(ncr_mean));    
    [~, mi] = min(dd);
       
    [max_delta_array(a,p-1), mi] = max(ncr_delta(1:mi+1));    
    prim_delta = imgaussfilt(primary_mean,1)-imgaussfilt(negative_mean,1);
    prim_delta_vec(a) = max(prim_delta(1:mi+1));
  end
 
 
   
end


cmap1 = brewermap([],'set2');
close all

max_fig = figure;
hold on
lgd_str = {};
p = [];
for i = 1:numel(ncr_project_indices)
    plot(1000*act1_vec(1:end-2),max_delta_array(:,i),'Color',cmap1(i,:))
    p = [p scatter(1000*act1_vec(1:end-2),max_delta_array(:,i),'MarkerFaceColor',cmap1(i,:),'MarkerEdgeAlpha',0)];
    project = master_struct(i+1).project;
    lgd_str = [lgd_str{:} {project(end-7:end)}];
end
plot(1000*act1_vec(1:end-2),prim_delta_vec,'Color','black');
p = [p scatter(1000*act1_vec(1:end-2),prim_delta_vec,'MarkerFaceColor','black','MarkerEdgeAlpha',0)];
lgd_str = [lgd_str{:} {'primary only'}];

xlabel('primary activator concentration (pM)')
ylabel('difference over background (au)')
legend(p,lgd_str{:},'Interpreter', 'none')
set(gca,'Fontsize',12)
set(gca, 'XDir','reverse')
set(gca,'XScale','log')
grid on
xlim([5e-6 5e2])
saveas(max_fig,[FigPath 'titration_summary.png'])
