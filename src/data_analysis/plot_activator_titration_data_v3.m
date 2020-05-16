clear
close all
addpath('../utilities')

% specify which project we're dealing with
DataName = '051420_NCRData_round2';
DataRoot = '../../out/data_analysis/';
load([DataRoot DataName '_processed_data.mat'])

ncr_project_indices = 2:6;
primary_project_indices = [1 7];


for p = ncr_project_indices    
  project = master_struct(p).project;
  % load data
  x_lim = [0 20];
  
  % make path for figures
  FigRoot = '../../fig/data_analysis/';
  FigPath = [FigRoot project '/'];
  mkdir(FigPath);

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
  % initialize max and time to max arrays
  max_delta_array = NaN(length(ncr_plot_struct)-1,2);
  t_max_array = NaN(length(ncr_plot_struct)-1,2);

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

    % plot individualdata sets
    full_plot = figure('Visible','off');
    hold on
    
    plot(time,ncr_mean,'Color',cmap(2,:)/1.2,'LineWidth',1.5);    
    plot(time,primary_mean,'Color',cmap(3,:)/1.2,'LineWidth',1.5);
    ncr = [];
    for n = 1:size(ncr_data,2)
      ncr = [ncr scatter(time,ncr_data(:,n),25,sym_cell{ncr_outliers(n)+1},'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','black','MarkerFaceAlpha',0.7)];
    end
    prim = [];
    for p = 1:size(primary_data,2)
      prim = [prim scatter(time,primary_data(:,p),25,sym_cell{primary_outliers(p)+1},'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','black','MarkerFaceAlpha',0.7)];
    end
    xlabel('time (minutes)')
    ylabel('fluorescence (au)')
    legend([ncr(1) prim(1)], ['ncr replicates (' num2str(act1_val) 'nM)'],['primary only replicates (' num2str(act1_val) 'nM)']','Location','southeast')
    set(gca,'Fontsize',12)
    xlim(x_lim)
    % save
    saveas(full_plot,[FigPath 'raw_data_scatter_' num2str(1e6*act1_val) 'fM.png'])

%     % mean and se plots
%     mean_plot = figure;
%     hold on  
% 
    ncr_ub = ncr_mean + ncr_std;
    ncr_lb = ncr_mean - ncr_std;
    prim_ub = primary_mean + primary_std;
    prim_lb = primary_mean - primary_std;
% 
%     fill([time fliplr(time)], [ncr_lb fliplr(ncr_ub)],cmap(2,:),'EdgeColor','black','FaceAlpha',0.5,'EdgeAlpha',0)
%     fill([time fliplr(time)], [prim_lb fliplr(prim_ub)],cmap(3,:),'EdgeColor','black','FaceAlpha',0.5,'EdgeAlpha',0)
% 
%     s1 = plot(time,ncr_mean,'Color',cmap(2,:)/1.2,'LineWidth',1.5);    
%     s2 = plot(time,primary_mean,'Color',cmap(3,:)/1.2,'LineWidth',1.5);    
% 
% 
%     xlabel('time (minutes)')
%     ylabel('fluorescence (au)')
%     legend([s1 s2], ['ncr (' num2str(act1_val) 'nM)'],['primary (' num2str(act1_val) 'nM)']','Location','southeast')
%     set(gca,'Fontsize',12)
%     xlim(x_lim)
%     % save
%     saveas(mean_plot,[FigPath 'mean_data_plot_' num2str(1e6*act1_val) 'fM.png'])    

    % comparison plots    

    % basic ref plot
    mean_ref_plot = figure('Visible','off');
    hold on  

    fill([time fliplr(time)], [ncr_lb fliplr(ncr_ub)],cmap(2,:),'EdgeColor','black','FaceAlpha',0.5,'EdgeAlpha',0)
    fill([time fliplr(time)], [prim_lb fliplr(prim_ub)],cmap(3,:),'EdgeColor','black','FaceAlpha',0.5,'EdgeAlpha',0)
    fill([time fliplr(time)], [cage_lb fliplr(cage_ub)],cmap(8,:),'EdgeColor','black','FaceAlpha',0.5,'EdgeAlpha',0)
    fill([time fliplr(time)], [negative_lb fliplr(negative_ub)],cmap(9,:),'EdgeColor','black','FaceAlpha',0.5,'EdgeAlpha',0)

    s1 = plot(time,ncr_mean,'Color',cmap(2,:)/1.2,'LineWidth',1.5);    
    s2 = plot(time,primary_mean,'Color',cmap(3,:)/1.2,'LineWidth',1.5); 
    s3 = plot(time,cage_mean,'Color',cmap(8,:)/1.2,'LineWidth',1.5);    
    s4 = plot(time,negative_mean,'Color',cmap(9,:)/1.2,'LineWidth',1.5); 


    xlabel('time (minutes)')
    ylabel('fluorescence (au)')
    legend([s1 s2 s3 s4], ['ncr (' num2str(act1_val) 'nM)'],['primary (' num2str(act1_val) 'nM)']','cage only','no activator','Location','southeast')
    set(gca,'Fontsize',12)
    xlim(x_lim)
    % save
    saveas(mean_ref_plot,[FigPath 'mean_ref_data_plot_' num2str(1e6*act1_val) 'fM.png']) 

    % difference plot
    diff_plot = figure('Visible','off');
    hold on  
    f_size = 1;
%     if a == 1
%       f_size = 1;
%     end
    ncr_delta_sm = imgaussfilt(ncr_mean-cage_mean,f_size);
    primary_delta_sm = imgaussfilt(primary_mean-negative_mean,f_size);

    s1 = plot(time,ncr_delta_sm,'Color',cmap(2,:),'LineWidth',2);    
    s2 = plot(time,primary_delta_sm,'Color',cmap(3,:),'LineWidth',2);  

    scatter(time,ncr_mean-cage_mean,25,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','black','MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',1)
    scatter(time,primary_mean-negative_mean,25,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','black','MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',1)        

    xlabel('time (minutes)')
    ylabel('relative fluorescence (au)')
    legend([s1 s2], ['ncr delta (' num2str(act1_val) 'nM)'],['primary delta (' num2str(act1_val) 'nM)']','Location','southeast')
    set(gca,'Fontsize',12)
    grid on
    xlim(x_lim)
    % save
    saveas(diff_plot,[FigPath 'diff_data_plot_' num2str(1e6*act1_val) 'fM_.png']) 

    % record metrics
    [max_delta_array(a,1), mi] = max(primary_delta_sm);
    t_max_array(a,1) = time(mi);

    [max_delta_array(a,2), mi] = max(ncr_delta_sm);
    t_max_array(a,2) = time(mi);
    
  end


  cmap1 = brewermap([],'set2');
  cmap2 = flipud(brewermap(9,'Blues'));
  cmap3 = flipud(brewermap(9,'Reds'));
  close all

  
  plot_indices = 1:size(max_delta_array,1);


  % check for 
  max_delta_array_plot = max_delta_array;
  max_delta_array_plot(max_delta_array_plot<=0) = realmin;

  shift_fig = figure;
  hold on
  plot(max_delta_array_plot(plot_indices,:)'./t_max_array(plot_indices,:)',max_delta_array_plot(plot_indices,:)','Color','black','LineWidth',1)
  s1 = [];
  s2 = [];
  for i = plot_indices
    s1 = [s1 scatter(max_delta_array_plot(i,1)./t_max_array(i,1),max_delta_array_plot(i,1),300/i,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','black','LineWidth',1)];
    s2 = [s2 scatter(max_delta_array_plot(i,2)./t_max_array(i,2),max_delta_array_plot(i,2),300/i,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','black','LineWidth',1)];
  end
  grid on
  xlabel('detection rate')
  ylabel('maximum difference')
  set(gca,'XScale','log')
  set(gca,'YScale','log')
  legend([s1(1) s2(1)],'primary only','NCR','Location','southeast')
  axis([5e-2 5e2 1e0 1e4])
  set(gca,'Fontsize',12)
  saveas(shift_fig,[FigPath 'ncr_shift_plot.png'])
  
  ncr_fig = figure;
  hold on
  for i = 1:size(ncr_titration_array,2)
    plot(time,ncr_titration_array(:,i),'Color',cmap3(i,:),'LineWidth',1.5)
  end
  
  grid on
  xlabel('time (minutes)')
  ylabel('fluorescence (au)')
  set(gca,'XScale','log')
  set(gca,'YScale','log')
%   legend([s1(1) s2(1)],'primary only','NCR','Location','southeast')
  xlim(x_lim)
  set(gca,'Fontsize',12)
  saveas(ncr_fig,[FigPath 'ncr_titration_plot.png'])
  
  primary_fig = figure;
  primary_titration_array(primary_titration_array<0) = 0;
  hold on
  for i = 1:size(ncr_titration_array,2)
    plot(time,primary_titration_array(:,i),'Color',cmap2(i,:),'LineWidth',1.5)
  end
  
  grid on
  xlabel('time (minutes)')
  ylabel('fluorescence (au)')
  set(gca,'XScale','log')
  set(gca,'YScale','log')
%   legend([s1(1) s2(1)],'primary only','NCR','Location','southeast')
  xlim(x_lim)
  set(gca,'Fontsize',12)
  saveas(primary_fig,[FigPath 'primary_titration_plot.png'])
   
end