clear
close all
addpath('../utilities')

% specify which project we're dealing with
project_cell = {'050720_NCR018'};
activator_rows = 3:10;
for p = 1:length(project_cell)
  project = project_cell{p};
  % load data
  load(['../../out/data_analysis/' project '.mat']);
  
  % make path for figures
  FigRoot = '../../fig/data_analysis/';
  FigPath = [FigRoot project '/'];
  mkdir(FigPath);

%   activator_rows = plot_struct(1).activator_rows;
  % extract reference data
  cage_data = plot_struct(end).ncr_data;
  cage_outliers = plot_struct(end).ncr_outlier_flags;
  negative_data = plot_struct(end).primary_data;
  negative_outliers = plot_struct(end).primary_outlier_flags;

  negative_mean = plot_struct(end).primary_mean';
  negative_std = plot_struct(end).primary_std';

  cage_mean = plot_struct(end).ncr_mean';
  cage_std = plot_struct(end).ncr_std';

  cage_ub = cage_mean + cage_std;
  cage_lb = cage_mean - cage_std;
  negative_ub = negative_mean + negative_std;
  negative_lb = negative_mean - negative_std;

  % generate color map
  cmap = brewermap([],'Set2');
  % plot symbol shapes
  sym_cell = {'o','s'};
  % initialize max and time to max arrays
  max_delta_array = NaN(numel(activator_rows),2);
  t_max_array = NaN(numel(activator_rows),2);

  for a = [1 3 4 6]%1:numel(activator_rows)
    close all
    % get activator val
    ac_val = plot_struct(a).ac_val;
    % extract data
    ncr_data = plot_struct(a).ncr_data;
    ncr_outliers = plot_struct(a).ncr_outlier_flags;
    primary_data = plot_struct(a).primary_data;
    primary_outliers = plot_struct(a).primary_outlier_flags;

    primary_mean = plot_struct(a).primary_mean';
    primary_std = plot_struct(a).primary_std';

    ncr_mean = plot_struct(a).ncr_mean';
    ncr_std = plot_struct(a).ncr_std';

    time = plot_struct(a).time';

    % plot individualdata sets
    full_plot = figure;
    hold on
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
    legend([ncr(1) prim(1)], ['ncr replicates (' num2str(ac_val) 'nM)'],['primary only replicates (' num2str(ac_val) 'nM)']','Location','southeast')
    set(gca,'Fontsize',12)
    xlim([0 max(time)])
    % save
    saveas(full_plot,[FigPath num2str(1e3*ac_val) 'pM_raw_data_scatter.png'])

    % mean and se plots
    mean_plot = figure;
    hold on  

    ncr_ub = ncr_mean + ncr_std;
    ncr_lb = ncr_mean - ncr_std;
    prim_ub = primary_mean + primary_std;
    prim_lb = primary_mean - primary_std;

    fill([time fliplr(time)], [ncr_lb fliplr(ncr_ub)],cmap(2,:),'EdgeColor','black','FaceAlpha',0.5,'EdgeAlpha',0)
    fill([time fliplr(time)], [prim_lb fliplr(prim_ub)],cmap(3,:),'EdgeColor','black','FaceAlpha',0.5,'EdgeAlpha',0)

    s1 = plot(time,ncr_mean,'Color',cmap(2,:)/1.2,'LineWidth',1.5);    
    s2 = plot(time,primary_mean,'Color',cmap(3,:)/1.2,'LineWidth',1.5);    


    xlabel('time (minutes)')
    ylabel('fluorescence (au)')
    legend([s1 s2], ['ncr (' num2str(ac_val) 'nM)'],['primary (' num2str(ac_val) 'nM)']','Location','southeast')
    set(gca,'Fontsize',12)
    xlim([0 max(time)])
    % save
    saveas(mean_plot,[FigPath num2str(1e3*ac_val) 'pM_mean_data_plot.png'])    

    % comparison plots
    if a <= numel(activator_rows)

      % basic ref plot
      mean_ref_plot = figure;
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
      legend([s1 s2 s3 s4], ['ncr (' num2str(ac_val) 'nM)'],['primary (' num2str(ac_val) 'nM)']','cage only','no activator','Location','southeast')
      set(gca,'Fontsize',12)
      xlim([0 max(time)])
      % save
      saveas(mean_ref_plot,[FigPath num2str(1e3*ac_val) 'pM_mean_ref_data_plot.png']) 

      % difference plot
      diff_plot = figure;
      hold on  
      f_size = 10;
      if a == 1
        f_size = 1;
      end
      ncr_delta_sm = imgaussfilt(ncr_mean-cage_mean,f_size);
      primary_delta_sm = imgaussfilt(primary_mean-negative_mean,f_size);

      s1 = plot(time,ncr_delta_sm,'Color',cmap(2,:),'LineWidth',2);    
      s2 = plot(time,primary_delta_sm,'Color',cmap(3,:),'LineWidth',2);  
      
      scatter(time,ncr_mean-cage_mean,25,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','black','MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',1)
      scatter(time,primary_mean-negative_mean,25,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','black','MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',1)        

      xlabel('time (minutes)')
      ylabel('relative fluorescence (au)')
      legend([s1 s2], ['ncr delta (' num2str(ac_val) 'nM)'],['primary delta (' num2str(ac_val) 'nM)']','Location','southeast')
      set(gca,'Fontsize',12)
      grid on
      xlim([0 max(time)])
      % save
      saveas(diff_plot,[FigPath num2str(1e3*ac_val) 'pM_diff_data_plot.png']) 

      % record metrics
      [max_delta_array(a,1), mi] = max(primary_delta_sm);
      t_max_array(a,1) = time(mi);

      [max_delta_array(a,2), mi] = max(ncr_delta_sm);
      t_max_array(a,2) = time(mi);
    end
  end


  cmap1 = brewermap([],'set2');
  cmap2 = brewermap(9,'Blues');
  close all

  if strcmp(project,'NCR018')
    plot_indices = [1 2 3 4 6];
  else  
    plot_indices = 1:size(max_delta_array,1);
  end

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
end