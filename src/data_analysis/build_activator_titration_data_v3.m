% script to read in tecan data, group columns, and plot
clear
close all
addpath('../utilities')
% set path to data
DataFolder = 'C:\Users\nlamm\Dropbox (Personal)\rna_detection\np_data\';

% DataName = '050720_NCR018.xlsx';
% SheetName = 'NCR018';
DataName = '051420_NCRData_round2';

DataPath = [DataFolder DataName '.xlsx'];

% make path for cleaned figure data
WriteRoot = '../../out/data_analysis/';

% get sheet info
[~,SheetNames]=xlsfinfo(DataPath);
SheetIndex = find(ismember(SheetNames,'DataSheet'));
MapIndex = find(contains(SheetNames, 'PlateMap'));

% read sheet
[~,~,SheetCell] = xlsread(DataPath,SheetIndex);

% get index of data header row
HeaderRowIndex = find(strcmp(SheetCell(:,1),{'Cycle Nr.'}));
LastRowIndex = find(~cellfun(@isempty,SheetCell(:,8)),1,'last')-4;
LastColIndex = find(~cellfun(@isempty,SheetCell(HeaderRowIndex,:)),1,'last');

% generate trucnated sheet
ExpCell = SheetCell(HeaderRowIndex:LastRowIndex,1:LastColIndex);
HeaderRow = ExpCell(1,:);

%%%%%%%%%%%%%%
% Group columns based on position within 384 well plate
%%%%%%%%%%%%%%
col_num_vec = [];
row_num_vec = [];
row_string_cell = {};
for e = 4:length(HeaderRow)
    string = HeaderRow{e};   
    row_string_cell = [row_string_cell {string(1)}];
end

row_string_index = unique(row_string_cell);
for e = 4:length(HeaderRow)
    string = HeaderRow{e};
    row_ind = find(strcmp(string(1),row_string_index));
    col_num_vec = [col_num_vec str2double(string(2:end))];
    row_num_vec = [row_num_vec row_ind];    
end

row_vec_coarse = round(row_num_vec/2);
col_vec_coarse = round(col_num_vec/2);
group_id_vec = [NaN(1,3) (row_vec_coarse-1)*max(col_vec_coarse) + col_vec_coarse];


%%%%%%%%%%%%%%
% Extract Plate Map
%%%%%%%%%%%%%%

% read sheet
[~,~,MapCell] = xlsread(DataPath,MapIndex);

%%%%%%%%%%%%%%
% Read in plate map
%%%%%%%%%%%%%%
WritePath = [WriteRoot '/'];
mkdir(WritePath);

activator1_rows = find(contains(MapCell(:,1),'NCR_001'));
activator2_row = find(contains(MapCell(:,1),'NCR_003'));
guide1_row = find(contains(MapCell(:,1),'NCR_004'));
exp_rows = max(activator1_rows)+2:size(MapCell,1)-2;
exp_cell = MapCell(exp_rows,1);%{'100nM NCR_018','100nM NCR_266','100nM NCR_061'};

master_struct = struct;

for e = 1:length(exp_cell)

  project = [DataName '_' exp_cell{e}];    

  %%%%%%%%%%%%%%
  % Compile data to plot
  %%%%%%%%%%%%%%

  % initialize structure to store plotting info
  plot_struct = struct;  
  
  project_rows = [activator2_row find(contains(MapCell(:,1),exp_cell{e}))'];  
  
  % check which groups are NCR
  project_flag_vec = max(isnan(cell2mat(MapCell(project_rows,3:end))))==0;
  activator_col_vec = find(project_flag_vec);  
  activator_col_vec = [activator_col_vec max(activator_col_vec)+1];
  
  for n = 1:numel(activator_col_vec)
      plot_struct(n).project = project;
      % get basic info
      act_col = activator_col_vec(n);
      project_map_col = [NaN ; cell2mat(MapCell(2:end,2+act_col))];
      act1_ind = min(activator1_rows) + find(~isnan(project_map_col(activator1_rows)))-1;
      if ~isempty(act1_ind)
        act1_val = MapCell{act1_ind,2};
      else
        act1_val = 0;
      end
      act2_val = project_map_col(activator2_row)*50;
      guide1_val = project_map_col(guide1_row)*50;
      if contains(project,'primary')
          guide2_val = 0;
      else
          guide2_val = project_map_col(project_rows(2))*50;
      end
      if isnan(guide2_val)
          guide2_val = 0;
      end
      plot_struct(n).act_groups = act_col;           
      % record
      plot_struct(n).act1_val = act1_val;
      plot_struct(n).act2_val = act2_val;
      plot_struct(n).guide1_val = guide1_val;
      plot_struct(n).guide2_val = guide2_val;
      % pull raw values to plot
      plot_struct(n).cols_orig = find(ismember(group_id_vec,plot_struct(n).act_groups));
      plot_struct(n).raw_data = cell2mat(ExpCell(2:end,plot_struct(n).cols_orig));
      plot_struct(n).raw_data = plot_struct(n).raw_data - plot_struct(n).raw_data(1,:);
      
   
      % flag potential outlier series
      plot_struct(n).outlier_flags = prctile(isoutlier(plot_struct(n).raw_data,2),33)==1;
      
      % calculate summary statistics
      plot_struct(n).activator_mean = nanmean(plot_struct(n).raw_data(:,~plot_struct(n).outlier_flags),2);
      plot_struct(n).activator_std = nanstd(plot_struct(n).raw_data(:,~plot_struct(n).outlier_flags),[],2);

      % time vector
      plot_struct(n).time = cell2mat(ExpCell(2:end,2))/60;
  end
  % save
%   save([WritePath project '.mat'],'plot_struct') 
  master_struct(e).plot_struct = plot_struct;
  master_struct(e).project = project;
end
save([WritePath DataName '_processed_data.mat'],'master_struct') 