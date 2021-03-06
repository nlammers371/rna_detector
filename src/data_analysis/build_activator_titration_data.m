% script to read in tecan data, group columns, and plot
clear
close all
addpath('../utilities')
% set path to data
DataFolder = 'C:\Users\nlamm\Dropbox (Personal)\rna_detection\np_data\';

% DataName = '050720_NCR018.xlsx';
% SheetName = 'NCR018';
DataName = '051220_NCRData';

DataPath = [DataFolder DataName '.xlsx'];

% make path for cleaned figure data
WriteRoot = '../../out/data_analysis/';

% get sheet info
[~,SheetNames]=xlsfinfo(DataPath);
SheetIndex = find(ismember(SheetNames,'DataSheet'));
MapIndices = find(contains(SheetNames, 'map'));
n_maps = numel(MapIndices);

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

% Iterate through map tabs and generate data sets
% specify whether neg control is in separate tab
if strcmp(DataName,'051220_NCRData')
  neg_sep_flag = true;
  neg_map_ind = 5;
  neg_group_id = 38;
else
  neg_sep_flag = false;
  neg_map_ind = NaN;
end

activator_rows = 4:8;
NCR_rows = 9:10;
ExpTabIndices = MapIndices(~ismember(MapIndices,neg_map_ind));
group_offset = 0;

for e = ExpTabIndices

  sheet_name = SheetNames{e};  
  underscore_ind = strfind(sheet_name,'_');
  project = [DataName '_' sheet_name(underscore_ind+1:end)];
  
  %%%%%%%%%%%%%%
  % Read in plate map
  %%%%%%%%%%%%%%
  WritePath = [WriteRoot '/'];
  mkdir(WritePath);

  % read sheet
  [~,~,MapCell] = xlsread(DataPath,e);

  %%%%%%%%%%%%%%
  % Compile data to plot
  %%%%%%%%%%%%%%

  % initialize structure to store plotting info
  plot_struct = struct;
  full_table = table;
  stats_table = table;

  % check which groups are NCR
  ncr_flag_vec = max(isnan(cell2mat(MapCell(NCR_rows,3:end))))==0;
  
  for a = 1:numel(activator_rows)+1
      % get basic info
      if a <= numel(activator_rows)      
        ac_ind = activator_rows(a);
        ac_val = MapCell{ac_ind,2};
        plot_struct(a).act_ncr_groups = find(ncr_flag_vec&~isnan(cell2mat(MapCell(ac_ind,3:end))))+group_offset;
        plot_struct(a).act_prim_groups = find(~ncr_flag_vec&~isnan(cell2mat(MapCell(ac_ind,3:end))))+group_offset;
      else      
        ac_val = -1;
        ac_ind = NaN;
        plot_struct(a).act_ncr_groups = find(ncr_flag_vec&min(isnan(cell2mat(MapCell(activator_rows,3:end))))==0)+group_offset;
        if neg_sep_flag
          plot_struct(a).act_prim_groups = neg_group_id;
        else
          plot_struct(a).act_prim_groups = find(~ncr_flag_vec&min(isnan(cell2mat(MapCell(activator_rows,3:end))))==0)+group_offset;
        end
      end
      % record
      plot_struct(a).activator_rows = activator_rows;
      plot_struct(a).ac_ind = ac_ind;
      plot_struct(a).ac_val = ac_val;
      % pull raw values to plotd
      plot_struct(a).ncr_cols_orig = find(ismember(group_id_vec,plot_struct(a).act_ncr_groups));
      plot_struct(a).ncr_data = cell2mat(ExpCell(2:end,plot_struct(a).ncr_cols_orig));
      plot_struct(a).ncr_data = plot_struct(a).ncr_data - plot_struct(a).ncr_data(1,:);
      plot_struct(a).primary_cols_orig = find(ismember(group_id_vec,plot_struct(a).act_prim_groups));
      plot_struct(a).primary_data = cell2mat(ExpCell(2:end,plot_struct(a).primary_cols_orig));
      plot_struct(a).primary_data = plot_struct(a).primary_data - plot_struct(a).primary_data(1,:);
      % flag potential outlier series
      plot_struct(a).ncr_outlier_flags = prctile(isoutlier(plot_struct(a).ncr_data,2),33)==1;
      plot_struct(a).primary_outlier_flags = prctile(isoutlier(plot_struct(a).primary_data,2),33)==1;
      % calculate summary statistics
      plot_struct(a).ncr_mean = nanmean(plot_struct(a).ncr_data(:,~plot_struct(a).ncr_outlier_flags),2);
      plot_struct(a).ncr_std = nanstd(plot_struct(a).ncr_data(:,~plot_struct(a).ncr_outlier_flags),[],2);
      plot_struct(a).primary_mean = nanmean(plot_struct(a).primary_data(:,~plot_struct(a).primary_outlier_flags),2);
      plot_struct(a).primary_std = nanstd(plot_struct(a).primary_data(:,~plot_struct(a).primary_outlier_flags),[],2);
      % time vector
      plot_struct(a).time = cell2mat(ExpCell(2:end,2))/60;
  end
  % save
  save([WritePath project '.mat'],'plot_struct')
  % increment
  group_offset = group_offset + size(MapCell(activator_rows,3:end),2);
end
