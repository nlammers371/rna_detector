clear
close all
addpath('../utilities')

% specify which project we're dealing with
DataName = '051420_NCRData_round2';
DataRoot = '../../out/data_analysis/';
load([DataRoot DataName '_processed_data.mat'])

primary_indices = [1 7];

data_array = [master_struct(1).plot_struct(1).time];
header_strings = {'time'};
activator_vec = [NaN];
for m = primary_indices
    plot_struct = master_struct(m).plot_struct;
    activator_val_vec = [plot_struct.act1_val];
    for i = 1:numel(activator_val_vec)
        % add date 
        data_array = [data_array plot_struct(i).raw_data];
        % make stringd
        ac_val = activator_val_vec(i);        
        for j = 1:size(plot_struct(i).raw_data,2)
            iter = sum(activator_vec==ac_val)+1;
            header_strings = [header_strings{:} {['activator_' num2str(ac_val) 'nM_rep' sprintf('%02d',iter)]}];
            activator_vec = [activator_vec ac_val];
        end
    end
end

primary_only_table = array2table(data_array, 'VariableNames', header_strings);
metadata_struct = struct;
metadata_struct.guide_id = 'NCR004';
metadata_struct.RNP_nM = 10;
metadata_struct.activator_id = 'NCR001';
writetable(primary_only_table,[DataRoot DataName '_primary_only_titrations.csv'])
save([DataRoot DataName '_primary_only_metadata.mat'],'metadata_struct')