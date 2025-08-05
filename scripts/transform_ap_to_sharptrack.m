%% Load file from AP_histology
file_path = 'E:\Projects\SSA\data\brainslice\NP_probe_trajectory\52N\52N\slices\';
probe_ccf_filepath = [file_path, '\probe_ccf.mat'];
load(probe_ccf_filepath, 'probe_ccf');

%% Create depth table
% Iterate over each probe
for i = 1:size(probe_ccf, 1)
    this_probe_info = {};
    trajectory_areas = probe_ccf(i).trajectory_areas;
    for range_index = 1:size(trajectory_areas, 1)
        this_probe_info{range_index, 1} = trajectory_areas.trajectory_depth(range_index, 1);
        this_probe_info{range_index, 2} = trajectory_areas.trajectory_depth(range_index, 2);
        this_probe_info{range_index, 3} = trajectory_areas.acronym(range_index);
        this_probe_info{range_index, 4} = trajectory_areas.safe_name(range_index);
        this_probe_info{range_index, 5} = trajectory_areas.sphinx_id(range_index);        
    end
    % Sphinx_id
    ProbeInfo{i} = this_probe_info;
end
ProbeInfo{size(probe_ccf, 1) + 1} = file_path;

save([file_path, '\depth_table.mat'], "ProbeInfo")


%% Create depth table