function view_probe_trajectory(~,~,histology_toolbar_gui)
% Part of AP_histology toolbox
%
% View histology slices with overlaid aligned CCF areas

% Initialize guidata
gui_data = struct;

% Store toolbar handle
gui_data.histology_toolbar_gui = histology_toolbar_gui;

% Load atlas
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
if isempty(allen_atlas_path)
    error('No CCF atlas found (add CCF atlas to path)')
end
disp('Loading Allen CCF atlas...')
gui_data.tv = readNPY(fullfile(allen_atlas_path,'template_volume_10um.npy'));
gui_data.av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
gui_data.st = ap_histology.loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));
disp('Done.')

% Get images (from path in toolbar GUI)
histology_toolbar_guidata = guidata(histology_toolbar_gui);
gui_data.save_path = histology_toolbar_guidata.save_path;

slice_dir = dir(fullfile(gui_data.save_path,'*.tif'));
slice_fn = natsortfiles(cellfun(@(path,fn) fullfile(path,fn), ...
    {slice_dir.folder},{slice_dir.name},'uni',false));

gui_data.slice_im = cell(length(slice_fn),1);
for curr_slice = 1:length(slice_fn)
    gui_data.slice_im{curr_slice} = imread(slice_fn{curr_slice});
end

% Load corresponding CCF slices
ccf_slice_fn = fullfile(gui_data.save_path,'histology_ccf.mat');
load(ccf_slice_fn);
gui_data.histology_ccf = histology_ccf;

% Load histology/CCF alignment
ccf_alignment_fn = fullfile(gui_data.save_path,'atlas2histology_tform.mat');
load(ccf_alignment_fn);
gui_data.histology_ccf_alignment = atlas2histology_tform;

% Initialize probe points
lines_colormap = lines(7);
probe_colormap = [lines_colormap;lines_colormap(:,[2,3,1]);lines_colormap(:,[3,1,2])];

gui_data.probe_color = probe_colormap;
gui_data.probe_points_histology = cell(length(gui_data.slice_im),1);
gui_data.probe_lines = gobjects(1);

% Plot probe trajectories
save_fn = fullfile(gui_data.save_path,'probe_ccf.mat');
load(save_fn,'probe_ccf');
plot_probe(gui_data,probe_ccf);

end

%%%
% Copied from annotate_neuropixels.m
%%%
function plot_probe(gui_data,probe_ccf)

% Plot probe trajectories
figure('Name','Probe trajectories');
axes_atlas = axes;
[~, brain_outline] = plotBrainGrid([],axes_atlas);
set(axes_atlas,'ZDir','reverse');
hold(axes_atlas,'on');
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[ap_max,dv_max,ml_max] = size(gui_data.tv);
xlim([-10,ap_max+10])
ylim([-10,ml_max+10])
zlim([-10,dv_max+10])
h = rotate3d(gca);
h.Enable = 'on';

for curr_probe = 1:length(probe_ccf)
    % Plot points and line of best fit
    r0 = mean(probe_ccf(curr_probe).points,1);
    xyz = bsxfun(@minus,probe_ccf(curr_probe).points,r0);
    [~,~,V] = svd(xyz,0);
    histology_probe_direction = V(:,1);
    % (make sure the direction goes down in DV - flip if it's going up)
    if histology_probe_direction(2) < 0
        histology_probe_direction = -histology_probe_direction;
    end

    line_eval = [-1000,1000];
    probe_fit_line = bsxfun(@plus,bsxfun(@times,line_eval',histology_probe_direction'),r0);
    plot3(probe_ccf(curr_probe).points(:,1), ...
        probe_ccf(curr_probe).points(:,3), ...
        probe_ccf(curr_probe).points(:,2), ...
        '.','color',gui_data.probe_color(curr_probe,:),'MarkerSize',20);
    line(probe_fit_line(:,1),probe_fit_line(:,3),probe_fit_line(:,2), ...
        'color',gui_data.probe_color(curr_probe,:),'linewidth',2)
end

% Plot probe areas
figure('Name','Trajectory areas');
for curr_probe = 1:length(probe_ccf)

    curr_axes = subplot(1,length(probe_ccf),curr_probe);

    trajectory_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
        probe_ccf(curr_probe).trajectory_areas.color_hex_triplet,'uni',false)),[1,3,2]);

    trajectory_areas_boundaries = probe_ccf(curr_probe).trajectory_areas.trajectory_depth;
    trajectory_areas_centers = mean(trajectory_areas_boundaries,2);

    trajectory_areas_image_depth = 0:0.01:max(trajectory_areas_boundaries,[],'all');
    trajectory_areas_image_idx = interp1(trajectory_areas_boundaries(:,1), ...
        1:height(probe_ccf(curr_probe).trajectory_areas),trajectory_areas_image_depth, ...
        'previous','extrap');
    trajectory_areas_image = trajectory_areas_rgb(trajectory_areas_image_idx,:,:);

    image([],trajectory_areas_image_depth,trajectory_areas_image);
    yline(unique(trajectory_areas_boundaries(:)),'color','k','linewidth',1);
    set(curr_axes,'XTick',[],'YTick',trajectory_areas_centers, ...
        'YTickLabels',probe_ccf(curr_probe).trajectory_areas.acronym);
    set(curr_axes,'XTick',[]);
    title(['Probe ' num2str(curr_probe)]);

end

end