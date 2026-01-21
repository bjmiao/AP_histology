replay_action_record_path = 'C:\Users\benji\Documents\harini\8506C_ket\overlay\slices';
image_path = 'C:\Users\benji\Documents\harini\8506C_ket\ppdh';
save_dir = [image_path, '\slices_replay'];
mkdir(save_dir);
%% Load the image path and the replay file

% Get and sort image files
image_path_dir = dir(fullfile(image_path,'*.tif'));
% (error if none found)
if isempty(image_path_dir)
    error('No TIFF images found in %s',image_path)
end
im_fn = natsortfiles(cellfun(@(path,fn) fullfile(path,fn), ...
    {image_path_dir.folder},{image_path_dir.name},'uni',false));

n_im = length(im_fn);

% Load the image replay
load([replay_action_record_path, '\preprocessing_replay.mat'], "store_replay_steps");

%% For each image replay step, act the replay step
num_prep_steps = size(store_replay_steps, 2);
for step_i = 1:num_prep_steps
    operation = store_replay_steps{step_i};
    if strcmp(operation.type, 'create_slice_images')
        replay_create_slice_images(im_fn, replay_action_record_path, save_dir, operation);
    elseif strcmp(operation.type, 'rotate_center')
        replay_rotate_center(save_dir, operation);
    elseif strcmp(operation.type, 'reorder_slices')
        replay_reorder_slices(save_dir, operation);
    elseif strcmp(operation.type, 'flip_slices')
        replay_flip_slices(save_dir, operation);    
    end
end


function replay_create_slice_images(im_fn, replay_action_record_path, save_dir, operation)
    disp('replay_create_slice_images');
    if operation.slice_images == 1
        % Each image is a slice

    elseif operation.slice_images == 0
        % we have multiple images in a slice
        load([replay_action_record_path, '\slice_slide_locations.mat'], 'slice_slide_locations');
        assert(length(im_fn) == length(slice_slide_locations));
        assert(operation.downsample_factor == 1); % We only support dowmsample factor = 1 for now

        % Write all slice images to separate files
        curr_slice = 0;
        for curr_im = 1:length(im_fn)
            im_rgb = imread(im_fn{curr_im});
            for curr_slice_in_image = 1:length(slice_slide_locations{curr_im})
                curr_slice = curr_slice + 1;
                mask_y = slice_slide_locations{curr_im}{curr_slice_in_image}{1};
                mask_x = slice_slide_locations{curr_im}{curr_slice_in_image}{2};
                im_slice = im_rgb(mask_y, mask_x, :);
                curr_fn = fullfile(save_dir,sprintf('slice_%d.tif',curr_slice));
                imwrite(im_slice,curr_fn,'tif');
            end
        end
    end
    return;
end

function replay_rotate_center(save_dir, operation)
    disp('replay_rotate_center');

    slice_dir = dir(fullfile(save_dir,'*.tif'));
    slice_fn = natsortfiles(cellfun(@(path,fn) fullfile(path,fn), ...
        {slice_dir.folder},{slice_dir.name},'uni',false));
    
    slice_im = cell(length(slice_fn),1);
    for curr_slice = 1:length(slice_fn)
       slice_im{curr_slice} = imread(slice_fn{curr_slice});  
    end
    
    % Pad all slices centrally to the largest slice and make matrix
    slice_size_max = max(cell2mat(cellfun(@size,slice_im,'uni',false)),[],1);
    slice_im_pad = ...
        cell2mat(cellfun(@(x) x(1:slice_size_max(1),1:slice_size_max(2),:), ...
        reshape(cellfun(@(im) padarray(im, ...
        [ceil((slice_size_max(1) - size(im,1))./2), ...
        ceil((slice_size_max(2) - size(im,2))./2)],0,'both'), ...
        slice_im,'uni',false),1,1,1,[]),'uni',false));


    align_axis = operation.align_axis;
    % Get angle for all axes
    align_angle = squeeze(atan2d(diff(align_axis(:,1,:),[],1),diff(align_axis(:,2,:),[],1)));
    align_center = permute(nanmean(align_axis,1),[2,3,1]);
    
    % Set target angle as the nearest multiple of 90
    target_angle = round(nanmean(align_angle)/90)*90;
    
    % Set target position as the average center of the reference lines
    target_position = nanmean(align_center,2);
    
    im_aligned = zeros(size(slice_im_pad),class(slice_im_pad));
    for curr_im = 1:length(slice_im)
        angle_diff = target_angle - align_angle(curr_im);
        x_diff = target_position(1) - align_center(1,curr_im);
        y_diff = target_position(2) - align_center(2,curr_im);
        im_aligned(:,:,:,curr_im) = ...
            imrotate(imtranslate(slice_im_pad(:,:,:,curr_im), ...
            [x_diff,y_diff]),angle_diff,'bilinear','crop');
        
    end
    
    % Overwrite old images with new ones
    for curr_im = 1:size(im_aligned,4)
        imwrite(im_aligned(:,:,:,curr_im),slice_fn{curr_im},'tif');
    end
end

function replay_reorder_slices(save_dir, operation)
    disp('replay_reorder_slices');
    slice_idx = operation.slide_idx;
    
    slice_dir = dir(fullfile(save_dir,'*.tif'));
    slice_fn = natsortfiles(cellfun(@(path,fn) fullfile(path,fn), ...
        {slice_dir.folder},{slice_dir.name},'uni',false));
    
    % Make re-ordered filenames (with '_reorder to avoid overwriting)
    reordered_source_filenames = slice_fn(slice_idx);
    reorder_target_filenames = strrep(slice_fn,'.tif','_reorder.tif');
    
    % Rename files (with '_reorder')
    for curr_im = 1:length(slice_fn)
        movefile(reordered_source_filenames{curr_im},reorder_target_filenames{curr_im});
    end
    
    % Rename files (with original filenames)
    for curr_im = 1:length(slice_fn)
        movefile(reorder_target_filenames{curr_im},slice_fn{curr_im});
    end

end

function replay_flip_slices(save_dir, operation)
    disp('replay_flip_slices');
    is_fliplr = operation.is_fliplr;
    is_flipud = operation.is_flipud;
        
    slice_dir = dir(fullfile(save_dir,'*.tif'));
    slice_fn = natsortfiles(cellfun(@(path,fn) fullfile(path,fn), ...
        {slice_dir.folder},{slice_dir.name},'uni',false));
    
    assert(length(slice_fn) == length(is_fliplr));
    assert(length(slice_fn) == length(is_flipud));

    for curr_im = 1:length(slice_fn)
       im = imread(slice_fn{curr_im});
       if is_fliplr(curr_im) im = fliplr(im); end
       if is_flipud(curr_im) im = flipud(im); end
       imwrite(im, slice_fn{curr_im});
    end
       
end