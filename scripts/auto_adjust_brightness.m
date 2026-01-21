

% Script to process all TIFF images in a folder
% Usage: process_folder_tiff_images(folder_path)
% If folder_path is not provided, it will prompt for user input

folder_path = "C:\Users\benji\Documents\harini\8504A_iso\cfos\batch1\slices";

% Check if folder exists
if ~isfolder(folder_path)
    error('Folder does not exist: %s', folder_path);
end

% Find all TIFF files (both .tif and .tiff extensions)
tiff_files = [dir(fullfile(folder_path, '*.tif')); dir(fullfile(folder_path, '*.tiff'))];

if isempty(tiff_files)
    fprintf('No TIFF files found in folder: %s\n', folder_path);
    return;
end

fprintf('Found %d TIFF file(s). Processing...\n', length(tiff_files));

% Process each TIFF file
for i = 1:length(tiff_files)
    file_path = fullfile(folder_path, tiff_files(i).name);
    fprintf('Processing: %s\n', tiff_files(i).name);
    
    try
        % Read the image
        img_raw = imread(file_path);
        
        % Convert to double if needed (for proper quantile calculation)
        if ~isa(img_raw, 'double')
            img_raw = double(img_raw);
            if max(img_raw(:)) > 1
                img_raw = img_raw / max(img_raw(:)); % Normalize to [0, 1] range
            end
        end
        
        % Apply brightness adjustment
        img_adjusted = adjust_brightness(img_raw);
        
        % Convert back to appropriate data type for saving
        if max(img_adjusted(:)) <= 1
            % If image is in [0, 1] range, convert to uint16 for TIFF
            img_adjusted = uint16(img_adjusted * 65535);
        else
            img_adjusted = uint16(img_adjusted);
        end
        
        % Save the adjusted image, replacing the original
        imwrite(img_adjusted, file_path, 'tif');
        fprintf('  Saved: %s\n', tiff_files(i).name);
        
    catch ME
        fprintf('  Error processing %s: %s\n', tiff_files(i).name, ME.message);
    end
end

fprintf('Processing complete!\n');



% adjust the contrast and brightness of an image adaptively to make it human eye friendly
function img = adjust_brightness(img_raw)
    img = img_raw;
    img_hist = img_raw(:);
    
    % Display quantiles for debugging
    % fprintf('10th percentile: %f\n', quantile(img_hist(img_hist > 0), 0.1));
    % fprintf('90th percentile: %f\n', quantile(img_hist(img_hist > 0), 0.9));
    
    % this is not good for visualization, try to normalize the image so that the quantile 0.1 and 0.9 are 0 and 1
    img = (img - quantile(img_hist(img_hist > 0), 0.6)) / quantile(img_hist(img_hist > 0), 0.99) * 1;
    img(img < 0) = 0;
    img(img > 1) = 1;
    % img = normalize(img, 0, 255, 'range'); % Alternative normalization if needed
end