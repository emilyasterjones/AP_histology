function AP_process_histology(im_path)
% AP_process_histology(im_path,resize_factor,slice_images)
%
% im_path - path with images of slides (tif/tiff/ome.tiff)
% resize_factor (if not ome.tiff) - resizing factor for saving images (e.g.
% 1/10 scales slides to 1/10 size). Note: if ome.tiff, microns per pixel is
% grabbed from image headers and the images are resized to match the CCF
% scaling (which is 10 microns per pixel)
% slice_images - true/false: images are individual slices, skip slice
% choosing step (false by default)
%
% Resize and white balance histology images and extract images of each slice
% Andy Peters (peters.andrew.j@gmail.com)

% Get and sort image files
im_path_dir = dir([im_path filesep '*5x.jpg*']);
im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
    {im_path_dir.folder},{im_path_dir.name},'uni',false));
im_info = imfinfo(im_fn{1});

% Set resize factor to match to Allen CCF
im_um2px = 0.908;
allen_um2px = 10; % Allen CCF: 10 um/voxel
resize_factor = im_um2px/allen_um2px;

%% Load and resize images
n_im = length(im_fn);
h = waitbar(0,'Loading and resizing images...');
im_rgb = cell(n_im,1);
for curr_im = 1:n_im
    im_rgb{curr_im} = imresize(imread(im_fn{curr_im}),resize_factor);
    waitbar(curr_im/n_im,h,['Loading and resizing images (' num2str(curr_im) '/' num2str(n_im) ')...']);
end
close(h)

%% set up GUI to pick slices on slide to extract 
slice_fig = figure('KeyPressFcn',@slice_keypress);

% Initialize data
slice_data = struct;
slice_data.im_path = im_path;
slice_data.im_fn = im_fn;
slice_data.im_rescale_factor = resize_factor;
slice_data.im_rgb = im_rgb;
slice_data.curr_slide = 0;
slice_data.slice_mask = cell(0,0);
slice_data.slice_rgb = cell(0,0);

% Update gui data
guidata(slice_fig, slice_data);

% Update slide
update_slide(slice_fig);
    
end

function slice_click(slice_fig,eventdata)
% On slice click, mark to extract

slice_data = guidata(slice_fig);

if eventdata.Button == 1
    
    selected_slice_bw = bwselect(slice_data.mask,eventdata.IntersectionPoint(1),eventdata.IntersectionPoint(2));
    
    % If the selected slice is already part of a user mask, delete that ROI
    if size(slice_data.user_masks,3) > 0
        clicked_mask = false(size(slice_data.mask));
        clicked_mask(round(eventdata.IntersectionPoint(2)),round(eventdata.IntersectionPoint(1))) = true;
        overlap_roi = any(clicked_mask(:) & reshape(slice_data.user_masks,[],size(slice_data.user_masks,3)),1);
        if any(overlap_roi)
            % Clear overlapping mask
            slice_data.user_masks(:,:,overlap_roi) = [];
            
            % Delete and clear bounding box
            delete(slice_data.user_rectangles(overlap_roi));
            slice_data.user_rectangles(overlap_roi) = [];
            
            % Update gui data
            guidata(slice_fig, slice_data);
            return
        end
    end
    
    % If left button pressed, create new slice ROI
    roi_num = size(slice_data.user_masks,3) + 1;
    
    % Make new mask with object
    slice_data.user_masks(:,:,roi_num) = selected_slice_bw;
    
    % Draw bounding box around object
    box_x = find(any(slice_data.user_masks(:,:,roi_num),1),1);
    box_y = find(any(slice_data.user_masks(:,:,roi_num),2),1);
    box_w = find(any(slice_data.user_masks(:,:,roi_num),1),1,'last') - box_x;
    box_h = find(any(slice_data.user_masks(:,:,roi_num),2),1,'last') - box_y;
    slice_data.user_rectangles(roi_num) = ...
        rectangle('Position',[box_x,box_y,box_w,box_h],'EdgeColor','w');
    
elseif eventdata.Button == 3
    % If right button pressed, manually draw rectangle ROI
    roi_num = size(slice_data.user_masks,3) + 1;
    
    % Draw ROI
    manual_roi = imrect;
    
    % Make new mask with object
    slice_data.user_masks(:,:,roi_num) = manual_roi.createMask;
    
    % Draw bounding box
    slice_data.user_rectangles(roi_num) = ...
        rectangle('Position',manual_roi.getPosition,'EdgeColor','w');
    
    % Delete the ROI
    manual_roi.delete;
    
end

% Update gui data
guidata(slice_fig, slice_data);

end

function slice_keypress(slice_fig,eventdata)
% Move to next slide with spacebar

if strcmp(eventdata.Key,'space')
    update_slide(slice_fig)
end

end

function update_slide(slice_fig)
% Find slices on slide by over-threshold objects of a large enough size

slice_data = guidata(slice_fig);

% Pull the images from selected slices (not during initialization)
if slice_data.curr_slide > 0
    extract_slice_rgb(slice_fig);
    slice_data = guidata(slice_fig);
end

% After the last slice, save the images and close out
if slice_data.curr_slide == length(slice_data.im_rgb)
    save_slice_rgb(slice_fig);
    close(slice_fig);
    return
end

slice_data.curr_slide = slice_data.curr_slide + 1;

% Minimum slice size
min_slice = 1000; %(1000/10)^2; % (um/10(CCF units))^2

% Estimate slice white threshold
curr_im_bw = nanmean(slice_data.im_rgb{slice_data.curr_slide},3);
[im_hist,im_hist_edges] = histcounts(curr_im_bw, ...
    linspace(min(curr_im_bw(:)),max(curr_im_bw(:)),100));
im_hist_deriv = [0;diff(smooth(im_hist,3))];
[~,bg_down] = min(im_hist_deriv);
bg_signal_min = find(im_hist_deriv(bg_down:end) > 0,1) + bg_down;
slice_threshold = im_hist_edges(bg_signal_min);

slice_mask = imfill(bwareaopen(mean( ...
    slice_data.im_rgb{slice_data.curr_slide},3) > slice_threshold,min_slice),'holes');
slice_conncomp = bwconncomp(slice_mask);

im_handle = imshow(slice_data.im_rgb{slice_data.curr_slide});
set(im_handle,'ButtonDownFcn',@slice_click);
title('Finding slice boundaries...');
drawnow;

slice_boundaries = bwboundaries(slice_mask);
slice_lines = gobjects(length(slice_boundaries),1);
for curr_slice = 1:length(slice_boundaries)
    slice_lines(curr_slice) = line(slice_boundaries{curr_slice}(:,2), ...
        slice_boundaries{curr_slice}(:,1),'color','w','linewidth',2,'LineSmoothing','on','linestyle','--');
end
title('Click to save/remove (left = auto, right = manual), spacebar to finish slide');

slice_data.im_h = im_handle;
slice_data.mask = slice_mask;
slice_data.lines = slice_lines;
slice_data.user_masks = zeros(size(slice_mask,1),size(slice_mask,2),0,'logical');
slice_data.user_rectangles = gobjects(0);

% Update gui data
guidata(slice_fig, slice_data);

end


function extract_slice_rgb(slice_fig)
% When changing slide, extract the selected slice images

slice_data = guidata(slice_fig);

n_slices = size(slice_data.user_masks,3);
curr_slice_mask = cell(n_slices,1);
curr_slice_rgb = cell(n_slices,1);
for curr_slice = 1:n_slices
    % Pull a rectangular area, exclude spaces (e.g. between torn piece)
    dilate_size = 30;
    curr_mask = imdilate(logical(any(slice_data.user_masks(:,:,curr_slice),2).* ...
        any(slice_data.user_masks(:,:,curr_slice),1)),ones(dilate_size));
    
    curr_rgb = reshape(slice_data.im_rgb{slice_data.curr_slide}( ...
        repmat(curr_mask,1,3)),sum(any(curr_mask,2)),sum(any(curr_mask,1)),3);
    
    curr_slice_mask{curr_slice} = curr_mask;
    curr_slice_rgb{curr_slice} = curr_rgb;
    
end

% Store the image and mask for each slice
slice_data.slice_mask{slice_data.curr_slide} = curr_slice_mask;
slice_data.slice_rgb{slice_data.curr_slide} = curr_slice_rgb;

% Update gui data
guidata(slice_fig, slice_data);

end


function save_slice_rgb(slice_fig)
% After the last slide, save the slice images

slice_data = guidata(slice_fig);

% Set save directory as subdirectory within original
save_dir = [slice_data.im_path filesep 'atlas_alignment'];
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end

% Concatenate all slice images
slice_rgb_cat = vertcat(slice_data.slice_rgb{:});

% Write all slice images to separate files
for curr_im = 1:length(slice_rgb_cat)
    curr_fn = [save_dir filesep num2str(curr_im) '.tif'];
    imwrite(slice_rgb_cat{curr_im},curr_fn,'tif');
end

% Get rows and columns for each slice corresponding to full size image
slice_slide_locations = cell(size(slice_data.slice_mask));
for curr_slide = 1:length(slice_data.slice_mask)
    for curr_slice = 1:length(slice_data.slice_mask{curr_slide})
        
        curr_mask = slice_data.slice_mask{curr_slide}{curr_slice};
        
        mask_x = find(interp1(1:size(curr_mask,2),+any(curr_mask,1), ...
            linspace(1,size(curr_mask,2), ...
            round(size(curr_mask,2)/slice_data.im_rescale_factor)),'nearest'));
        mask_y = find(interp1(1:size(curr_mask,1),+any(curr_mask,2), ...
            linspace(1,size(curr_mask,1), ...
            round(size(curr_mask,1)/slice_data.im_rescale_factor)),'nearest'));
        
        slice_slide_locations{curr_slide}{curr_slice} = ...
            {mask_y,mask_x};
        
    end
end

slice_slide_locations_fn = [save_dir filesep 'slice_slide_locations.mat'];
save(slice_slide_locations_fn,'slice_slide_locations');

disp(['Slices saved in ' save_dir]);

end



