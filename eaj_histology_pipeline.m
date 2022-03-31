addpath(genpath('C:\Users\emily\OneDrive - Stanford\GitHub\External\AP_histology'))
addpath(genpath('C:\Users\emily\OneDrive - Stanford\GitHub\External\npy-matlab'))
addpath('C:\Users\emily\OneDrive - Stanford\GitHub\GiocomoLab\allenCCF\Browsing Functions')

%% 1) Load CCF and set paths for slide and slice images

% Load CCF atlas
allen_atlas_path = 'C:\Users\emily\OneDrive - Stanford\GitHub\External\AP_histology';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

% Set paths for histology images and directory to save slice/alignment
im_path = 'Z:\WT Sequences\2021_pilot\Histology\StickPin_DAPI_DID_5x';
slice_path = [im_path filesep 'atlas_alignment'];

%% 2) Preprocess slide images to produce slice images
% Set white balance and resize slide images, extract slice images

% Preprocess images
AP_process_histology(im_path);

% (optional) Rotate, center, pad, flip slice images
AP_rotate_histology(slice_path);

%% 3) Align CCF to slices

% Find CCF slices corresponding to each histology slice
AP_grab_histology_ccf(tv,av,st,slice_path);

% Align CCF slices and histology slices
% (first: automatically, by outline)
AP_auto_align_histology_ccf(slice_path);
% (second: curate manually)
AP_manual_align_histology_ccf(tv,av,st,slice_path);


%% 4) Map probe and save for further processing in Python

% Get probe trajectory from histology, convert to CCF coordinates
AP_get_probe_histology(tv,av,st,slice_path);

% Add relevant variables and save (EAJ addition, 20220329)
% map area indices to regions and layers
probe_ccf.areas = st.safe_name(probe_ccf.trajectory_areas);
probe_ccf.regions = char(strcat(extractBefore(probe_ccf.areas, 'area'), 'area'));
probe_ccf.layers = char(extractAfter(probe_ccf.areas, 'layer '));

% shift to bregma, convert to microns
probe_ccf.trajectory_bregma = (probe_ccf.trajectory_coords-allenCCFbregma())*10;
% shift DV so 0=brain surface (instead of skull surface at bregma)
% while this technically isn't bregma space, it's more useful
probe_ccf.trajectory_bregma(:,2) = ...
    probe_ccf.trajectory_bregma(:,2)-probe_ccf.trajectory_bregma(1,2);
% sites behind bregma in AP should be negative
probe_ccf.trajectory_bregma(:,1) = -1*probe_ccf.trajectory_bregma(:,1);
% re-order CCF APDVML to bregma MLAPDV
probe_ccf.trajectory_bregma = [probe_ccf.trajectory_bregma(:,3)...
    probe_ccf.trajectory_bregma(:,1)...
    probe_ccf.trajectory_bregma(:,2)];
% calculate angle of insertion
probe_ccf.angle = atan2d(probe_ccf.trajectory_bregma(end,3)...
    ,probe_ccf.trajectory_bregma(end,2)-probe_ccf.trajectory_bregma(1,2))-90;
save('probe_ccf.mat', 'probe_ccf')
