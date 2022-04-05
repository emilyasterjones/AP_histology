% each line must be run individually, as the functions 'end' while the GUI persists
% start by loading necessary paths
addpath(genpath('C:\Users\emily\OneDrive - Stanford\GitHub\External\AP_histology'))
addpath(genpath('C:\Users\emily\OneDrive - Stanford\GitHub\External\npy-matlab'))
addpath('C:\Users\emily\OneDrive - Stanford\GitHub\GiocomoLab\allenCCF\Browsing Functions')

%% 1) Load CCF and set paths for slide and slice images

% Load CCF atlas
% Download these files from http://data.cortexlab.net/allenCCF/
% Too large to be backed up to the Github repo
allen_atlas_path = 'C:\Users\emily\OneDrive - Stanford\GitHub\External\AP_histology';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

% Set paths for histology images and directory to save slice/alignment
im_path = 'Z:\WT Sequences\2022_winter\Histology\StrappyReeboks_DAPI_DiD_5x';
slice_path = [im_path filesep 'atlas_alignment'];

%% 2) Preprocess slide images to produce slice images
% You can run this, you'll get a better result if you just crop your slides
% into slices and rotate them how you like in image editing software. 
% For sagittal sections, orient the rostral end pointing left and the
% dental granule cell layer parallel to the side of the image (the same way
% it looks in the mouse atlas). Not required, but makes registration 100x
% easier.

% % Set white balance and resize slide images, extract slice images
% AP_process_histology(im_path);
% 
% % (optional) Rotate, center, pad, flip slice images
% AP_rotate_histology(slice_path);

%% 3) Align CCF to slices

% Find CCF slices corresponding to each histology slice
AP_grab_histology_ccf(tv,av,st,slice_path);

% Align CCF slices and histology slices
% (first: automatically, by outline)
% this has yet to work for me for sagittal slices, so I always run the next line
AP_auto_align_histology_ccf(slice_path);
% (second: curate manually)
% Turn off the overlay view. Click a point (1) on the image and (2) on the
% atlas. Repeat. Check overlay every so often until your target region at
% least matches your slide.
AP_manual_align_histology_ccf(tv,av,st,slice_path);


%% 4) Map probe and save for further processing in Python

% Get probe trajectory from histology, convert to CCF coordinates
% No need to be precise about the ends of the probe trajectory; this script
% just draws a straight line through the whole brain volume, averaging
% across slices.
AP_get_probe_histology(tv,av,st,slice_path);

% Add relevant variables and save (EAJ addition, 20220329)
% map area indices to regions and layers
load([slice_path filesep 'probe_ccf.mat'])
probe_ccf.areas = st.safe_name(probe_ccf.trajectory_areas);
for i=1:length(probe_ccf.areas)
    if contains(probe_ccf.areas(i), 'layer')
        probe_ccf.regions(i,:) = extractBefore(probe_ccf.areas(i), ' layer');
    else
        probe_ccf.regions(i,:) = probe_ccf.areas(i);
    end
end
% cast to char so it's readable when loaded into Python
probe_ccf.regions = char(probe_ccf.regions);
probe_ccf.layers = char(extractAfter(probe_ccf.areas, 'layer '));

% shift to bregma, convert to microns
bregma = allenCCFbregma();
probe_ccf.trajectory_bregma = (probe_ccf.trajectory_coords-bregma)*10;
% shift DV so 0=brain surface (instead of skull surface at bregma)
% while this technically isn't bregma space, it's more useful
probe_ccf.trajectory_bregma(:,2) = ...
    probe_ccf.trajectory_bregma(:,2)-probe_ccf.trajectory_bregma(1,2);
% sites behind bregma in AP should be negative
% I keep registering sites to the RH by accident when my insertions are all
% in LH, so this fixes that
probe_ccf.trajectory_bregma(:,1) = -1*probe_ccf.trajectory_bregma(:,1);
probe_ccf.trajectory_coords(:,3) = bregma(3) - (probe_ccf.trajectory_coords(:,3)-bregma(3));
% re-order CCF APDVML to bregma MLAPDV
probe_ccf.trajectory_bregma = [probe_ccf.trajectory_bregma(:,3)...
    probe_ccf.trajectory_bregma(:,1)...
    probe_ccf.trajectory_bregma(:,2)];

% calculate angle of insertion
% note this is much shallower or even in the opposite direction than your
% surgery angle, as the Allen atlas coordinates are relative to their
% sectioning plane, not to bregma space
% so not terribly useful
probe_ccf.angle = atan2d(probe_ccf.trajectory_bregma(end,3)...
    ,probe_ccf.trajectory_bregma(end,2)-probe_ccf.trajectory_bregma(1,2))-90;
save([slice_path filesep 'probe_ccf.mat'], 'probe_ccf')

%% 5. Save brain mesh for figures
addpath(genpath('C:\Users\emily\OneDrive - Stanford\GitHub\External\export_fig'))
% 1. load each probe_ccf into the cell array all_probes.mat (excludes CuffLinks)
% 2. subset trajectory_coords to stop at end of probe
% 3. run this script to generate the wire mesh plot with all probes plotted
plot_probe_3D(tv, probes)
export_fig all_probes.eps