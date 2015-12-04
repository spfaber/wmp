% This script shows how to use dtiExportRoiToNifti function to create nifti
% from matlab ROI (mrDiffusion)

% Set up paths to matlab ROI and t1 reference image
path = '/N/dc2/projects/lifebid/HCP/Sam/matlab_code/pestillilab_projects/sam_faber/optic_radiation/mrtrix_track_between_rois/118730';
lhRoi = fullfile(path,'/ROIs','lh_LGN.mat');
rhRoi = fullfile(path,'/ROIs','rh_LGN.mat');
lhMergeRois = fullfile(path,'/ROIs','lh_LGN2V1.mat');
rhMergeRois = fullfile(path,'/ROIs','rh_LGN2V1.mat');
refImg = fullfile(path,'T1w_acpc_dc_restore_1p25.nii.gz');

% dtiRoiNiftiFromMat(matRoi,refImg, 1) % This function serves a similar
% purpose but is not always compatible with ROIs

% Export the .mat to .nii
dtiExportRoiToNifti(lhRoi, refImg, fullfile(path,'/ROIs','lh_LGN.nii.gz'));
dtiExportRoiToNifti(rhRoi, refImg, fullfile(path,'/ROIs','rh_LGN.nii.gz'));
dtiExportRoiToNifti(lhMergeRois, refImg, fullfile(path, '/ROIs','lh_LGN2V1.nii.gz'));
dtiExportRoiToNifti(rhMergeRois, refImg, fullfile(path, '/ROIs','rh_LGN2V1.nii.gz'));

