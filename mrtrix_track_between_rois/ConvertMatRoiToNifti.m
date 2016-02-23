% This script shows how to use dtiExportRoiToNifti function to create nifti
% from matlab ROI (mrDiffusion)

% Set up paths to matlab ROI and t1 reference image
path = '/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/mrtrix_track_between_rois/KW';
lhRoi = fullfile(path,'/ROIs','lh_LGN.mat');
rhRoi = fullfile(path,'/ROIs','rh_LGN.mat');
%lhMergeRois = fullfile(path,'/ROIs','lh_LGN2V1.mat');
%rhMergeRois = fullfile(path,'/ROIs','rh_LGN2V1.mat');
refImg = fullfile(path,'t1.nii.gz');

% If you want to merge ROIs computationally
% lhLGN = load(fullfile(path,'/ROIs','lh_LGN.mat'));
% lhV1 = load(fullfile(path,'/ROIs','lh_V1_label_smooth3mm_ROI.mat'));
% MergeLGNV1 = [lhRoi lhV1];

% dtiRoiNiftiFromMat(matRoi,refImg, 1) % This function serves a similar
% purpose but is not always compatible with ROIs

% Export the .mat to .nii
dtiExportRoiToNifti(lhRoi, refImg, fullfile(path,'/ROIs','lh_LGN.nii.gz'));
dtiExportRoiToNifti(rhRoi, refImg, fullfile(path,'/ROIs','rh_LGN.nii.gz'));
dtiExportRoiToNifti(lhMergeRois, refImg, fullfile(path, '/ROIs','lh_LGN2V1.nii.gz'));
dtiExportRoiToNifti(rhMergeRois, refImg, fullfile(path, '/ROIs','rh_LGN2V1.nii.gz'));

