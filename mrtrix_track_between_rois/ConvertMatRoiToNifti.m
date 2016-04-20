% This script shows how to use dtiExportRoiToNifti function to create nifti
% from matlab ROI (mrDiffusion)

% Set up paths to matlab ROI and t1 reference image
path = '/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/simulator/105115/';
Arcdors = fullfile(path,'/ROIs','rArc_dorsal_endpointROI.mat');
Arcvent = fullfile(path,'/ROIs','rArc_ventral_endpointROI.mat');
Cortsup = fullfile(path,'/ROIs','rCort_sup_endpointROI.mat');
Cortinf = fullfile(path,'/ROIs','rCort_inf_endpointROI.mat');


%rhMergeRois = fullfile(path,'/ROIs','rh_LGN2V1.mat');
anatpath = '/N/dc2/projects/lifebid/HCP/Sam/105115/anatomy/';
refImg = fullfile(anatpath,'T1w_acpc_dc_restore_1p25.nii.gz');

% If you want to merge ROIs computationally
% lhLGN = load(fullfile(path,'/ROIs','lh_LGN.mat'));
% lhV1 = load(fullfile(path,'/ROIs','lh_V1_label_smooth3mm_ROI.mat'));
% MergeLGNV1 = [lhRoi lhV1];

% dtiRoiNiftiFromMat(matRoi,refImg, 1) % This function serves a similar
% purpose but is not always compatible with ROIs

% Export the .mat to .nii
dtiExportRoiToNifti(Arcdors, refImg, fullfile(path,'/ROIs','rArc_dorsal_endpointROI.nii.gz'));
dtiExportRoiToNifti(Arcvent, refImg, fullfile(path,'/ROIs','rArc_ventral_endpointROI.nii.gz'));
dtiExportRoiToNifti(Cortsup, refImg, fullfile(path, '/ROIs','rCort_sup_endpointROI.nii.gz'));
dtiExportRoiToNifti(Cortinf, refImg, fullfile(path, '/ROIs','rCort_inf_endpointROI.nii.gz'));

