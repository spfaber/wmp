function s_dev_fe_make_wm_mask_hcp
%
% This script makes the white-matter mask used to track the connectomes in
% Pestilli et al., LIFE paper.

subjects = {''};
SUBJECTS_DIR = '/N/dc2/projects/lifebid/HCP/Sam/';
setenv('SUBJECTS_DIR', SUBJECTS_DIR);
anatomypath = getenv('SUBJECTS_DIR');

%savepath = '/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/mrtrix_track_between_rois/KW/mrtrix_results';

wmMaskFile = fullfile(SUBJECTS_DIR,subjects{1},'wm_mask_SAM.nii.gz');

% bOpen the segmentation file (the MGZ) save it to disk as a nifti (nii.gz)
% make sure that the file is oriented with the mrDiffusion conventions
% (RAS=Right-Anterior-Superior)
fs_wm = fullfile(anatomypath,subjects{1},'mri','aseg.mgz');
eval(sprintf('!mri_convert  --out_orientation RAS %s %s', fs_wm, wmMaskFile));

% We load the NIFTI file just created, this file indicates the structure of
% assignement for each voxel. Some voxels are assigned to be WM others
% ventricles, etc.
wm = niftiRead(wmMaskFile);

% We want to select voxels that have indices corresponding to the WM
% structures. These are the indices we like.
invals  = [2 41 16 17 28 60 51 53 12 52 13 18 54 50 11 251 252 253 254 255 10 49 46 7];

% FInally we are going to select only the voxles with indices of WM. USe
% these to make a new NIFTI file containin only WM voxels. This is our
% mask.
origvals = unique(wm.data(:));
fprintf('\n[%s] Converting voxels... ',mfilename);

wmCounter=0;noWMCounter=0;
for ii = 1:length(origvals);
    if any(origvals(ii) == invals)
        wm.data( wm.data == origvals(ii) ) = 1;
        wmCounter=wmCounter+1;
    else
        wm.data( wm.data == origvals(ii) ) = 0;
        noWMCounter = noWMCounter + 1;
    end
end
fprintf('converted %i regions to White-matter (%i regions left outside of WM)\n\n',wmCounter,noWMCounter);
niftiWrite(wm);

end % 