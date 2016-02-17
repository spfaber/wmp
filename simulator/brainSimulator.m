%function [predDwiFile] = brainSimulator(dwiFile, bvecs, bvals, fgFile, noiseDist, diffModel)

% This simulator code is designed to:
%  - A - Take a fiber group/connectome.
%  - B - Build and fit the life model.
%  - C - Extract the model-predicted diffusion data and write to disk. 
%  - D - Pick major tracts to keep the signal.
%  - E - Track the new signal based on the major tracts.

%  2016/01/27 Sam Faber

% Build the file names for the diffusion data, bvecs, bval, fiber group, noise distribution, and model. 
subj = '105115';
rootDir   = fullfile('/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/mrtrix_track_between_rois/',subj);
saveDir   = fullfile('/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/simulator/',subj);
dwiFile   = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.nii.gz');
t1File    = fullfile(rootDir, 't1.nii.gz');
bvecsFile = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.bvecs');
bvalsFile = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.bvals');
bvecs     = dlmread(bvecsFile);
bvals     = dlmread(bvalsFile);


% read in fiber group if cleaning needs to be done first
fgr = fgRead(fullfile(rootDir,'mrtrix_results','right_optic_radiation_PCSD_updated.mat'));


% clean fibers or load in cleaned fiber group
disp('Cleaning fibers...');
[fgrClean, fibersToKeep] = mbaComputeFibersOutliers(fgr, 4, 4);
disp('done');
fgFileName    = fgrClean;
feFileName    = strcat('fe_rOR_simulator_',subj);
savedir = saveDir;


% Run LiFE (lines 39-54) or load in LiFE fe structure (line 55)
% initialize LiFE model 'fe' structure
disp('Building life model...');
N = 360; % Discretization parameter
fe = feConnectomeInit(dwiFile,fgFileName,feFileName,savedir,dwiFile,t1File,N,[1,0],0);
disp('done');


% fit the model
disp('Fitting life model...');
tic
fe = feSet(fe,'fit',feFitModel(feGet(fe,'model'),feGet(fe,'dsigdemeaned'),'bbnnls'));
toc
disp('done');


% save fe structure
save(fullfile(saveDir,feFileName),'-struct','fe','-v7.3'); 
fe = load(fullfile(saveDir,feFileName));


% read in the diffusion nifti
nii = niftiRead(dwiFile);


% rename diffusion file and set data to Nan everywhere (we only want signal
% from our fiber group)
niiOut = nii;


% Get the coordinates of the nodes in each voxel of the VOI (fiber
% group)
coords = fe.roi.coords;


% Get the measured diffusion signal in the voxels of the fiber group
% returns #bvecs X #voxels array
vals = feGet(fe,'dsiinvox',coords);   % dw signal in tract


% Keep the bvals that are not zero (images that are diffusion weighted)
% since this is the dimension of vals we get above for bvecs
indexes = find(bvals~=0);
b0indexes = find(bvals==0);


% b0_data has to be a matrix nB0 X nVoxels
b0_data = nan(size(b0indexes,2),size(coords,1));
for ivx = 1:size(coords,1)
    b0_data(:,ivx) = niiOut.data(coords(ivx,1),coords(ivx,2),coords(ivx,3),b0indexes);
end
niiOut.data = nan(size(niiOut.data));

% cooment me well
niiOut.data = feReplaceImageValues(niiOut.data,b0_data,coords,b0indexes);


% Replace Nans in voxels of fiber group with vals
niiOut.data = feReplaceImageValues(niiOut.data,vals,coords,indexes);

% Create WM mask of fiber group for tracking - 1's everywhere the fibers
% are, 0's everywhere else
niiWM = niiOut;
niiWM.data = zeros(size(niiWM.data));
valsMask = ones(size(vals));
niiWM.data = feReplaceImageValues(niiWM.data, valsMask, coords, indexes);


% Get original isotropic diffusion signal in each voxel
iso = (fe.life.diffusion_signal_img(:) -feGet(fe,'dsigdemeaned') );


% Set iso to the median isotropic signal in every OR voxel 
mu = 0;
sd = 0.001; % Use these three types of noise: 0.001, 1, 0.01, 0.1
iso = iso + (mu + sd .* randn(size(iso)));
setIso = median(iso);


% Generate predicted (demeaned) signal
pSig = feGet(fe,'pSig fiber');


% Add the isotropic component of the original signal to the predicted signal
fullPredOrigIso = pSig+iso;


% Add the standardized isotropic signal to the predicted signal
fullPredSetIso = pSig+setIso;


% Reshape the full signal so that the dimensions are (#bvals X #voxels)
fullPredOrigIso = reshape(fullPredOrigIso,size(coords,1),length(indexes));
fullPredSetIso  = reshape(fullPredSetIso, size(coords,1),length(indexes));


% Replace the nifti image values with the diffusion weighted bvals, coords of voxels, and full signal
niiOrigIso      = niiOut;
niiSetIso       = niiOut;
niiOrigIso.data = feReplaceImageValues(niiOrigIso.data,fullPredOrigIso',coords,indexes);
niiSetIso.data  = feReplaceImageValues(niiSetIso.data,fullPredSetIso',coords, indexes);


% plot all data in same slice containing the OR to see difference
map = 'bone';

figure(1)
imagesc(nii.data(:,:,15))
colormap(map)
title('Original Diffusion Signal');

figure(2)
imshow(niiOut.data(:,:,37))
colormap(map)
title('Original Diffusion Signal only in OR');

figure(3)
imagesc(niiOrigIso.data(:,:,37))
colormap(map)
title('Full Predicted Signal w/ original Iso only in OR');

figure(4)
imagesc(niiSetIso.data(:,:,37))
colormap(map)
title('Full Predicted Signal w/ set Iso only in OR');


% check WM data
figure(5)
imagesc(niiWM.data(:,:,37))
title('Fiber White Matter Mask');


% check using mbaDisplayBrainSlice
figure(6)
mbaDisplayBrainSlice(nii, [0 0 15]);
title('Original Signal');


% Make sure data type is correct and rename new nifti files to be used for tracking
niiWM.data       = int16(niiWM.data);
niiOrigIso.data  = int16(niiOrigIso.data);
niiSetIso.data   = int16(niiSetIso.data);
niiWM.fname      = fullfile(saveDir,'diffusion',sprintf('WM_mask_rOR_%s.nii.gz',subj));
niiOrigIso.fname = fullfile(saveDir,'diffusion',sprintf('sim_orig_aniso_diff_rOR_%s.nii.gz',subj));
niiSetIso.fname  = fullfile(saveDir,'diffusion',sprintf('sim_set_aniso_diff_rOR_%s.nii.gz',subj));


% Write the new nifti files to disk
niftiWrite(niiWM);
niftiWrite(niiOrigIso);
niftiWrite(niiSetIso);

% Check that the nifti files have the right data type - THEY DO NOT
checkWM      = niftiRead(fullfile(saveDir,'diffusion',sprintf('WM_mask_rOR_%s.nii.gz',subj)));
checkOrigIso = niftiRead(fullfile(saveDir,'diffusion',sprintf('sim_orig_aniso_diff_rOR_%s.nii.gz',subj)));
checkSetIso  = niftiRead(fullfile(saveDir,'diffusion',sprintf('sim_set_aniso_diff_rOR_%s.nii.gz',subj)));



