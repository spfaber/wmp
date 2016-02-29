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
rootDir   = fullfile('/N/dc2/projects/lifebid/HCP/Sam/',subj);
saveDir   = fullfile('/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/simulator/',subj);
TRKDIR = fullfile(saveDir, 'mrtrix_results');
dwiFile   = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.nii.gz');
t1File    = fullfile(rootDir,'anatomy','T1w_acpc_dc_restore_1p25.nii.gz');
bvecsFile = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.bvecs');
bvalsFile = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.bvals');
bvecs     = dlmread(bvecsFile);
bvals     = dlmread(bvalsFile);


% read in fiber group if cleaning needs to be done first
fgArcuate = fgRead(fullfile(saveDir,'fascicles','fg_right_arc_500000_segmented.pdb'));
fgCortico = fgRead(fullfile(saveDir,'fascicles','fg_right_cort_500000_segmented.pdb'));


% clean fibers or load in cleaned fiber group
disp('Cleaning fibers...');
[fgArcClean, fibersToKeep1] = mbaComputeFibersOutliers(fgArcuate, 3.5, 3.5);
disp('done');

disp('Cleaning fibers...');
[fgCortClean, fibersToKeep2] = mbaComputeFibersOutliers(fgCortico, 3.5, 3.5);
disp('done');

fgArcCortClean = fgMerge(fgArcClean,fgCortClean);
%fgWrite(fgArctCortClean);

%fgWrite(fgrClean, fullfile(TRKDIR,strcat(fgrClean.name,'_cleaned')),'mat');

fgFileName    = fgArcCortClean;
feFileName    = strcat('fe_rArCort_simulator_',subj);
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
save(fullfile(saveDir,'fe',feFileName),'-struct','fe','-v7.3'); 

fe = load('fe_rOR_simulator_105115.mat');
% read in the diffusion nifti
nii = niftiRead(dwiFile);

% rename diffusion file and set data to Nan everywhere (we only want signal
% from our fiber group)
niiOut = nii;


% Get the coordinates of the nodes in each voxel of the connectome
coords = fe.roi.coords;


% Get the measured diffusion signal in the voxels of the fiber group
% returns #bvecs X #voxels array
dw_vals = feGet(fe,'dsiinvox',coords);   % dw signal in tract

% Keep the bvals that are not zero (images that are diffusion weighted)
% since this is the dimension of vals we get above for bvecs
indexes = find(bvals~=0);
b0indexes = find(bvals==0);


% check b0 signal distributions
figure('color','w')
[y,x] = hist(b0_data(1,:),100);
plot(x,y,'k')
hold on
[y,x] = hist(b0_data(2,:),100);
plot(x,y,'r')
[y,x] = hist(b0_data(3,:),100);
plot(x,y,'g')


% make b0_data a matrix nB0 X nVoxels
b0_data = nan(size(b0indexes,2),size(coords,1));
for ivx = 1:size(coords,1)
    b0_data(:,ivx) = niiOut.data(coords(ivx,1),coords(ivx,2),coords(ivx,3),b0indexes);
end
niiOut.data = nan(size(niiOut.data));

% Replace Nans in voxels of fiber group with b0_data 
niiOut.data = feReplaceImageValues(niiOut.data,b0_data,coords,b0indexes);


% Replace Nans in voxels of fiber group with dw_vals
niiOut.data = feReplaceImageValues(niiOut.data,dw_vals,coords,indexes);

% Create WM mask of fiber group for tracking - 1's everywhere the fibers
% are, 0's everywhere else
niiWM = niiOut;
niiWM.data = zeros(size(niiWM.data));
b0_valsMask = ones(size(b0_data));
dw_valsMask = ones(size(dw_vals));
niiWM.data = feReplaceImageValues(niiWM.data, dw_valsMask, coords, indexes);
niiWM.data = feReplaceImageValues(niiWM.data, b0_valsMask, coords, b0indexes);


% Get original isotropic diffusion signal in each voxel
iso = (fe.life.diffusion_signal_img(:) -feGet(fe,'dsigdemeaned') );


% Set iso to the median isotropic signal in every OR voxel 
setIso = median(iso).*ones(size(iso)); % to use this, you must change MRtrix
% threshold in order to compute csdeconv

% Add random noise to isotropic signal in every OR voxel
mu = 0;
sd = 250; % Use these three types of noise: 0.001, 1, 0.01, 0.1
setRandIso = iso + (mu + sd.*randn(size(iso)));


% Generate predicted (demeaned) signal
pSig = feGet(fe,'pSig fiber');


% Add the isotropic component of the original signal to the predicted signal
% use feGet(fe, 'pSig full')
fullPredOrigIso    = pSig+iso;
fullPredSetRandIso = pSig+setRandIso;



% Add the standardized isotropic signal to the predicted signal
fullPredSetIso     = pSig+setIso;
fullPredSetRandIso = reshape(fullPredSetRandIso, size(coords,1),length(indexes));

% Reshape the full signal so that the dimensions are (#bvals X #voxels)
fullOrigSig        = fe.life.diffusion_signal_img;
fullPredOrigIso    = reshape(fullPredOrigIso,size(coords,1),length(indexes));
fullPredSetIso     = reshape(fullPredSetIso, size(coords,1),length(indexes));
fullPredSetRandIso = reshape(fullPredSetRandIso, size(coords,1),length(indexes));

% Replace the nifti image values with the diffusion weighted bvals, coords of voxels, and full signal
niiOrigIso      = niiOut;
niiSetIso       = niiOut;
niiSetRandIso   = niiOut;
niiOrigIso.data    = feReplaceImageValues(niiOrigIso.data,fullPredOrigIso',coords,indexes);
niiSetIso.data     = feReplaceImageValues(niiSetIso.data,fullPredSetIso',coords, indexes);
niiSetRandIso.data = feReplaceImageValues(niiSetRandIso.data, fullPredSetRandIso',coords, indexes);


% plot all data in same slice containing the OR to see difference
map = 'bone';

figure(2)
imagesc(niiOut.data(:,:,37))
colormap(map)
title('Original Diffusion Signal');

figure(2)
imshow(niiWM.data(:,:,37))
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
imshow(niiWM.data(:,:,37))
title('Fiber White Matter Mask');


% check using mbaDisplayBrainSlice
figure(6)
mbaDisplayBrainSlice(nii, [0 0 15]);
title('Original Signal');


% Make sure data type is correct and rename new nifti files to be used for tracking
niiOut.data        = int16(niiOut.data);
niiWM.data         = int16(niiWM.data);
niiOrigIso.data    = int16(niiOrigIso.data);
niiSetIso.data     = int16(niiSetIso.data);
niiSetRandIso.data = int16(niiSetRandIso.data);
niiOut.fname     = fullfile(saveDir,'diffusion',sprintf('measured_diff_sig_rOR_%s.nii.gz',subj));
niiWM.fname      = fullfile(saveDir,'diffusion',sprintf('WM_mask_rOR_%s.nii.gz',subj));
niiOrigIso.fname = fullfile(saveDir,'diffusion',sprintf('sim_orig_aniso_diff_rOR_%s.nii.gz',subj));
niiSetIso.fname  = fullfile(saveDir,'diffusion',sprintf('sim_set_aniso_diff_rOR_%s.nii.gz',subj));
niiSetRandIso.fname = fullfile(saveDir, 'diffusion',...
    sprintf('sim_set_rand_aniso_diff_rOR_%s.nii.gz',subj));

% Write the new nifti files to disk
niftiWrite(niiWM);
niftiWrite(niiOrigIso);
niftiWrite(niiSetIso);
niftiWrite(niiSetRandIso);

% Check that the nifti files have the right data type - THEY DO NOT
checkWM      = niftiRead(fullfile(saveDir,'diffusion',sprintf('WM_mask_rOR_%s.nii.gz',subj)));
checkOrigIso = niftiRead(fullfile(saveDir,'diffusion',sprintf('sim_orig_aniso_diff_rOR_%s.nii.gz',subj)));
checkSetIso  = niftiRead(fullfile(saveDir,'diffusion',sprintf('sim_set_aniso_diff_rOR_%s.nii.gz',subj)));

%% Plot the new fiber group tracked in MRtrix using the simulated signal

% convert .tck to .pdb
mrtrix_tck2pdb(fullfile(TRKDIR,'right_OR_sim_orig_aniso_PCSD.tck'),fullfile(TRKDIR,'right_OR_sim_orig_aniso_PCSD.pdb'))
mrtrix_tck2pdb(fullfile(TRKDIR,'right_OR_sim_set_aniso_PCSD.tck'),fullfile(TRKDIR,'right_OR_sim_set_aniso_PCSD.pdb'))

% read in .pdb and save as .mat
fgOrig = fgRead(fullfile(TRKDIR, 'right_OR_sim_orig_aniso_PCSD.pdb'));
fgWrite(fgOrig, fullfile(TRKDIR,fgOrig.name),'mat'); 

fgSet = fgRead(fullfile(TRKDIR, 'right_OR_sim_set_aniso_PCSD.pdb'));
fgWrite(fgSet, fullfile(TRKDIR,fgSet.name),'mat'); 

figure (1)
for ii = 1:length(fgArcClean.fibers);
    plot3(fgArcClean.fibers{ii}(1,:),fgArcClean.fibers{ii}(2,:),fgArcClean.fibers{ii}(3,:),'b'); hold on
end
hold on

for ii = 1:length(fgCortClean.fibers);
    plot3(fgCortClean.fibers{ii}(1,:),fgCortClean.fibers{ii}(2,:),fgCortClean.fibers{ii}(3,:),'r'); hold on
end
view(90,0)

figure (3)
for ii = 1:length(fgSet.fibers);
    plot3(fgSet.fibers{ii}(1,:),fgSet.fibers{ii}(2,:),fgSet.fibers{ii}(3,:)); hold on
end

% plot different signals together 

iso_sig = reshape(iso,size(coords,1), length(indexes));

figure(1)
plot(fullPredSetIso(1:20,1:10),'c');
hold on

plot(fullPredSetRandIso(1:20,1:10),'g');
hold on
plot(iso_sig(1:20,1:10),'b');

plot(fullPredOrigIso(1:20,1:10),'b');
plot(fullOrigSig(1:20,1:10),'r');


hold on 
legend('Original Full Sig','Full Pred Sig Set Iso',...
    'Full Pred Sig Orig Iso','Original Iso Sig');

mrtrix_tck2pdb(fullfile(TRKDIR,'right_OR_sim_orig_aniso_PCSD.tck'),fullfile(TRKDIR,'right_OR_sim_orig_aniso_PCSD.pdb'))
mrtrix_tck2pdb(fullfile(TRKDIR,'right_OR_sim_set_rand_aniso_PCSD.tck'),fullfile(TRKDIR,'right_OR_sim_set_rand_aniso_PCSD.pdb'))

fgOrig = fgRead(fullfile(TRKDIR,'right_optic_radiation_PCSD_updated_cleaned.mat'));
fgSimOrigIso = fgRead(fullfile(TRKDIR,'right_OR_sim_orig_aniso_PCSD.pdb'));
fgSimSetRandISo = fgRead(fullfile(TRKDIR,'right_OR_sim_set_rand_aniso_PCSD.pdb'));

figure (2)
for ii = 1:length(fgOrig.fibers);
    plot3(fgOrig.fibers{ii}(1,:),fgOrig.fibers{ii}(2,:),fgOrig.fibers{ii}(3,:),'g'); hold on
end


figure (3)
for ii = 1:length(fgSimOrigIso.fibers);
    plot3(fgSimOrigIso.fibers{ii}(1,:),fgSimOrigIso.fibers{ii}(2,:),fgSimOrigIso.fibers{ii}(3,:),'c'); hold on
end
view(90,0)

figure (4)
for ii = 1:length(fgSimSetRandISo.fibers);
    plot3(fgSimSetRandISo.fibers{ii}(1,:),fgSimSetRandISo.fibers{ii}(2,:),fgSimSetRandISo.fibers{ii}(3,:),'y'); hold on
end

