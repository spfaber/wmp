% This simulator code is designed to:
%  - A - Take a fiber group/connectome.
%  - B - Build and fit the life model.
%  - C - Extract the model-predicted diffusion data and write to disk. 

%  2016/02/23 Sam Faber

% Set up paths and build the file names for the diffusion data, bvecs, bval, and fiber group. 
subj = '105115';
rootDir   = fullfile('/N/dc2/projects/lifebid/HCP/Sam/',subj);
saveDir   = fullfile('/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/simulator/',subj);
dwiFile   = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.nii.gz');
t1File    = fullfile(rootDir,'anatomy','T1w_acpc_dc_restore_1p25.nii.gz');
bvecsFile = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.bvecs');
bvalsFile = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.bvals');
bvecs     = dlmread(bvecsFile);
bvals     = dlmread(bvalsFile);

%% load the fe structure
fe = load(fullfile(saveDir,'fe','fe_rArCort_simulator_105115.mat'));
% load in fiber groups
fgArcClean     = fgRead(fullfile(saveDir,'fascicles', strcat('fgrArcClean',subj,'.mat')));
fgCortClean    = fgRead(fullfile(saveDir,'fascicles', strcat('fgrCortClean',subj,'.mat')));
fgArcCortClean = fgRead(fullfile(saveDir,'fascicles', strcat('fgrArcCortClean',subj,'.mat')));


% find fiber weights
weights = feGet(fe,'fiber weights');
weights_arc = weights(1:length(fgArcClean.fibers));
weights_cort = weights(1+length(fgArcClean.fibers):end);
%[y] = M_times_w(M,w);


% change weights
alpha = 0.999;
beta  = 1-alpha;
new_weight_arc = alpha*weights_arc;
new_weight_cort = beta*weights_cort;

% now multiply new weights times fibers to get new signal
new_weights = [new_weight_arc; new_weight_cort];

%new_signal = zeros(size(fe.life.M.Phi,1),size(fe.life.M.Phi,2),1);
new_signal = M_times_w( fe.life.M, new_weights);

% need to get new signal for just arc and just cort too
new_weight_arc  = [new_weight_arc; zeros(size(new_weight_cort))];
new_weight_cort = [zeros(size(weights_arc)); new_weight_cort];
new_sig_arc     = M_times_w(fe.life.M, new_weight_arc);
new_sig_cort    = M_times_w(fe.life.M, new_weight_cort);
new_sig_arc  = reshape(new_sig_arc, size(coords,1), length(indexes));
new_sig_cort = reshape(new_sig_cort, size(coords,1), length(indexes));

% read in the diffusion nifti
nii = niftiRead(dwiFile);

% rename diffusion file and set data to Nan everywhere (we only want signal
% from our fiber group)
niiOut  = nii;
niiArc  = nii;
niiCort = nii;

% get coords of voxels with overlapping fibers
% first, transform to image coords
fg_img = feGet(fe,'fg img');
fgArc_img = fg_img;
fgCrt_img = fg_img;
fgArc_img.fibers = fg_img.fibers(1:size(fgArcClean.fibers));
fgCort_img.fibers = fg_img.fibers(1+size(fgArcClean.fibers):end);

% then find unique coords of each fiber group
arc_coords  = fgGet(fgArc_img,'unique image coords');
cort_coords = fgGet(fgCort_img, 'unique image coords');

% then find the shared group of coords
test_int = intersect(arc_coords, cort_coords, 'rows','stable');
coords = fe.roi.coords;
common1Coords = intersect(coords,cort_coords,'rows','stable');
common2Coords = intersect(coords,arc_coords,'rows','stable');
commonAllCoords = intersect(common1Coords,common2Coords,'rows','stable');


% Get the measured diffusion signal in the voxels of the fiber group
% returns #bvecs X #voxels array
dw_vals = feGet(fe,'dsiinvox',commonAllCoords);   % dw signal in common voxels
%dw_arc  = feGet(fe,'dsiinvox',arc_coords);        % dw signal in arcuate
%dw_cort = feGet(fe,'dsiinvox',cort_coords);       % dw signal in cortico

% Keep the bvals that are not zero (images that are diffusion weighted)
% since this is the dimension of vals we get above for bvecs
indexes = find(bvals~=0);
b0indexes = find(bvals==0);


% make b0_data a matrix nB0 X nVoxels
b0_data = nan(size(b0indexes,2),size(coords,1));
for ivx = 1:size(coords,1)
    b0_data(:,ivx) = niiOut.data(coords(ivx,1),coords(ivx,2),coords(ivx,3),b0indexes);
end
niiOut.data  = nan(size(niiOut.data));
niiArc.data  = nan(size(niiArc.data));
niiCort.data = nan(size(niiCort.data));

% check b0 signal distributions
figure('color','w')
[y,x] = hist(b0_data(1,:),100);
plot(x,y,'k')
hold on
[y,x] = hist(b0_data(2,:),100);
plot(x,y,'r')
[y,x] = hist(b0_data(3,:),100);
plot(x,y,'g')




% Replace Nans in voxels of fiber group with b0_data 
niiOut.data  = feReplaceImageValues(niiOut.data,b0_data,commonAllCoords,b0indexes);
niiArc.data  = feReplaceImageValues(niiOut.data,b0_data,coords,b0indexes);
niiCort.data = feReplaceImageValues(niiOut.data,b0_data,coords,b0indexes);

% Replace Nans in voxels of fiber group with dw_vals
niiOut.data  = feReplaceImageValues(niiOut.data,dw_vals,commonAllCoords,indexes);
niiArc.data  = feReplaceImageValues(niiArc.data,new_sig_arc',coords,indexes);
niiCort.data = feReplaceImageValues(niiCort.data,new_sig_cort',coords,indexes);

% THIS IS THE ORIGINAL NII DATA ONLY IN THE VOXELS OF INTEREST

% Create WM mask of fiber group for tracking - 1's everywhere the fibers
% are, 0's everywhere else
niiWM = niiOut;
niiWM.data = zeros(size(niiWM.data));
b0_valsMask = ones(size(b0_data));
dw_valsMask = ones(size(dw_vals));
niiWM.data = feReplaceImageValues(niiWM.data, dw_valsMask, commonAllCoords, indexes);
niiWM.data = feReplaceImageValues(niiWM.data, b0_valsMask, commonAllCoords, b0indexes);

map = 'bone';
 % to check WM mask, plot slices 43 - 51 or so, 47 looks like it has most
 % voxels
for ii = 1:size(niiWM.data,3)
    figure(ii)
imagesc(niiWM.data(:,:,ii))
colormap(map)
end


figure(1) 
imagesc(niiWM.data(:,:,40))
colormap(map)


figure(2)
imagesc(niiArc.data(:,:,50));
colormap(map)

figure(3)
imagesc(niiCort.data(:,:,50));
colormap(map)


figure(4)
mbaDisplayBrainSlice(niiArc, [ 0 0 37]);
title('Original Signal');

% Get original isotropic diffusion signal in each voxel
iso = (fe.life.diffusion_signal_img(:) -feGet(fe,'dsigdemeaned') );


% Generate predicted (demeaned) signal
pSig = feGet(fe,'pSig fiber');


% Add the isotropic component of the original signal to the predicted signal
% use feGet(fe, 'pSig full')
fullPredSig      = pSig+iso;
fullPredSetWtSig = new_signal + iso; 

% Reshape the full signal so that the dimensions are (#bvals X #voxels)
fullOrigSig      = fe.life.diffusion_signal_img; 
fullPredSig      = reshape(fullPredSig, size(coords,1),length(indexes));
fullPredSetWtSig = reshape(fullPredSetWtSig, size(coords,1),length(indexes));
  
%ArcSetWtSig     = ;
%CortSetWtSig    = ;
%ArcCortSetWtSig = ;

% Replace the nifti image values with the diffusion weighted bvals, coords of voxels, and full signal
niiOrigIso         = niiOut;
niiOrigIso.data    = feReplaceImageValues(niiOrigIso.data,fullPredOrigIso',coords,indexes);

% Make sure data type is correct and rename new nifti files to be used for tracking
%niiWM.data         = int16(niiWM.data);
niiArc.data = int16(niiArc.data);
niiCort.data = int16(niiCort.data);

niiOrigIso.data    = int16(niiOrigIso.data);
%niiWM.fname      = fullfile(saveDir,'diffusion',sprintf('WM_mask_rArCort_%s.nii.gz',subj));
niiOrigIso.fname = fullfile(saveDir,'diffusion',sprintf('sim_orig_aniso_diff_rOR_%s.nii.gz',subj));
niiArc.fname = fullfile(saveDir, 'diffusion' ,sprintf('setwt_Rarc_sig',subj));
niiCort.fname = fullfile(saveDir, 'diffusion' ,sprintf('setwt_Rcort_sig',subj));

% Write the new nifti files to disk
%niftiWrite(niiWM);
niftiWrite(niiCort);
niftiWrite(niiArc);
niftiWrite(niiOrigIso);
niftiWrite(niiSetIso);
niftiWrite(niiSetRandIso);