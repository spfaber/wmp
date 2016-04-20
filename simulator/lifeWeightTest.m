% LiFE weight test
% Can we recover LiFE weights after running LiFE a second time ? 

% set up paths to data
subj      = '105115';
rootDir   = fullfile('/N/dc2/projects/lifebid/HCP/Sam/',subj);
saveDir   = fullfile('/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/simulator/',subj);
dwiFile   = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.nii.gz');
t1File    = fullfile(rootDir,'anatomy','T1w_acpc_dc_restore_1p25.nii.gz');
bvecsFile = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.bvecs');
bvalsFile = fullfile(rootDir,'dwi_data_b2000_aligned_trilin.bvals');
bvecs     = dlmread(bvecsFile);
bvals     = dlmread(bvalsFile);

% load fe struct and get weights
feOrig = load(fullfile(saveDir,'fe',strcat('fe_rArCort_simulator_105115','.mat')));
weights_orig = feGet(feOrig,'fiber weights');

% get LiFE predicted signal 
pSig    = feGet(feOrig, 'pSig fiber');
iso = (feOrig.life.diffusion_signal_img(:) -feGet(feOrig,'dsigdemeaned'));
fullPredSig = pSig + iso;

coords    = feOrig.roi.coords;  % get coords of nodes in voxels
indexes   = find(bvals~=0); % get diffusion weighted indexes
b0indexes = find(bvals==0); % get b0 indexes

fullPredSig = reshape(fullPredSig, size(coords,1),length(indexes));

nii = niftiRead(dwiFile); % read in nifti
niiOut = nii;             % rename 

% get b0 data as matrix of # voxels by # b0's 
b0_data = nan(size(b0indexes,2),size(coords,1));
for ivx = 1:size(coords,1)
    b0_data(:,ivx) = niiOut.data(coords(ivx,1),coords(ivx,2),coords(ivx,3),b0indexes);
end

niiOut.data = nan(size(niiOut.data));

% Replace nans with b0_data and fullPredSig
niiOut.data = feReplaceImageValues(niiOut.data,b0_data,coords,b0indexes);
niiOut.data = feReplaceImageValues(niiOut.data,fullPredSig',coords,indexes);

% set correct nifti data type
niiOut.data = int16(niiOut.data);

% chane nifti name
niiOut.fname = fullfile(saveDir,'diffusion','ArcCort',strcat('lifeSig_lifeWeightTest5_',subj));
dlmwrite([niiOut.fname,'.bvecs'],bvecs);
dlmwrite([niiOut.fname,'.bvals'],bvals);


% write out nifti
niftiWrite(niiOut)

% read in newly created dwi file
dwiFile = fullfile(saveDir, 'diffusion','ArcCort',strcat('lifeSig_lifeWeightTest5_',subj,'.nii.gz'));

% read in original fg 
fgArcCortClean = fgRead(fullfile(saveDir,'fascicles', strcat('fgrArcCortClean',subj,'.mat')));

% set LiFE params
fgFileName    = fgArcCortClean;
feFileName    = strcat('fe_rArcCort_sim_lifeWeightTest5_',subj);
savedir = saveDir;


% build LiFE model
disp('Building life model...');
N = 360; % Discretization parameter
fe = feConnectomeInit(dwiFile,fgFileName,feFileName,savedir,[],t1File,N,[1,0],0);
disp('done');

% fit LiFE model
disp('Fitting life model...');
tic
fe = feSet(fe,'fit',feFitModel(feGet(fe,'model'),feGet(fe,'dsigdemeaned'),'bbnnls'));
toc
disp('done');

% save fe struct
save(fullfile(saveDir,'fe',strcat(feFileName,'.mat')),'-struct','fe','-v7.3');


% get new LiFE weights
weights_new = feGet(fe,'fiber weights');

% load fe struct and get weights
feOrig = load(fullfile(saveDir,'fe',strcat('fe_rArCort_simulator_105115','.mat')));
weights_orig = feGet(feOrig,'fiber weights');


% We load several FE each one was a re-fit of the LiFE model on the signal
% (Y) predicted by an original LiFE model
weights = nan(length(weights_orig),5);
for iiLif = 1:5
    if iiLif == 1
        fe = load('fe_rArcCort_sim_lifeWeightTest_105115.mat');
    else
        fe = load(sprintf('fe_rArcCort_sim_lifeWeightTest%i_105115.mat',iiLif));
    end
    
    weights(:,iiLif) = feGet(fe,'fiber weights');
end



% plot weight distributions
%weight_bins = log10(logspace(min(new_weights),max(new_weights),30));
bins = linspace(-6, 0, 60);
figure('name','LiFE Weights Test','color','w'); 
hold on
[y, x] = hist(log10(weights_orig(weights_orig>0)), bins);
plot(x,y,'-','linewidth',2,'color',[.8 .6 .4])
hold on
for iiLif = 1:size(weights_orig,2)
    w_temp = weights(:,iiLif);
    [y, x] = hist(log10(w_temp(w_temp > 0)), bins);   
    plot(x, y,'o',  ...
        'markerfacecolor',[.3 .5 .9],  ...
        'markeredgecolor',[.3 .5 .9]);
    hold on
end
hold on

legend({'Dist of weights on Y Orig';'Dist of weights on $\hat{Y}$ Pred'})
set(legend,'Interpreter','latex','fontsize',12)
title('LiFE Weights Test')
ylabel('Number of Fascicles','fontsize',15)
xlabel('log_{10}(Fiber Weight)','fontsize',15)
set(gca,'xlim',[-6 0], 'xtick',[-6 -3 0],'ylim',[ 0 100], 'ytick',...
    [0 50 100], 'tickdir', 'out', 'box', 'off', 'fontsize', 15, 'visible', 'on')

