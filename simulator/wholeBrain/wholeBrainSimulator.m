% 20 major tracts simulator code
% Sam Faber
% 20160412


%% Set up paths and build the file names for the diffusion data, bvecs, bvals, and fiber group. 
subj      = '105115';
rootDir   = fullfile('/N/dc2/projects/lifebid/HCP/Sam/',subj);
saveDir   = fullfile('/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/simulator/wholeBrain');
dwiFile   = fullfile(rootDir,'diffusion_data','dwi_data_b2000_aligned_trilin.nii.gz');
t1File    = fullfile(rootDir,'anatomy','T1w_acpc_dc_restore_1p25.nii.gz');
dtFile    = fullfile(rootDir,'diffusion_data','/dt6_b2000trilin','dt6.mat');
bvecsFile = fullfile(rootDir,'diffusion_data','dwi_data_b2000_aligned_trilin.bvecs');
bvalsFile = fullfile(rootDir,'diffusion_data','dwi_data_b2000_aligned_trilin.bvals');
bvecs     = dlmread(bvecsFile);
bvals     = dlmread(bvalsFile);


%% fe struct for 100,000 fibers before segmentation

fe = load(fullfile(saveDir,'fe',strcat('fe_b2000_aligned_trilin_100000_csd_',subj,'.mat')));

weights = feGet(fe,'fiber weights');
zeroInd = find(weights==0);
fe.fg.fibers(zeroInd)=[];
wholeBrainConnectome = feGet(fe, 'fibers acpc'); 


%% Segment connectome into 20 major fiber groups
fg_orig = fgRead(fullfile(saveDir,'mrtrix_results','orig_20maj_masked.pdb'));
fg_weighted = fgRead(fullfile(saveDir,'mrtrix_results','weighted_20maj_masked.pdb'));

% original segmentation of tracked fibers
% once segmented, unclassified fibers cannot be segmented
[fg_classified_or, fg_unclassified_or, classification_or, fg_or] = AFQ_SegmentFiberGroups(dtFile,fg_orig);
[fg_classified_wt, fg_unclassified_wt, classification_wt, fg_wt] = AFQ_SegmentFiberGroups(dtFile,fg_weighted);


% check number of total fibers that have been classified
for i = 1:20
    numClassFib(i) = length(fg_classified_nz(i).fibers);
end
totClassFib = sum(numClassFib)


for i = 1:20
    numUnClassFib_or(i) = length(fg_classified_or(i).fibers);
end
sum(numUnClassFib_or)

% load 20 maj fibers
load(fullfile(saveDir,'segmentation','AFQ_wb_nz_segmented_classified_fg.mat'));

fgMerge1 = fgMerge(fg_classified_nz(1),fg_classified_nz(2),'merge1');
fgMerge2 = fgMerge(fgMerge1,fg_classified_nz(3),'merge1');
fgMerge3 = fgMerge(fgMerge2,fg_classified_nz(4),'merge1');
fgMerge4 = fgMerge(fgMerge3,fg_classified_nz(5),'merge1');
fgMerge5 = fgMerge(fgMerge4,fg_classified_nz(6),'merge1');
fgMerge6 = fgMerge(fgMerge5,fg_classified_nz(7),'merge1');
fgMerge7 = fgMerge(fgMerge6,fg_classified_nz(8),'merge1');
fgMerge8 = fgMerge(fgMerge7,fg_classified_nz(9),'merge1');
fgMerge9 = fgMerge(fgMerge8,fg_classified_nz(10),'merge1');
fgMerge10 = fgMerge(fgMerge9,fg_classified_nz(11),'merge1');
fgMerge11 = fgMerge(fgMerge10,fg_classified_nz(12),'merge1');
fgMerge12 = fgMerge(fgMerge11,fg_classified_nz(13),'merge1');
fgMerge13 = fgMerge(fgMerge12,fg_classified_nz(14),'merge1');
fgMerge14 = fgMerge(fgMerge13,fg_classified_nz(15),'merge1');
fgMerge15 = fgMerge(fgMerge14,fg_classified_nz(16),'merge1');
fgMerge16 = fgMerge15;
fgMerge16.fibers = vertcat(fgMerge15.fibers,fg_classified_nz(17).fibers);
fgMerge17 = fgMerge(fgMerge16,fg_classified_nz(18),'merge1');
fgMerge18 = fgMerge(fgMerge17,fg_classified_nz(19),'merge1');
fg20maj   = fgMerge(fgMerge18,fg_classified_nz(20));
fg20maj.name = '20_major_fibergroups_merged';

% save 20 major fiber groups
save(fullfile(saveDir,'fascicles','fg_nz_20major.mat'),'fg20maj');
% plot fibers

c = {'r','r','m','m','b','b','g','g','y','c',[1 0.6 0], [1 0.6 0], [0 0.5 0],...
    [0 0.5 0],[0.9 0.75 0], [0.9 0.75 0],[0.5 0.5 1] , [0.5 0.5 1], ...
    [0.3 0.8 0],[0.3 0.8 0]};


figure('name','Classified Weighted Fibers','color','w')
for i = 1:20;
for ii = 1:length(fg_classified_nz(i).fibers);
    plot3(fg_classified_nz(i).fibers{ii}(1,:),fg_classified_nz(i).fibers{ii}(2,:),fg_classified_nz(i).fibers{ii}(3,:),'color',c{i}); hold on
end
end
hold on 
%legend('Left TR','Right TR','Left CST','Right CST','Left CC','Right CC',...
%    'Left CH','Right CH','CFMaj','CFMin','Left IFOF','Right IFOF','Right ILF', ...
%    'Left ILF','Left SLF','Right SLF','Left Uncinate','Right Uncinate','Left AF','Right AF');
axes off


figure('name','Classified Original Fibers','color','w')
for i = 1:20;
for ii = 1:length(fg_classified_or(i).fibers);
    plot3(fg_classified_or(i).fibers{ii}(1,:),fg_classified_or(i).fibers{ii}(2,:),fg_classified_or(i).fibers{ii}(3,:)); hold on
end
end


figure('name','UnClassified Original Fibers','color','w')
for ii = 1:100 % length(fg_unclassified_or.fibers);
    plot3(fg_unclassified_or.fibers{ii}(1,:),fg_unclassified_or.fibers{ii}(2,:),fg_unclassified_or.fibers{ii}(3,:)); hold on
end

figure('name','UnClassified Weighted Fibers','color','w')
for ii = 1:length(fg_unclassified_wt.fibers);
    plot3(fg_unclassified_wt.fibers{ii}(1,:),fg_unclassified_wt.fibers{ii}(2,:),fg_unclassified_wt.fibers{ii}(3,:)); hold on
end




% for jj = 1:length(fg_unclassified_or.fibers{ii});
%     x = fg_unclassified_or.fibers{ii}(1,:);
%     y = fg_unclassified_or.fibers{ii}(2,:);
%     z = fg_unclassified_or.fibers{ii}(3,:);
% d_withinfib(jj) = sqrt((x(jj+1)-x(jj+2))^2 + (y(jj+1)-y(jj+2))^2 + (z(jj+1) - z(jj+2))^2);
% 
% end
%%
% for ii = 1:100;
%     numCols = length(fg_unclassified_or.fibers{ii});
%     T{ii} = nchoosek(1:numCols,2); % pairwise indexes for combinations of columns
%     
%     for k = 1:length(T{ii});
%          d(k) = norm(fg_unclassified_or.fibers{ii}(:,T{ii}(k,1))-fg_unclassified_or.fibers{ii}(:,T{ii}(k,2)));
%     end
%     [~ ,minIndex] = min(d);
%     clear d
%     T{ii}(minIndex,:);
%     centroids{ii} = (fg_unclassified_or.fibers{ii}(:,T{ii}(minIndex,1))+...
%         fg_unclassified_or.fibers{ii}(:,T{ii}(minIndex,2)))/2;
% end

% get centroids for each fiber
for ii = 1:length(fg_unclassified_or.fibers);
numCols = length(fg_unclassified_or.fibers{ii});
if mod(numCols,2)==0
    centPts = [numCols/2 numCols/2+1];
    centroids{ii} = (fg_unclassified_or.fibers{ii}(:,centPts(:,1))+ ...
        fg_unclassified_or.fibers{ii}(:,centPts(:,2)))/2;
elseif mod(numCols,2)==1
    centInd = numCols/2 + 1/2;
    centroids{ii} = fg_unclassified_or.fibers{ii}(:,centInd);
end 
end

clear d k ii T numCols mins minIndex centPts centInd icent


% test figure
figure('name','UnClassified Original Fibers','color','w')
for ii = 1:10:10000 % length(fg_unclassified_or.fibers);
    plot3(fg_unclassified_or.fibers{ii}(1,:),fg_unclassified_or.fibers{ii}(2,:),fg_unclassified_or.fibers{ii}(3,:)); hold on
end
hold on
for icent = 1:10:10000 %1:length(centroids);
 scatter3(centroids{icent}(1,1),centroids{icent}(2,1),centroids{icent}(3,1),'k'); hold on
end
view(0,90)


% test figure (15 mm threshold)
figure
for i = [1 ,2, 3, 4, 5, 7, 8, 9, 10];
plot3(fg_unclassified_or.fibers{i}(1,:),fg_unclassified_or.fibers{i}(2,:),fg_unclassified_or.fibers{i}(3,:)); hold on
scatter3(centroids{i}(1,1),centroids{i}(2,1),centroids{i}(3,1)); 
end

%%
% do pairwise distance on centroids (centroid index is fiber index)
numFib =100;
for ii = 1:numFib;
   pairInd = nchoosek(1:numFib,2); % pairwise indexes for combinations of fibers
    for k = 1:length(pairInd);
         d(k) = norm(centroids{pairInd(k,1)}-centroids{pairInd(k,2)});
    end
end
% Find pairwise indexes of fiber distances -create groups
% based on distances and angles
[~ ,minIndex] = find(d);
pairInd(minIndex,:);


% for each centroid, look at all distances from that centroid to other centroids,
% if distances are <10mm and <20 degrees away, keep them, if not, do not
% keep them and move on to the next centroid
d = [];
for i = 1:10; %length(centroids)
    for j = 1:10 ;%length(centroids);
    d(i,j) = norm(centroids{i} - centroids{j}); % matrix of distances between centroids
    end
    
end
centMat = cell2mat(centroids);
centMat = centMat';

k = 10;
[idxbest, Cbest, sumDbest, Dbest] = kmeans(centMat, k,'Distance','cosine');
gcf 
hold on
Cbest = Cbest';

figure

for i = 1:length(Cbest);
scatter3(Cbest(i,1),Cbest(i,2),Cbest(i,3),15,'k*'); hold on
end

numClust = [40 30 20 10];
avgDist2Clust = [0.8821 0.8707 0.8655 0.8531];

% convert to non-repeating, ordered column-vector of fiber indices

%fibers_5mm = fg_unclassified_or.fibers{pairInd(minIndex)};


clear d k ii numFib pairInd

% group centroids/fibers based on distance threshold and angle between
% centroids (get unit vector at each of the closest points (divide points by magnitude)
% and subtract these to get unit vector in direction of fiber at centroid.
% Do this for both centroids and then compute their dot product and do
% inverse cos to get angle (convert to degrees). Threshold at 10 mm and 20 degrees to
% start out with.

% get vector in direction of fiber at centroid
for ii = 1:10; %length(fg_unclassified_or.fibers);
numCols = length(fg_unclassified_or.fibers{ii});
if mod(numCols,2)==0
    centPts = [numCols/2 numCols/2+1];
    centUnits{ii} = [(fg_unclassified_or.fibers{ii}(:,centPts(:,1)))/norm(fg_unclassified_or.fibers{ii}(:,centPts(:,1)))...
        (fg_unclassified_or.fibers{ii}(:,centPts(:,2)))/norm(fg_unclassified_or.fibers{ii}(:,centPts(:,2))) ];
   fibDir{ii} = centUnits{ii}(:,1)-centUnits{ii}(:,2);
elseif mod(numCols,2)==1
    centInds = [numCols/2 + 1/2 numCols/2+1.5];
    centUnits{ii} = [(fg_unclassified_or.fibers{ii}(:,centInds(:,1)))/norm(fg_unclassified_or.fibers{ii}(:,centInds(:,1)))...
        (fg_unclassified_or.fibers{ii}(:,centInds(:,2)))/norm(fg_unclassified_or.fibers{ii}(:,centInds(:,2))) ];
   fibDir{ii} = centUnits{ii}(:,1)-centUnits{ii}(:,2); 
end 
end

% look at vectors
fibDirMat = cell2mat(fibDir);
gcf 
hold on
plotv(fibDirMat, '-')

figure
for i = 1:length(fibDir);
quiver3(centroids{i}(1,1),centroids{i}(2,1),centroids{i}(3,1),fibDir{i}(1,1),...
    fibDir{i}(2,1),fibDir{i}(3,1),0);
end

% look into literature about angles between fibers
angles = [];
for i = 1:10; %length(centroids)
    for j = 1:10 ;%length(centroids);
    angles(i,j) = radtodeg(acos(dot(fibDir{i},fibDir{j}))); % matrix of distances between centroids
    end 
end

clear numCols fibDir centPts centUnits centInd ii




save(fullfile(saveDir,'segmentation', 'AFQ_wb_segmented_classified_or_fg.mat'),'fg_classified_or','-v7.3');
save(fullfile(saveDir, 'segmentation', 'AFQ_wb_segmented_unclassified_or_fg.mat'),'fg_unclassified_or','-v7.3');
save(fullfile(saveDir,'segmentation', 'AFQ_wb_segmented_classified_wt_fg.mat'),'fg_classified_wt','-v7.3');
save(fullfile(saveDir, 'segmentation', 'AFQ_wb_segmented_unclassified_wt_fg.mat'),'fg_unclassified_wt','-v7.3');


%% Run LIFE on segmented fibers to get weights and signal in voxels
fgFileName    = fg20maj;
feFileName    = strcat('fe_nz_20maj_b2000_aligned_trilin_100000_csd_',subj);

% run LiFE
disp('Building life model...');
N = 360; % Discretization parameter
fe = feConnectomeInit(dwiFile,fgFileName,feFileName,saveDir,[],t1File,N,[1,0],0);
disp('done');

% fit the model
disp('Fitting life model...');
tic
fe = feSet(fe,'fit',feFitModel(feGet(fe,'model'),feGet(fe,'dsigdemeaned'),'bbnnls'));
toc
disp('done');

% save LiFE fe struct
save(fullfile(saveDir,'fe',strcat(feFileName,'.mat')),'-struct','fe','-v7.3'); 

%% get original signal in voxels of 20 major fiber groups and save out
% make wm_mask of signal
% then track on new signal

% read in nifti file
nii = niftiRead(dwiFile);

% rename diffusion files to keep track of different signals
niiOrig    = nii;
niiNew     = nii;
niiWM      = nii;

% get coords_acpc of nodes in voxels
coords = fe.roi.coords;  
dw_vals = feGet(fe,'dsiinvox',coords);

% find indexes of nonzero and zero b-values
indexes = find(bvals~=0);
b0indexes = find(bvals==0);

% make b0_data a matrix (# b0 X # voxels) 
b0_data = nan(size(b0indexes,2),size(coords,1));
for ivx = 1:size(coords,1)
    b0_data(:,ivx) = niiOrig.data(coords(ivx,1),coords(ivx,2),coords(ivx,3),b0indexes);
end

% set signal to zero everywhere outside AF and CST
niiOrig.data    = zeros(size(niiOrig.data));
niiWM.data      = zeros(size(niiWM.data));

% replace image values to get original signal only in AF and CST
niiOrig.data    = feReplaceImageValues(niiOrig.data,b0_data,coords,b0indexes);
niiOrig.data    = feReplaceImageValues(niiOrig.data,dw_vals,coords,indexes);
niiOrig.data    = int16(niiOrig.data);
b0_valsMask = ones(size(b0_data));
dw_valsMask = ones(size(dw_vals));
niiWM.data = feReplaceImageValues(niiWM.data, dw_valsMask, coords, indexes);
niiWM.data = feReplaceImageValues(niiWM.data, b0_valsMask, coords, b0indexes);
niiWM.data      = int16(niiWM.data);

niiOrig.fname  = fullfile(saveDir, 'diffusion',...
    sprintf('orig_sig_20maj_%s.nii.gz',subj));
niiWM.fname    = fullfile(saveDir, 'diffusion',...
    sprintf('wm_mask_%s.nii.gz',subj));

% write out original signal only in AF and CST
niftiWrite(niiOrig);
niftiWrite(niiWM);

%% get weighed/predicted signal in voxels of 20 major fiber groups and save out
% then track on new signal
% also plot weight distributions of all fiber groups

predSig = feGet(fe,'pSig fiber');
iso = (fe.life.diffusion_signal_img(:) -feGet(fe,'dsigdemeaned'));
fullPredSig = pSig + iso;

fullPredSig = reshape(fullPredSig, size(coords,1),length(indexes));

% Replace nans with b0_data and fullPredSig
niiNew.data = zeros(size(niiNew.data));
niiNew.data = feReplaceImageValues(niiNew.data,b0_data,coords,b0indexes);
niiNew.data = feReplaceImageValues(niiNew.data,fullPredSig',coords,indexes);

% set correct nifti data type
niiNew.data = int16(niiNew.data);

% set name of new nifti and write it out
niiNew.fname = fullfile(saveDir,'diffusion', ...
    sprintf('weighted_sig_20maj_%s.nii.gz',subj));
niftiWrite(niiNew);

