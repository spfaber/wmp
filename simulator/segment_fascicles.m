
subDir = '/N/dc2/projects/lifebid/2t1/HCP/105115/fibers';
fibers = 'dwi_data_b2000_aligned_trilin_csd_lmax10_dwi_data_b2000_aligned_trilin_brainmask_dwi_data_b2000_aligned_trilin_wm_prob-500000.pdb';
SAVEDIR = '/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/simulator/110411/fascicles';

% point to dt6 file
dtFile = fullfile(subDir,'/dt6_b2000trilin','dt6.mat');
%wholeBrainConnectome = fgRead(fullfile(subDir, '/fibers',fibers));

feDir = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results';
feFile = fullfile(feDir, 'fe_structure_110411_ETC_500000.mat');
load(feFile);

% remove zero-weighted fibers
weights = feGet(fe,'fiber weights');
zeroInd = find(weights==0);
fe.fg.fibers(zeroInd)=[];
wholeBrainConnectome = feGet(fe, 'fibers acpc'); 

% Segment the fascicles using AFQ
[fascicles,classification,fg,fg_classified] = feAfqSegment(dtFile,wholeBrainConnectome);

save(fullfile(SAVEDIR, 'AFQ_segmented_fascicles.mat'),'fascicles','-v7.3'); %3,4,19,20 (left/right cortico, left/right arcuate)


for i = 1:20
    numFib(i) = length(fascicles(i).fibers);
end

fgWrite(fg_classified(3), fullfile(SAVEDIR,'fg_left_cort_500000_segmented.pdb'));
fgWrite(fg_classified(4), fullfile(SAVEDIR,'fg_right_cort_500000_segmented.pdb'));

fgWrite(fg_classified(19), fullfile(SAVEDIR,'fg_left_arc_500000_segmented.pdb'));
fgWrite(fg_classified(20), fullfile(SAVEDIR,'fg_right_arc_500000_segmented.pdb'));

fgDir = '/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/simulator/105115/fascicles';
fg = fgRead(fullfile(fgDir, 'fgrArcCortClean105115.mat'));


% try with retracking
[fg_classified,fg_unclassified,classification,fg] = AFQ_SegmentFiberGroups(dtFile);

% try with original fibers (just Arcuate and Cortico)
[fg_classified,fg_unclassified,classification,fg] = AFQ_SegmentFiberGroups(dtFile, fg);


fg_Orig_arc = fgRead(fullfile(fgDir,'fgrArcClean105115.mat'));

fg_Orig_cort = fgRead(fullfile(fgDir,'fgrCortClean105115.mat'));

fe = load('fe_rArCort_simulator_105115.mat');
weights = feGet(fe,'fiber weights');
weights_arc = weights(1:length(fg_Orig_arc.fibers));
weights_cort = weights(1+length(fg_Orig_arc.fibers):end);

fe_arc = fe;
fe_cort = fe;

fe_arc.fg.fibers = fe.fg.fibers(1:length(fg_Orig_arc.fibers));
fe_cort.fg.fibers = fe.fg.fibers(1+length(fg_Orig_arc.fibers):end);


fe_arc.fg.fibers = fe_arc.fg.fibers(weights_arc>0);
fe_cort.fg.fibers = fe_cort.fg.fibers(weights_cort>0);

%% plots

fg_arc_img = feGet(fe_arc, 'fg acpc' );
fg_cort_img = feGet(fe_cort, 'fg acpc');

figure (7)
for ii = 1:length(fg_arc_img.fibers);
    plot3(fg_arc_img.fibers{ii}(1,:),fg_arc_img.fibers{ii}(2,:),fg_arc_img.fibers{ii}(3,:),'b'); hold on
end


for ii = 1:length(fg_cort_img.fibers);
    plot3(fg_cort_img.fibers{ii}(1,:),fg_cort_img.fibers{ii}(2,:),fg_cort_img.fibers{ii}(3,:),'r'); hold on
end
view(90,0)

figure(8)
for ii = 1:length(fg_Orig_arc.fibers);
    plot3(fg_Orig_arc.fibers{ii}(1,:),fg_Orig_arc.fibers{ii}(2,:),fg_Orig_arc.fibers{ii}(3,:),'b'); hold on
end


for ii = 1:length(fg_Orig_cort.fibers);
    plot3(fg_Orig_cort.fibers{ii}(1,:),fg_Orig_cort.fibers{ii}(2,:),fg_Orig_cort.fibers{ii}(3,:),'r'); hold on
end
view(90,0)