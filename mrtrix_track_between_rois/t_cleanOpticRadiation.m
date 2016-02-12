% This is an exmample of how to load a previously computed optic radiation
% from an individual brain. I will show how to display the optic radiation
% and plot the profile of the optic radiation showing properties of the
% white matter.
% 10-27-2015
% Franco Pestilli and Sam Faber 

% set up the path to the data 
BASEDIR = '/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/mrtrix_track_between_rois/KK/'; %105115, 113619, 115320, 117122, 118730, HT, FP, JW, KK
TRKDIR = fullfile(BASEDIR, 'mrtrix_results');


% 1. We will load the fiber group
% This file was first read into MATLAB as a .pdb using fgRead. 
% It was then written out as .mat using fgWrite so that it could be 
% visualized and initially cleaned in mrDiffusion.
% It was then saved (in mrDiffusion) as an updated version and then read in
% again to be cleaned further.

mrtrix_tck2pdb(fullfile(TRKDIR,'whole_brain_ORs_excluded.tck'),fullfile(TRKDIR,'whole_brain_ORs_excluded.pdb'))
mrtrix_tck2pdb(fullfile(TRKDIR,'left_optic_radiation_PCSD.tck'),fullfile(TRKDIR,'left_optic_radiation_PCSD.pdb'))
mrtrix_tck2pdb(fullfile(TRKDIR,'right_optic_radiation_PCSD.tck'),fullfile(TRKDIR,'right_optic_radiation_PCSD.pdb'))

fg = fgRead(fullfile(TRKDIR, 'left_optic_radiation_PCSD.pdb'));
fgWrite(fg, fullfile(TRKDIR,fg.name),'mat'); % load this into mrDiffusion for initial cleaning

fg = fgRead(fullfile(TRKDIR, 'right_optic_radiation_PCSD.pdb'));
fgWrite(fg, fullfile(TRKDIR,fg.name),'mat'); % load this into mrDiffusion for initial cleaning

% Save initially cleaned file as _updated.mat and then read it in to be
% cleaned
fg1 = fgRead(fullfile(TRKDIR,'left_optic_radiation_PCSD_updated.mat'));
fg2 = fgRead(fullfile(TRKDIR,'right_optic_radiation_PCSD_updated.mat'));

% 2. We will clean the fiber group
tic, [fiberClean1, fibersToKeep] = mbaComputeFibersOutliers(fg1, 3, 3);
toc
    
tic, [fiberClean2, fibersToKeep] = mbaComputeFibersOutliers(fg2, 3, 3);
toc

%tic, [fiberClean3, fibersToKeep] = mbaComputeFibersOutliers(fiberClean2, 3, 3);
%toc


% % 3.a. We will measure the volume of the optic radiation.
% fibers = fiberClean1.fibers;
% for i = 1:size(fibers)
%     volOptRad = unique(transpose(fibers{i}),'rows');
% end
% length(volOptRad) % Number of unique voxels => volume of optic radiation
% % Unique voxels = 103 => volume = 103 mm^3
% 
% 
% % 3.b. We will measure the length of the optic radiation.
% for j = 1:size(fibers)
%     nodesPerFiber(j) = length(fibers{j});
% end

%mean(diff(nodesPerFiber))
%mean(nodesPerFiber) % average length of fiber  % 101.9091

%fiberLength = mean(fefgGet(fg,'length'));      % 98.1134


% 4. Plot a tract profile.
dt = dtiLoadDt6( fullfile(BASEDIR,'dt6.mat') );
[fa, md, rd, ad, cl, SuperFiber] = dtiComputeDiffusionPropertiesAlongFG( fiberClean1, dt,[],[],200);
 

% 5.0 plot the tract profile for FA and MD
h.tpfig = figure('name', 'OpticRadiation_TP','color', 'w');
plot(fa(50:151),'color', [0.2 0.2 0.9],'linewidth',4)
set(gca, 'fontsize',20, 'box','off', 'TickDir','out', ...
    'xticklabel',{'LGN','V1'},'xlim',[0 100],'ylim',[0.25 .75],'Ytick',[0 .25 .5 .75],'Xtick',[0 100])
title('Tract Profile')
xlabel('Location')
ylh = ylabel('Fractional Anisotropy');
 
feSavefig(h.tpfig,'verbose','yes', ...
        'figName','OpticRadiation_TP', ...
        'figDir',TRKDIR, ...
        'figType','jpg');

% 6. Plot the fiber group:
% 6.1 as a full fibers
h.fig = figure('name', 'OpticRadiation','color', 'k');
t1 = niftiRead(fullfile(BASEDIR, 't1.nii.gz'));
slices      = {[6 0 0],[0 1 0],[0 0 -15]}; 
hold on
h.fig  = mbaDisplayBrainSlice(t1, slices{1});
h.fig  = mbaDisplayBrainSlice(t1, slices{2});
h.fig  = mbaDisplayBrainSlice(t1, slices{3});
[h.fig, h.light] = mbaDisplayConnectome(fiberClean1.fibers, h.fig);
hold on 

% 6.2 as a super fiber with the FA map overlayed
%[h.fig, h.light] = mbaDisplayConnectome(SuperFiber.fibers,h.fig,[0 .25 .75],'single',[], [], 2);
axis([0 90 -110 0 -90 90])
view(0,90)
camlight('left')
 
feSavefig(h.fig,'verbose','yes', ...
        'figName','OpticRadiation', ...
        'figDir',TRKDIR, ...
        'figType','jpg');


%cMap = colormap(jet);
%zMap = zeros(101,2);
%faMap = [zMap fa(50:150)];
 % map of dimentions as requsted by mbaDisplayConnectome with values as returned inside FA
%[h.fig, h.light] = mbaDisplayConnectome(SuperFiber.fibers,h.fig, [0.9 0.2 0.1],'map',cMap, [], 1.5);
%colorbar

% Change light and view (orientation to optimize) 
% Print the image to jpg using fePrintFigure.m print with 2-3 orientations
% to choose from.
