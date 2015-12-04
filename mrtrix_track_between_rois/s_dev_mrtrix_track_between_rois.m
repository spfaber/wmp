function fibersPDB = s_dev_mrtrix_track_between_rois
%
% This functions shows how to track between two ROIS using mrtrix.
% This is very helpful for ideintifying some fiber groups, for example the
% optic radiation.
%
% This is how the code works.
% 1. We load two ROIs in the brain, for Shumpei's project for example we
%    will load the right-LGN and the right-Visual cortex
% 2. We create union ROI by combining these two ROIs. The union ROI is used
%    as seeding for the fibers. mrtrix will initiate and terminate fibers only
%    within the volume defined by the Union ROI.
% 3. We create a white matter mask. THis mask is generally a large portion
%    of the white matter. A portion that contains both union ROIs. For example
%    the right hemisphere.
% 4. We use mrtrix to track between the right-LGN and righ-visual cortex.
% mrtrix will initiate fibers by seeding within the UNION ROI and it will
% only keep fibers that have paths within the white matter masks.
%
% The final result of this script is to generate lot's of candidate fibers 
% that specifically end and start from the ROI of interest. Thisis an 
% approach similar to Contrack. 
%
% INPUTS: none
% OUTPUTS: the finela name of the ROI created at each iteration
%
% Written by Franco Pestilli (c) Stanford University Vistasoft
% edited by Sam Faber 10/12/15
%LD_LIBRARY_PATH = '/N/soft/rhel6/gcc/4.9.2/lib64:/N/soft/rhel6/gsl/1.15/lib/:LD_LIBRARY_PATH'; 

%setenv('LD_LIBRARY_PATH', LD_LIBRARY_PATH);

old_lib_path = getenv('LD_LIBRARY_PATH');
new_lib_path = [ '/N/soft/rhel6/gcc/4.9.2/lib64:/N/soft/rhel6/gsl/1.15/lib/:' old_lib_path];
setenv('LD_LIBRARY_PATH', new_lib_path)

%Now the LD_LIBRARY_PATH is updated
!echo $LD_LIBRARY_PATH

%Get back old LD_LIBRARY_PATH 
setenv('LD_LIBRARY_PATH', old_lib_path)
%You should expect to see the old LD_LIBRARY_PATH 
!echo $LD_LIBRARY_PATH


baseDir = '/N/dc2/projects/lifebid/HCP/Sam';
codeDir = '/matlab_code/mrtrix_track_between_rois/';


dtFile = fullfile(baseDir, '/105115/diffusion_data/dt6_b2000trilin/dt6.mat');
refImg = fullfile(baseDir, '/105115/anatomy/T1w_acpc_dc_restore_1p25.nii.gz');
fibersFolder = fullfile(baseDir, codeDir, 'fibers');

% We want to track the optic radiation (LGN -> V1)
fromRois = load(fullfile(baseDir, codeDir, 'ROIs/rh_LGN.mat'));
toRois   = load(fullfile(baseDir, codeDir, 'ROIs/rh_V1_thresh_label_smooth3mm.mat'));

% Set up the MRtrix traking parameters
trackingAlgorithm = {'prob'};
lmax    = [10]; % The appropriate value depends on # of directions. For 32, use lower #'s like 4 or 6. For, 6 or 10 is good [10];
nSeeds  = 500; % 10000; 
nFibers = 500; %1000000;
wmMask  = [];

% Make an (include) white matter mask ROI. This mask is the smallest
% set of white matter that contains both ROIS (fromRois and toRois)
%
% We use a nifti ROi to select the portion of the White matter to use for
% seeding
wmMaskName = fullfile(baseDir, '/105115/anatomy/wm_mask_Sam'); 


% Then transform the niftis into .mif
[p,f,e] = fileparts(wmMaskName);
wmMaskMifName    = fullfile(p,sprintf('%s.mif',f)); 
wmMaskNiftiName  = sprintf('%s.nii.gz',wmMaskName);
mrtrix_mrconvert(wmMaskNiftiName, wmMaskMifName);    

% This first step initializes all the files necessary for mrtrix.
% This can take a long time.
files = mrtrix_init(dtFile,lmax,fibersFolder,wmMaskNiftiName);

% Some of the following steps only need to be done once for each ROI,
% so we want to do some sort of unique operation on the from/toRois
individualRois = unique([fromRois, toRois]);

% Convert the ROIs from .mat or .nii.gz to .mif format.
for i_roi = 1:length(individualRois)
    
    if   exist(fullfile(p, [individualRois{i_roi}, '.mat']),'file')
         thisroi = fullfile(p, [individualRois{i_roi}, '.mat']);
    end
    
    mrtrix_roi2mif(thisroi,refImg);
end
    
% Create joint from/to Rois to use as a mask
for nRoi = 1:length(fromRois)
    % MRTRIX tracking between 2 ROIs template.
    roi{1} = fullfile(baseDir, codeDir, '/ROIs/', fromRois{nRoi});
    roi{2} = fullfile(baseDir, codeDir, '/ROIs/', toRois{nRoi});
    
    roi1 = dtiRoiFromNifti([roi{1} '.nii.gz'],[],[],'.mat');
    roi2 = dtiRoiFromNifti([roi{2} '.nii.gz'],[],[],'.mat');
    
    % Make a union ROI to use as a seed mask:
    % We will generate as many seeds as requested but only inside the volume
    % defined by the Union ROI.
    %
    % The union ROI is used as seed, fibers will be generated starting ONLY
    % within this union ROI.
    roiUnion        = roi1; % seed union roi with roi1 info
    roiUnion.name   = ['union of ' roi1.name ' and ' roi2.name]; % r LGN calcarine';
    roiUnion.coords = vertcat(roiUnion.coords,roi2.coords);
    roiName         = fullfile(baseDir, codeDir, '/ROIs/',[roi1.name '_' roi2.name '_union']);
    [~, seedMask]   = dtiRoiNiftiFromMat(roiUnion,refImg,roiName,1);
    seedRoiNiftiName= sprintf('%s.nii.gz',seedMask);
    seedRoiMifName  = sprintf('%s.mif',seedMask); 
    
    % Transform the niftis into .mif
    mrtrix_mrconvert(seedRoiNiftiName, seedRoiMifName);
        
    % We cd into the folder where we want to save the fibers.
    cd(fibersFolder);
    
    % We generate and save the fibers in the current folder.
    [fibersPDB{nRoi}, status, results] = mrtrix_track_roi2roi(files, [roi{1} '.mif'], [roi{2} '.mif'], ...
        seedRoiMifName, wmMaskMifName, trackingAlgorithm{1}, ...
        nSeeds, nFibers);
    
    % fgWrite(fibersPDB,['fibername'],'pwd')
end

return
