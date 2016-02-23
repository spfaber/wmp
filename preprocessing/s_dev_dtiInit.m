% This script is an example of how to preprocess data using dtiInit. 

% Decide conditions that we are going to run (we will run different
% parameters set with dtiInit to test how parameters affect preprocessing).
% dtiInit_folderNames = {'plusBvec_Rx0_Can0', 'plusBvec_Rx1_Can0', 'plusBvec_Rx0_Can1', ...
%    'plusBvec_Rx1_Can1','minusBvec_Rx0_Can0', 'minusBvec_Rx1_Can0','minusBvec_Rx0_Can1', 'minusBvec_Rx1_Can1'};

dwi_dir  = '/N/dc2/projects/lifebid/HCP/Sam/dtiInit_test/takemura_fix/raw_diffusion/';
save_dir = '/N/dc2/projects/lifebid/HCP/Sam/dtiInit_test/preprocessing_test';
anat_dir = '/N/dc2/projects/lifebid/HCP/Sam/dtiInit_test/takemura_fix/anatomy/';


% Initialize the dtiInit parameters
dwParams = dtiInitParams;
dwParams.clobber = 0;
%dwParams.phaseEncodeDir = 2;


% make bvals vector from text file
bvals_filename = fullfile(dwi_dir, 'bvals_takemura.txt');
bvals_takemura = dlmread(bvals_filename, '%s\n');
bvals_takemura( (bvals_takemura == 10) )= 0;
dlmwrite(fullfile(dwi_dir,'dti_2mm_b1000_2000_ap_2_reform.bval'),bvals_takemura);

% for i = 1:length(dtiInit_folderNames)
    % reformat bvecs & bvals from .txt files
    dwParams.outDir      = fullfile(save_dir, 'plusBvec_Rx0_Can1');%dtiInit_folderNames{i}) ;
    dwParams.dt6BaseName = fullfile(save_dir, 'plusBvec_Rx0_Can1') ; 
    
    % make bvecs matrix from text file
    bvecs_filename = fullfile(dwi_dir, 'bvecs_takemura.txt');
    bvecs_takemura = dlmread(bvecs_filename);
    bvecs_takemura( 1, (bvals_takemura == 0) ) = 0;
    bvecs_takemura( 2, (bvals_takemura == 0) ) = 0;
    bvecs_takemura( 3, (bvals_takemura == 0) ) = 0;
    
    % Depending on the condition we will change the BVECS file or NOT
    %switch ( dtiInit_folderNames{i})
     %   case {'plusBvec_Rx0_Can0'}
       %     dwParams.rotateBvecsWithCanXform = 0;
       %     dwParams.rotateBvecsWithRx = 0;
            
%        case {'plusBvec_Rx1_Can0'}
%             dwParams.rotateBvecsWithCanXform = 0;
%             dwParams.rotateBvecsWithRx = 1;
%             
%        case {'plusBvec_Rx0_Can1'}
            dwParams.rotateBvecsWithCanXform = 1;
            dwParams.rotateBvecsWithRx = 0;
            
%         case {'plusBvec_Rx1_Can1'}
%             dwParams.rotateBvecsWithCanXform = 1;
%             dwParams.rotateBvecsWithRx = 1;
%             
%           
%         case {'minusBvec_Rx0_Can0'}
%             dwParams.rotateBvecsWithCanXform = 0;
%             dwParams.rotateBvecsWithRx = 0;
%             bvecs_takemura(2,:) = -bvecs_takemura(2,:);
%             
%         case {'minusBvec_Rx1_Can0'}
%             dwParams.rotateBvecsWithCanXform = 0;
%             dwParams.rotateBvecsWithRx = 1;
%             bvecs_takemura(2,:) = -bvecs_takemura(2,:);
%             
%         case {'minusBvec_Rx0_Can1'}
%             dwParams.rotateBvecsWithCanXform = 1;
%             dwParams.rotateBvecsWithRx = 0;
%             bvecs_takemura(2,:) = -bvecs_takemura(2,:);
%             
%         case {'minusBvec_Rx1_Can1'}
%             dwParams.rotateBvecsWithCanXform = 1;
%             dwParams.rotateBvecsWithRx = 1;
%             bvecs_takemura(2,:) = -bvecs_takemura(2,:);
%            
%         %otherwise
%          %   keyboard
%             
%     end

    % We always write the file to disk by adding the TAG from the current
    % condition we are testing.
    bvecs_fname = sprintf('dti_2mm_b1000_2000_ap_2_reform_%s.bvec','plusBvec_Rx0_Can1');
    dlmwrite(fullfile(dwi_dir,bvecs_fname),bvecs_takemura)
    
    % Assign reformatted bvecs & bvals
    dwParams.bvecsFile = fullfile(dwi_dir, bvecs_fname);
    dwParams.bvalsFile = fullfile(dwi_dir, 'dti_2mm_b1000_2000_ap_2_reform.bval');
    dwifile   = fullfile(dwi_dir,'dwi_b1000_2000_ap_2.nii.gz');
    anatfile  = fullfile(anat_dir,'t1.nii.gz');
     
    % Do the actual preprocessing given the current parameters
    [dt6FileName, outBaseDir] = dtiInit(dwifile, anatfile, dwParams);
% end
