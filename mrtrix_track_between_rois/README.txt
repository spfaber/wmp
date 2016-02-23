This folder contains code and data for using mrtrix to track between ROIs. 
There are many steps in this process. They are listed and described below.

1.a. Create ROIs which you wish to track between using mrDiffusion (.mat) or freesurfer labels (.nii)
        - If freesurfer labels have not already been made into nifti ROIs, see example code: s_fe_make_all_freesurfer_labels_into_rois.m. This code can be found in the life_scripts repository. 
                - To run this code, set environment in terminal before opening matlab
                        - type 'export SUBJECTS_DIR=/path/to/subjects/directory/' in command line 
                        - make sure fsl and freesurfer modules are loaded (module load fsl, module load freesurfer)
        - Save them in a known location for later use

  b. Merge two ROIs in mrDiffusion by loading them (menu: File > ROIs> Load ROIs) and then merging (menu: ROIs > Merge)    
        - we will need this merged ROI to be our seed region for tracking  

  c. If you obtain .mat ROIs from mrDiffusion, use vistasoft function 'dtiRoiNiftiFromMat.m' or 'dtiExportRoiToNifti.m' to         create .nii version
        - an example of this is in 'ConvertMatRoiToNifti.m'
        
2. Transform the .bvecs and .bvals into a .b file for mrtrix using function 'create_b_file.R' 
        - Save in a known location for later use
        - to run R script make sure 'create_b_file.R' is in folder with .bvecs and .bvals you want to transform
                - type 'Rscript create_b_file.R' in command line
        
3. Create a white matter mask using 's_dev_fe_make_wm_mask_hcp.m'
        - Save in a known location for later use

4. Run mrtrix bash shell script making sure to save your track file (.tck) to a known location
        - to run shell script type './filename.sh' in command line

5. Convert .tck file to .pdb using 'mrtrix_tck2pdb('path/to/.tck','path/to/.pdb')' in MATLAB
        - example in 't_cleanOpticRadiation.m'

6. Before Cleaning, write .pdb to .mat using fgWrite in MATLAB (example in 't_cleanOpticRadiation.m')
        - open .mat in mrDiffusion and do initial cleaning by excluding tracts using plane ROIs 
                - in mrDiffusion: menu: ROIs > New polygon (make polygon around area to exclude)
                  then menu: Fibers > NOT with current ROI
        
7. Clean tracts using 't_cleanOpticRadiation.m' as example script.
