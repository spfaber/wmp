## This text file summarizes the procedural pipeline for carrying out DWI preprocessing using FSL.

####  Pestilli Lab at Indiana University Bloomington
##### Written by Sam Faber and Jack Zhang
updated 2016-02-29

_______________________________

###### FSL can work with multi b-value data, like the HCP data. 
###### Processes can be run from the command line or the interactive GUI.

**Step 1:** Make sure you have MRIcron downloaded so that you can use the dcm2nii function.   
        Ex. $ dcm2nii -g -d /path/to/*.dcm     
*      *Note: the -g creates a .nii.gz file and -d places the date in the name of the nifti file.*
     
**Step 2:** Make sure that the .bvec and .bval have delimiters that are spaces only and are arranged in ROW (# of directions) X COL (# of dimensions).

Once you have completed the fsl preprocessing, you may wish to make a dt6.mat from the fsl output. Before doing this using 'dtiMakeDt6FromFsl.m' first make sure that the T1 (which you will use for aligning your newly made B0) is ac-pc aligned. This can be done using *'feAutoACPCalign.m'* OR *'mrAnatSetNiftiXform.m'* from the VISTASOFT anatomy toolbox.

       
          
