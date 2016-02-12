#!/bin/bash

## This shell script calls mrtrix functions to convert file types, create masks, and track between two ROIs (the LGN and V1).
## Information and examples can be found at: http://jdtournier.github.io/mrtrix-0.2/tractography/index.html .

## Sam Faber
## updated : 11-15-2015

## has a module now
## export PATH=$PATH:/N/soft/rhel6/mrtrix/0.2.12_test/bin


SUBJ=$1

TOPDIR=/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/mrtrix_track_between_rois/$SUBJ
OUTDIR=/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/mrtrix_track_between_rois/$SUBJ/mrtrix_results


## convert dwi's 
mrconvert $TOPDIR/dwi_data_b2000_aligned_trilin.nii.gz $OUTDIR/dwi_data_b2000_aligned_trilin_dwi.mif


## make mask from DWI data instead of importing mrDiffusion dt6
average $OUTDIR/dwi_data_b2000_aligned_trilin_dwi.mif -axis 3 - | threshold - - | median3D - - | median3D - $OUTDIR/dwi_data_b2000_aligned_trilin_brainmask.mif


## fit tensors
dwi2tensor $OUTDIR/dwi_data_b2000_aligned_trilin_dwi.mif -grad $TOPDIR/dwi_data_b2000_aligned_trilin.b $OUTDIR/dwi_data_b2000_aligned_trilin_dt.mif 


## create FA image
tensor2FA $OUTDIR/dwi_data_b2000_aligned_trilin_dt.mif - | mrmult - $OUTDIR/dwi_data_b2000_aligned_trilin_brainmask.mif $OUTDIR/dwi_data_b2000_aligned_trilin_fa.mif


## create eigenvector map
tensor2vector $OUTDIR/dwi_data_b2000_aligned_trilin_dt.mif - | mrmult - $OUTDIR/dwi_data_b2000_aligned_trilin_fa.mif $OUTDIR/dwi_data_b2000_aligned_trilin_ev.mif


## convert wm mask
mrconvert $TOPDIR/wm_mask.nii.gz $OUTDIR/dwi_data_b2000_aligned_trilin_wm.mif


## erodes brainmask - removes extreme artifacts (w/high FA) and creates an FA image, AND single fiber mask (with the artifactually high FA's)
erode $OUTDIR/dwi_data_b2000_aligned_trilin_brainmask.mif -npass 3 - | mrmult $OUTDIR/dwi_data_b2000_aligned_trilin_fa.mif - - | threshold - -abs 0.7 $OUTDIR/dwi_data_b2000_aligned_trilin_sf.mif


## estimates response function - estimates profile of diffusion wavelet in region of sf
estimate_response $OUTDIR/dwi_data_b2000_aligned_trilin_dwi.mif $OUTDIR/dwi_data_b2000_aligned_trilin_sf.mif -lmax 10 -grad $TOPDIR/dwi_data_b2000_aligned_trilin.b $OUTDIR/dwi_data_b2000_aligned_trilin_response.txt


## compute CSD 
csdeconv $OUTDIR/dwi_data_b2000_aligned_trilin_dwi.mif $OUTDIR/dwi_data_b2000_aligned_trilin_response.txt -lmax 10 -mask $OUTDIR/dwi_data_b2000_aligned_trilin_wm.mif $OUTDIR/CSD10.mif


## does CSD for creating fiber maps
csdeconv $OUTDIR/dwi_data_b2000_aligned_trilin_dwi.mif -grad $TOPDIR/dwi_data_b2000_aligned_trilin.b $OUTDIR/dwi_data_b2000_aligned_trilin_response.txt -lmax 10 -mask $OUTDIR/dwi_data_b2000_aligned_trilin_wm.mif $OUTDIR/CSD10.mif


## convert lh .nii ROIs to .mif
mrconvert $TOPDIR/ROIs/lh_LGN.nii.gz $OUTDIR/lh_LGN.mif
mrconvert $TOPDIR/ROIs/lh_V1_label_smooth3mm.nii.gz $OUTDIR/lh_V1.mif


## merge ROIs for creating seed region
mradd $OUTDIR/lh_LGN.mif $OUTDIR/lh_V1.mif $OUTDIR/lh_LGN2V1.mif


## convert rh .nii ROIs to .mif
mrconvert $TOPDIR/ROIs/rh_LGN.nii.gz $OUTDIR/rh_LGN.mif
mrconvert $TOPDIR/ROIs/rh_V1_label_smooth3mm.nii.gz $OUTDIR/rh_V1.mif


## merge ROIs for creating seed region
mradd $OUTDIR/rh_LGN.mif $OUTDIR/rh_V1.mif $OUTDIR/rh_LGN2V1.mif


## Perform multiple types of tracking between same ROIs 

## create a probabilistic lh tract between the two ROIs using CSD
streamtrack -seed $OUTDIR/lh_LGN2V1.mif -mask $OUTDIR/dwi_data_b2000_aligned_trilin_wm.mif -grad $TOPDIR/dwi_data_b2000_aligned_trilin.b -include $OUTDIR/lh_V1.mif -include $OUTDIR/lh_LGN.mif SD_PROB $OUTDIR/CSD10.mif $OUTDIR/left_optic_radiation_PCSD.tck -number 20000 -maxnum 20000000

## create a probabilistic rh tract between the two ROIs using CSD
streamtrack -seed $OUTDIR/rh_LGN2V1.mif -mask $OUTDIR/dwi_data_b2000_aligned_trilin_wm.mif -grad $TOPDIR/dwi_data_b2000_aligned_trilin.b -include $OUTDIR/rh_V1.mif -include $OUTDIR/rh_LGN.mif SD_PROB $OUTDIR/CSD10.mif $OUTDIR/right_optic_radiation_PCSD.tck -number 20000 -maxnum 20000000


## create a deterministic tract between the two ROIs using DT
#streamtrack -seed $OUTDIR/rh_LGN2V1.mif -mask $OUTDIR/dwi_data_b2000_aligned_trilin_wm.mif -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b -#include $OUTDIR/rh_V1.mif -include $OUTDIR/rh_LGN.mif DT_STREAM $OUTDIR/dwi_data_b2000_aligned_trilin_dwi.mif $OUTDIR/optic_tract_DDT.tck -#number 1000 -maxnum 1000


## create a probabilistic tract between the two ROIs using DT
#streamtrack -seed $OUTDIR/rh_LGN2V1.mif -mask $OUTDIR/dwi_data_b2000_aligned_trilin_wm.mif -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b -#include $OUTDIR/rh_V1.mif -include $OUTDIR/rh_LGN.mif DT_PROB $OUTDIR/dwi_data_b2000_aligned_trilin_dwi.mif $OUTDIR/optic_tract_PDT.tck 


## create a deterministic tract between the two ROIs using CSD
#streamtrack -seed $OUTDIR/rh_LGN2V1.mif -mask $OUTDIR/dwi_data_b2000_aligned_trilin_brainmask.mif -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b #-include $OUTDIR/rh_V1.mif -include $OUTDIR/rh_LGN.mif SD_STREAM CSD8.mif $OUTDIR/optic_tract_DCSD.tck -number 1000 -maxnum 1000000


# -seed 3.52,-74.73,4.84,4 (seed for V1)
# -seed 30.31,-13.18,-8.91,4 (seed for LGN)
