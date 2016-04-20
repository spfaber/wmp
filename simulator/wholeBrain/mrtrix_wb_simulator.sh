## Sam Faber
## 01-18-2016


# This script takes an input of a subject name or number
SUBJ=$1 

DWIDIR=/N/dc2/projects/lifebid/HCP/Sam/$SUBJ/diffusion_data
TOPDIR=/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/simulator/wholeBrain/
OUTDIR=/N/dc2/projects/lifebid/HCP/Sam/matlab_code/wmp/simulator/wholeBrain/mrtrix_results


## convert dwi's 
#mrconvert $DWIDIR/dwi_data_b2000_aligned_trilin.nii.gz $OUTDIR/dwi_data_b2000_aligned_trilin.mif


## make brainmask
#average $OUTDIR/dwi_data_b2000_aligned_trilin.mif -axis 3 - | threshold - - | median3D - - | median3D - $OUTDIR/dwi_data_b2000_aligned_trilin_brainmask.mif


## fit tensors
#dwi2tensor $OUTDIR/dwi_data_b2000_aligned_trilin.mif -grad $TOPDIR/dwi_data_b2000_aligned_trilin.b $OUTDIR/dwi_data_b2000_aligned_trilin_dt.mif 


## create FA image
#tensor2FA $OUTDIR/dwi_data_b2000_aligned_trilin_dt.mif - | mrmult - $OUTDIR/dwi_data_b2000_aligned_trilin_brainmask.mif $OUTDIR/dwi_data_b2000_aligned_trilin_fa.mif


## create eigenvector map
#tensor2vector $OUTDIR/dwi_data_b2000_aligned_trilin_dt.mif - | mrmult - $OUTDIR/dwi_data_b2000_aligned_trilin_fa.mif $OUTDIR/dwi_data_b2000_aligned_trilin_ev.mif


## convert wm mask
mrconvert $DWIDIR/iso_vox_wm_mask.nii.gz $OUTDIR/wm_mask.mif


## erodes brainmask - removes extreme artifacts (w/high FA) and creates an FA image, AND single fiber mask (with the artifactually high FA's)
#erode $OUTDIR/dwi_data_b2000_aligned_trilin_brainmask.mif -npass 3 - | mrmult $OUTDIR/dwi_data_b2000_aligned_trilin_fa.mif - - | threshold - -abs 0.7 $OUTDIR/dwi_data_b2000_aligned_trilin_sf.mif


## estimates response function - estimates profile of diffusion wavelet in region of sf
#estimate_response $OUTDIR/dwi_data_b2000_aligned_trilin.mif $OUTDIR/dwi_data_b2000_aligned_trilin_sf.mif -lmax 10 -grad $TOPDIR/dwi_data_b2000_aligned_trilin.b $OUTDIR/dwi_data_b2000_aligned_trilin_response.txt



## does CSD for creating fiber maps
#csdeconv $OUTDIR/dwi_data_b2000_aligned_trilin.mif -grad $TOPDIR/dwi_data_b2000_aligned_trilin.b $OUTDIR/dwi_data_b2000_aligned_trilin_response.txt -lmax 10 -mask $OUTDIR/dwi_data_b2000_aligned_trilin_brainmask.mif $OUTDIR/CSD10.mif


## Perform whole brain tracking
## create a probabilistic whole brain connectome with the left and right ORs excluded using CSD
streamtrack -seed $OUTDIR/wm_mask.mif -mask $OUTDIR/wm_mask.mif -grad $TOPDIR/dwi_data_b2000_aligned_trilin.b -number 100000 -maxnum 0 SD_PROB $OUTDIR/CSD10.mif $OUTDIR/whole_brain.tck 


