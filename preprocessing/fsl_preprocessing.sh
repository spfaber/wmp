#!/bin/bash

# This script preprocesses raw diffusion and anatomical data using fsl.

# Written by Sam Faber, Brent McPherson and Jack Zhang 
# Modified from original stream by Hu Chen
# Updated 2016-04-18

# Make sure we are using the right fsl module version
export PATH=$PATH:/N/soft/rhel6/fsl/5.0.6/bin


# Set paths to data
ROOTDIR=/N/dc2/projects/lifebid/HCP/Sam
RAWDIR=$ROOTDIR/sophia_project/pilot/diffusion/raw
PARAMDIR=$ROOTDIR/preprocessing

# output directory
OUTDIR=/N/dc2/projects/lifebid/HCP/Sam/sophia_project/pilot/diffusion/preprocessed

#STEP 1: Part A
## Make nifti of b0 volumes for AP acquisitions
## AP 
# do fslsplit to separate AP volumes
cd $RAWDIR/AP_images
fslsplit dwi_b1000_b2000_AP.nii.gz $RAWDIR/AP_images/

# use bvals to keep only b0 volumes
# read txt file
AP_bvalString=$(cat $RAWDIR/AP_images/dwi_b1000_b2000_AP.bval |tr "\n" " ")
# make array of textfile elements
AP_bval=($AP_bvalString)
#echo "${AP_bval[@]}"

# get indices of b0's
value='0'
for i in "${!AP_bval[@]}"; do
   if [[ "${AP_bval[$i]}" = "${value}" ]]; then
       AP_b0_inds+=("${i}"); 
   fi
done
#echo "${AP_b0_inds[@]}"

# create array of volumes files
padtowidth=4
for i in "${!AP_b0_inds[@]}"; do 
#vols+=$( printf "%0*d\n" $padtowidth "${AP_b0_inds[$i]}" |tr "\n" " ");
vols+=$( printf "%0*d " $padtowidth "${AP_b0_inds[$i]}");
done

# splits into an indexable array
read -a arr <<< $vols

# checking
#echo "${vols}" 
#echo "${arr[0]}"
#echo "${arr[5]}"

for x in "${arr[@]}" 
do
    b0files+=$(echo vol$x.nii.gz " ")
done

#echo $b0files  
# splits into an idexable array
read -a b0files <<< $b0files

echo ${b0files[@]}

# merge them 
fslmerge -t nodif_AP ${b0files[@]} 

# STEP 1: Part B
## Make nifti of b0 volumes for PA acquisitions
# do fslsplit to separate PA volumes
cd $RAWDIR/PA_images
fslsplit dwi_b1000_b2000_PA.nii.gz $RAWDIR/PA_images/

# use bvals to keep only b0 volumes
# read txt file
PA_bvalString=$(cat $RAWDIR/PA_images/dwi_b1000_b2000_PA.bval |tr "\n" " ")
# make array of textfile elements
PA_bval=($PA_bvalString)
#echo "${PA_bval[@]}" # check bval array

# get indices of b0's
value='0'
for i in "${!PA_bval[@]}"; do
   if [[ "${PA_bval[$i]}" = "${value}" ]]; then
       AP_b0_inds+=("${i}"); 
   fi
done
#echo "${PA_b0_inds[@]}" # check b0 Inds

# create array of volumes files
padtowidth=4
for i in "${!PA_b0_inds[@]}"; do 
#vols+=$( printf "%0*d\n" $padtowidth "${PA_b0_inds[$i]}" |tr "\n" " ");
vols+=$( printf "%0*d " $padtowidth "${PA_b0_inds[$i]}");
done

# splits into an indexable array
read -a arr <<< $vols

# checking
#echo "${vols}" 
#echo "${arr[0]}"
#echo "${arr[5]}"

for x in "${arr[@]}" 
do
    b0files+=$(echo vol$x.nii.gz " ")
done

#echo $b0files  
# splits into an idexable array
read -a b0files <<< $b0files

echo ${b0files[@]}

# merge them 
fslmerge -t nodif_PA ${b0files[@]} 


# Change back to code directory when finished making nodif files
cd $PARAMDIR

# STEP 2: Par A Perform eddy correction, average b0's, and do brain extraction for AP images
echo "AP image preprocessing..."
eddy_correct $RAWDIR/AP_images/nodif_AP $OUTDIR/nodif_AP_ec 0
fslmaths $OUTDIR/nodif_AP_ec -Tmean $OUTDIR/nodif_AP_mean
bet $OUTDIR/nodif_AP_mean $OUTDIR/nodif_AP_brain -f 0.4 -g 0 -m

# Step 2: PArt B Perform eddy correction, average b0's, and do brain extraction for PA images
echo "PA image preprocessing..."
eddy_correct $RAWDIR/PA_images/nodif_PA $OUTDIR/nodif_PA_ec 0
fslmaths $OUTDIR/nodif_PA_ec -Tmean $OUTDIR/nodif_PA_mean
bet $OUTDIR/nodif_PA_mean $OUTDIR/nodif_PA_brain -f 0.4 -g 0 -m

# Step 3: Merge AP/PA mean b0 images
fslmerge -t $OUTDIR/b0_images $OUTDIR/nodif_AP_mean $OUTDIR/nodif_PA_mean

# Step 4: Do topup on AP/PA images
echo "Performing topup..."
topup --imain=$RAWDIR/b0_images --datain=$PARAMDIR/acq_params.txt --config=$PARAMDIR/b02b0.cnf --out=$OUTDIR/topup_results --fout=$OUTDIR/topup_field --iout=$OUTDIR/topup_unwarped_images
echo "done."

# Step 5: Average unwarped b0's
fslmaths $OUTDIR/topup_unwarped_images -Tmean $OUTDIR/hifi_b0 (use fslmaths -add -div N (may need to separate first - use fslsplit))

# Step 6: Do averaged b0 brain extraction 
bet $OUTDIR/hifi_b0 $OUTDIR/hifi_b0_brain -m

# Step 7: Do eddy correction on all directions for AP and PA images
# A. Create AP index file for eddy function
indx=""
for ((i=1; i<=143; i+=1)); do indx="$indx 1"; done
echo $indx > $PARAMDIR/index.txt

# B. Perform AP eddy current correction
eddy --imain=$RAWDIR/dwi_b1000_b2000_AP.nii.gz --mask=$OUTDIR/hifi_b0_brain_mask --acqp=$PARAMDIR/acq_params.txt --index=$PARAMDIR/index.txt --bvecs=$RAWDIR/dwi_b1000_b2000_AP.bvec --bvals=$RAWDIR/dwi_b1000_b2000_AP.bval --topup=$OUTDIR/topup_results --out=$OUTDIR/eddy_corrected_data_AP

# C. Create PA index file
indx=""
for ((i=1; i<=143; i+=1)); do indx="$indx 2"; done
echo $indx > $PARAMDIR/index2.txt

# D. Do PA eddy 
eddy --imain=$RAWDIR/dwi_b1000_b2000_PA.nii.gz --mask=$OUTDIR/hifi_b0_brain_mask --acqp=$PARAMDIR/acq_params.txt --index=$PARAMDIR/index2.txt --bvecs=$RAWDIR/dwi_b1000_b2000_PA.bvec --bvals=$RAWDIR/dwi_b1000_b2000_PA.bval --topup=$OUTDIR/topup_results --out=$OUTDIR/eddy_corrected_data_PA

# Step 8: Merge eddy corrected AP and PA 
fslmerge -t $OUTDIR/data $OUTDIR/eddy_corrected_data_PA $OUTDIR/eddy_corrected_data_AP

# Make bvecs and bvals that are twice as long so we can compute tensors
paste $RAWDIR/dwi_b1000_b2000_AP.bval dwi_b1000_b2000_AP.bval > dwi_b1000_b2000_AP_PA.bval
paste $RAWDIR/dwi_b1000_b2000_AP.bvec dwi_b1000_b2000_AP.bvec > dwi_b1000_b2000_AP_PA.bvec

# Do b0 brain extraction 
bet $OUTDIR/topup_unwarped_images $OUTDIR/nodif_brain -f 0.3 -g 0 -m

# Fit tensors
dtifit --data=$OUTDIR/data --out=$OUTDIR/dti --mask=$OUTDIR/nodif_brain_mask --bvecs=$RAWDIR/dwi_b1000_b2000_AP_PA.bvec --bvals=$RAWDIR/dwi_b1000_b2000_AP_PA.bval


# should we be computing the tensors on data that is the same size as the original data? why do we merge eddy corrected topup output AP and PA data?

# ++++++++++++++++++++++++++++++++++++++++++
# Step 3: Fill Holes in Brain Mask
#echo "Filling holes in binary mask..."
#fslmaths dwi_data_b1000_scan1_ec_bm_mask.nii.gz -fillh dwi_data_b1000_scan1_ec_bm_mask_fh.nii.gz


# Step 4: Calculate diffusion tensors - for this step, make sure that the .bvec and .bval have delimiters that are spaces only and are arranged in ROW X COL
# To do this, read into MATLAB your .bval and .bvec files (using M = dlmread) and save them with a new name out specifying a delimiter of space ' '. E.g. dlmwrite('abc.dat',M,'delimiter',' ');
#echo "Calculating diffusion tensors..."
#dtifit -k $OUTDIR/dwi_data_b1000_scan1_ec -o $OUTDIR/dwi_data_b1000_scan1_ec -m $OUTDIR/dwi_data_b1000_scan1_ec_bm_mask_fh -r $RAWDIR/dwi_data_b1000_scan1.bvec -b $RAWDIR/dwi_data_b1000_scan1.bval
