#!/bin/bash

# This script performs basic preprocessing of raw diffusion data using fsl.

# Written by Sam Faber and Jack Zhang 
# Pestilli Lab at Indiana University Bloomington
# 2016-02-08


# Make sure we are using the right fsl module version
export PATH=$PATH:/N/soft/rhel6/fsl/5.0.8/bin


# Set path to data
RAWDIR=/N/dc2/projects/lifebid/HCP/Sam/kid_data/180743/raw


# Set output path
fslDIR=/N/dc2/projects/lifebid/HCP/Sam/kid_data/180743/fsl


# Step 1: Perform eddy current correction of the data using the b0 as the reference volume
echo Performing eddy current correction...
eddy_correct $RAWDIR/dwi_data_b1000_scan1 $fslDIR/dwi_data_b1000_scan1_ec 0


# Step 2: Perform brain extraction to create a binary mask of the brain to be used for fitting of the diffusion tensors.
#         Lowering the fractional intensity threshold 'f' will give more liberal brain outline estimates.
#         Setting the center of the brain 'c' will ensure that the masking is done from the correct origin.
echo Creating binary mask...
bet $fslDIR/dwi_data_b1000_scan1_ec $fslDIR/dwi_data_b1000_scan1_ec_bm -m -f 0.3 -c 48 48 30


# Step 3: Fill holes in brain mask
echo Filling holes in binary mask...
fslmaths dwi_data_b1000_scan1_ec_bm_mask.nii.gz -fillh dwi_data_b1000_scan1_ec_bm_mask_fh.nii.gz


# Step 4: Calculate diffusion tensors from the data using the binary brain mask, bvecs and bvals 
echo Calculating diffusion tensors...
dtifit -k $fslDIR/dwi_data_b1000_scan1_ec -o $fslDIR/dwi_data_b1000_scan1_ec -m $fslDIR/dwi_data_b1000_scan1_ec_bm_mask_fh -r $RAWDIR/dwi_data_b1000_scan1.bvec -b $RAWDIR/dwi_data_b1000_scan1.bval





