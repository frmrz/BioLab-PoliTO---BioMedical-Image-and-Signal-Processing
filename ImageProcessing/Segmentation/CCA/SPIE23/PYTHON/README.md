# SPIE2023

This directory contains the Python and MATLAB scripts used for the SPIE2023 conference (10.1117/12.2653394). In brief, we defined a new hybrid ground truth segmentation starting from manual and computerized annotations achieving better generalization performances with respect to training with a purely manual ground truth.

## File Descriptions PYTHON

- train.py - script for Unet training, both for whole image and patch-based 
- inference.py - Whole image inference using segmentation_models
- patchBasedInference.py - Patch-based inference using MONAI APIs and Unet trained with segmentation_models
- myClasses.py - support classes for train.py and inference.py
- extractPatch.py - Patch extraction using MONAI to be used during training

## File Descriptions MATLAB
- consensusGT.m - creates consensus profiles for the data in the multicenter test set
- cropAndResizeUSimage.m  - performs the first cropping step of the paper
- evaluateProfileMetricsCommonSupport.m  - calculates the metrics using the manual and computerized LI and MA profiles
- find_US_Image_area.m - get the ROI of the image data inside a US image with auxiliary info
- fm_autocrop.m - automatic cropping of ultrasound images in DICOM3 format.
- get_dataset_summary.m - get mean and std of matrics stored in structures
- getLIMAfromMask.m - Extracts LIMA and MA profiles from a binary mask.

- makeHybridGT.m - creates the HybridGT starting from the 2 manual and 2 computerized profile
- plotProfiles.m - simple scripts for plotting profiles on the image
- profile2Mask.m - creates the IM complex mask starting from LI and MA profiles

