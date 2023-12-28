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
- cropAndResizeUSimage.m - performs the first cropping step of the paper
- evaluateProfileMetricsCommonSupport.m - calculates the metrics using the manual and computerized LI and MA profiles
- find_US_Image_area.m - get the ROI of the image data inside a US image with auxiliary info
- fm_autocrop.m - automatic cropping of ultrasound images in DICOM3 format.
- get_dataset_summary.m - get mean and std of matrics stored in structures
- getLIMAfromMask.m - Extracts LIMA and MA profiles from a binary mask.
- km_CommonSupport.m - calculates the common support between two sets of points.
- LI_MA_interp.m - Interpolates LI and MA profiles
- LI_MA_stats.m - calculates various statistics for image segmentation evaluation.
- LI_MA_stats_light.m - calculates various statistics for LI and MA tracing precision assesment.
- LI_MA_stats_not_interp.m - calculates various statistics for LI and MA tracing precision assesment when profile are not interpolated
- makeHybridGT.m - creates the HybridGT for training starting from the 2 manual and 1 computerized profile
- padding_rectangular.m - pads the input image I to make it a square image by adding zeros on the sides.
- PolyDistMethod.m - Calculates the distance between two profiles B1 and B2 by using the Polyline Distance Method 
- plotProfiles.m - simple scripts for plotting profiles on the image
- profile2Mask.m - creates the IM complex mask starting from LI and MA profiles
- Struct_Empty_To_Nan.m - replaces empty values in a structure array with NaN.
- TurnColumn.m - Converts a vector into a row vector
- TurnRow.m - Converts the input vector x into a row vector X.
- write_txt_file.m - Writes a vector with numerical data to a file with a certain filename in a certain directory
