# SPIE2023

This directory contains the Python and MATLAB scripts used for the SPIE2023 conference (10.1117/12.2653394). In brief, we defined a new hybrid ground truth segmentation starting from manual and computerized annotations achieving better generalization performances with respect to training with a purely manual ground truth.

## File Descriptions PYTHON

- train.py - script for Unet training, both for whole image and patch-based 
- inference.py - Whole image inference using segmentation_models
- patchBasedInference.py - Patch-based inference using MONAI APIs and Unet trained with segmentation_models
- myClasses.py - support classes for train.py and inference.py
- extractPatch.py - Patch extraction using MONAI to be used during training

## File Descriptions PYTHON

