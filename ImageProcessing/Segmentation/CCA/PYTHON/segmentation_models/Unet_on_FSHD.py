# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 10:55:23 2021

@author: franc
"""

import torch
from torch.utils.data import DataLoader
from myClasses import Dataset, Dataset_test, visualize

import os
import segmentation_models_pytorch as smp
from tqdm import tqdm
import cv2

import albumentations as albu

MAIN_DIR = 'C://Users//franc//Desktop//POLITO//FSHD'
DATA_DIR = os.path.join(MAIN_DIR,'DATA', 'imageData_FSHD')

Experiment = 'Unet50_pretrained'

Images_dir = os.path.join(DATA_DIR)
Labels_dir = os.path.join(DATA_DIR)

def get_validation_augmentation():
    """Add paddings to make image shape divisible by 32"""
    test_transform = [
        # albu.PadIfNeeded(240, 256)
    ]
    return albu.Compose(test_transform)

def to_tensor(x, **kwargs):
    return x.transpose(2, 0, 1).astype('float32')


def get_preprocessing(preprocessing_fn):
    """Construct preprocessing transform
    
    Args:
        preprocessing_fn (callbale): data normalization function 
            (can be specific for each pretrained neural network)
    Return:
        transform: albumentations.Compose
    
    """
    
    _transform = [
        albu.Lambda(image=preprocessing_fn),
        albu.Lambda(image=to_tensor, mask=to_tensor),
    ]
    return albu.Compose(_transform)
    
# Lets look at data we have

# dataset = Dataset(Images_dir, Labels_dir, classes=['muscle'])

# image, mask = dataset[4] # get some sample
# # image=image.transpose(1,2,0)

# visualize(
#     # image=image.transpose(1,2,0), 
#     image=image.squeeze(), 
#     muscle_mask=mask.squeeze(),
# )
    
    
ENCODER         = 'resnet50'
ENCODER_WEIGHTS = 'imagenet'
ACTIVATION      = 'sigmoid'

preprocessing_fn = smp.encoders.get_preprocessing_fn(ENCODER, ENCODER_WEIGHTS) 

dataset = Dataset(Images_dir, Labels_dir,
                   augmentation=get_validation_augmentation(), 
                   preprocessing=get_preprocessing(preprocessing_fn),
                  classes=['muscle'])

# loader = DataLoader(dataset, batch_size=2, shuffle=False, num_workers=0, pin_memory=True)

best_model = torch.load(os.path.join(MAIN_DIR,'CODE','PYTHON','models','Unet.pth'))
loss       = smp.utils.losses.DiceLoss()
metrics    = [smp.utils.metrics.IoU(threshold=0.5)]
DEVICE     = 'cuda'

# runner = smp.utils.train.ValidEpoch(
#         model=best_model,
#         loss=loss,
#         metrics=metrics,
#         device=DEVICE,
#     ) 

# logs = runner.run(loader)
folder = os.path.join(MAIN_DIR,'RESULTS','Unet50_pretrained')
# output_folder = os.makedirs(os.path.join(DATA_DIR,'Unet_results',folder), exist_ok = True)

for i in tqdm(range(len(dataset))):
    
    image_vis = dataset[i][0].astype('uint8')
    image, gt_mask = dataset[i]
    
    gt_mask = gt_mask.squeeze()
    
    x_tensor = torch.from_numpy(image).to(DEVICE).unsqueeze(0)
    pr_mask = best_model.predict(x_tensor)
    pr_mask = (pr_mask.squeeze().cpu().numpy().round())
    pr_mask = pr_mask.astype('uint8')*255
        
    pt = dataset.images_fps[i]
    pz_id, temp_id = pt.split(os.sep)[3], pt.split(os.sep)[5]
    
    output_folder = os.path.join(folder, pz_id)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    cv2.imwrite(os.path.join(output_folder,temp_id), pr_mask)