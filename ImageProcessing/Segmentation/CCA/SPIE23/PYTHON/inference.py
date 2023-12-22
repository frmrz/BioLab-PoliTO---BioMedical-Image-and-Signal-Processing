# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 16:52:20 2021

@author: franc
"""

import torch
from myClasses import Dataset, Dataset_test, visualize
from torch.utils.data import DataLoader
from tqdm import tqdm
import random
import albumentations as albu


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

def run():
    torch.multiprocessing.freeze_support()
    print('loop')

if __name__ == '__main__':
    run()
    
    import os
    os.environ['CUDA_VISIBLE_DEVICES'] = '0'
    
    import numpy as np
    import cv2
    
    ROOT = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/'
    DATA_DIR = os.path.join(ROOT, 'DATA', 'DEVELOPMENT')
    RESULTS_DIR =  os.path.join(ROOT, 'RESULTS', 'PYTHON')
    
    ## WHEN PATCHES CHANGE PREPROCESSING
    Experiment_name = 'UNET-480x480-v4'    
    checkpoints_path = os.path.join(RESULTS_DIR,Experiment_name,'checkpoints')
    
    x_test_dir = os.path.join(DATA_DIR, 'TEST', 'IMAGES-RESIZED')
    y_test_dir = os.path.join(DATA_DIR, 'TEST', 'MASKS-RESIZED')
        
    # Lets look at data we have
    
    dataset = Dataset(x_test_dir, y_test_dir, classes=['imt'])
    
    image, mask = dataset[4] # get some sample
    
    visualize(
        image=image.squeeze(), 
        muscle_mask=mask.squeeze(),
    )
    
    from torch.utils.tensorboard import SummaryWriter
    import segmentation_models_pytorch as smp
    
    ENCODER         = 'resnet50'
    ENCODER_WEIGHTS = 'imagenet'
    CLASSES         = ['imt']
    ACTIVATION      = 'sigmoid' # could be None for logits or 'softmax2d' for multicalss segmentation
    DEVICE          = 'cuda'
    ATTENTION       = None
    ACTIVATION       = None
    
    # create segmentation model with pretrained encoder
    model = smp.Unet (
        encoder_name=ENCODER, 
        encoder_depth = 5,
        encoder_weights=ENCODER_WEIGHTS, 
        classes=len(CLASSES), 
        in_channels=3,
        decoder_use_batchnorm=True,
        decoder_channels=(512,256, 128, 64, 32),
        activation=ACTIVATION,
    )
    
    model.classification_head = None

    #open encoders and change mean and std accordingly to dataset in franc\.conda\envs\torch-cuda11\Lib\site-packages\segmentation_models_pytorch/encoders/__init__.py
    preprocessing_fn = smp.encoders.get_preprocessing_fn(ENCODER, ENCODER_WEIGHTS) 
    
    test_dataset = Dataset(
        x_test_dir, 
        y_test_dir, 
        augmentation=get_validation_augmentation(), 
        preprocessing=get_preprocessing(preprocessing_fn),
        classes=CLASSES,
    )
    
    test_loader = DataLoader(test_dataset, batch_size=1, shuffle=False, num_workers=16, pin_memory=True, drop_last=True)

    loss = smp.utils.losses.DiceLoss()

    metrics = [
        smp.utils.metrics.IoU(threshold=0.5)
    ]
    
    optimizer = torch.optim.Adam([ 
        dict(params=model.parameters(), lr=0.0002),
    ])
    
    test_epoch = smp.utils.train.ValidEpoch(
        model, 
        loss=loss, 
        metrics=metrics, 
        device=DEVICE,
        verbose=True,
    )
    
    writer = SummaryWriter()

    
    # load best saved checkpoint
    best_model = torch.load(os.path.join(checkpoints_path,'best_model.pth'))
        
    #### TEST SET
    
    # logs = test_epoch.run(test_loader)

        
    Test_output_folder = os.makedirs(os.path.join(RESULTS_DIR,Experiment_name,'Test_output'), exist_ok = True)

    for i in tqdm(range(len(test_dataset))):
        
        image_vis = test_dataset[i][0].astype('uint8')
        image, gt_mask = test_dataset[i]
        
        gt_mask = gt_mask.squeeze()
        
        x_tensor = torch.from_numpy(image).to(DEVICE).unsqueeze(0)
        pr_mask = best_model.predict(x_tensor)
        pr_mask = (pr_mask.squeeze().cpu().numpy().round())
            
        cv2.imwrite(os.path.join(RESULTS_DIR,Experiment_name,'Test_output',os.listdir(x_test_dir)[i]),
                    pr_mask*255)
           
        
        