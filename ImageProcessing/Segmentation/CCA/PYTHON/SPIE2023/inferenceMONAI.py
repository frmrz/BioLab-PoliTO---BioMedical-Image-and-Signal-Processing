#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 22:18:34 2022

@author: francesco
"""

import os
import tempfile

import matplotlib.pyplot as plt
import numpy as np
# from tqdm import tqdm
import torch 
# import skimage
# import cv2

# from monai.losses import DiceCELoss
from monai.inferers import sliding_window_inference
from monai.transforms import (
    AsDiscrete,
    # DataStatsd,
    AddChanneld,
    Compose,
    LoadImaged,
    Transposed,
    # EnsureChannelFirstd,
    ToTensord,
    RepeatChanneld,
    ScaleIntensityRanged,
    ScaleIntensityRangePercentilesd,
    EnsureType,

)

from monai.config import print_config
from monai.metrics import DiceMetric
# from monai.networks.nets import UNETR

from monai.data import (
    DataLoader,
    CacheDataset,
    load_decathlon_datalist,
    decollate_batch,
)

from PIL import Image
import segmentation_models_pytorch as smp
from skimage.measure import regionprops

from monai.utils.enums import BlendMode
from monai.utils import set_determinism

os.environ['CUDA_LAUNCH_BLOCKING']="1"

print_config()

directory = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/'

root_dir = tempfile.mkdtemp() if directory is None else directory
print(root_dir)

val_transforms = Compose(
    [
        LoadImaged(keys=["image", "label"]),
        AddChanneld(keys=["image", "label"]),
        
        ScaleIntensityRanged(
            keys=["image", "label"],
            a_min=0,
            a_max=255,
            b_min=0.0,
            b_max=1.0,
            clip=True,
        ),
        ScaleIntensityRangePercentilesd(
            keys=["image"],
            lower=1,
            upper=99,
            b_min=0.0,
            b_max=1.0,
            clip=True,
        ),        
        RepeatChanneld(keys=["image"],repeats=3),
        Transposed(keys=["image", "label"],indices=(0,2,1)),
        # DataStatsd(keys=['image', 'label'], data_value=False),

        ToTensord(keys=["image", "label"]),
    ]
)

SEED = 46  # can be anything
set_determinism(SEED)

data_dir = os.path.join(root_dir, 'DATA', 'DEVELOPMENT')
split_JSON = "/IMT-TEST-patchInference.json"
datasets = data_dir + split_JSON

val_files = load_decathlon_datalist(datasets, True, "validation")

wks = 16

val_ds = CacheDataset(
    data=val_files, transform=val_transforms, cache_num=128, cache_rate=1.0, num_workers=wks
)

val_loader = DataLoader(
    val_ds, batch_size=1, shuffle=False, num_workers=wks, pin_memory=True
)


case_num = 3

dataPlot = val_ds[case_num]

q=0

# loop through the length of tickers and keep track of index
for n in range(0,11,2):
    img = dataPlot["image"]
    label = dataPlot["label"]

    # add a new subplot iteratively
    ax = plt.subplot(2,1,1)
    plt.imshow(img[0, :, :].detach().cpu(), cmap="gray")

    # add a new subplot iteratively
    ax = plt.subplot(2,1,2)
    plt.imshow(label[0, :, :].detach().cpu(), cmap="gray")
    q=q+1

plt.show()


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

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
    decoder_channels=(512,256, 128, 64,32),
    activation=ACTIVATION,
)

model.classification_head = None

# post_pred = Compose([EnsureType(), AsDiscrete(argmax=True, to_onehot=1)])
# post_pred = Compose([EnsureType(), AsDiscrete(to_onehot=1)])
post_pred = Compose([EnsureType(), AsDiscrete(threshold=0.5)])
post_label = Compose([EnsureType()])
# post_label = Compose([EnsureType(), AsDiscrete(to_onehot=1)])
# ck = Compose([EnsureType(), DataStats(prefix='Data')])

dice_metric = DiceMetric(include_background=True, reduction="mean", get_not_nans=False)
global_step = 0
dice_val_best = 0.0
global_step_best = 0
epoch_loss_values = []
metric_values = []

## BBOX DATA
bboxDir = os.path.join(data_dir, 'TEST', 'BBOX')
resizedImgDir = os.path.join(data_dir, 'TEST', 'IMAGES-RESIZED')
resizedMskDir = os.path.join(data_dir, 'TEST', 'MASKS-RESIZED')
detectdMskDir = os.path.join(data_dir, 'TEST', 'MASKS-DETECT')


# experiments = ['UNET-64x64-v0','UNET-64x64-v1','UNET-64x64-v2','UNET-64x64-v3','UNET-64x64-v4']
experiments = ['UNET-96x96-v0','UNET-96x96-v1','UNET-96x96-v2','UNET-96x96-v3','UNET-96x96-v4']

for experiment in experiments:
    model = torch.load(os.path.join(root_dir,'RESULTS','PYTHON',
                                                  experiment,
                                                  'checkpoints',
                                                  "best_model.pth"))
    
    
    model.eval()
    with torch.no_grad():
        
        dice_vals = list()
    
        for i, val_data in enumerate(val_loader):
                
            val_inputs, val_labels = (val_data["image"].to(device),val_data["label"].to(device))
            
            filename = os.path.split(val_data["image_meta_dict"]["filename_or_obj"][0])[1]
            
            print(' Processing {}'.format(filename))
            
            roi_size = (96, 96)
            sw_batch_size = 128
            
            val_outputs = sliding_window_inference(val_inputs, roi_size,
                                                   sw_batch_size,
                                                   model,
                                                   overlap=0.75,
                                                   mode=BlendMode.GAUSSIAN,
                                                   sigma_scale=0.125)
            
            val_outputs = torch.sigmoid(val_outputs)
            
            val_outputs = [post_pred(j) for j in decollate_batch(val_outputs)]
            val_labels  = [post_label(j) for j in decollate_batch(val_labels)]        
            
            # compute metric for current iteration
            dice_metric(y_pred=val_outputs, y=val_labels)
            
            dice = dice_metric.aggregate().item()
            dice_vals.append(dice)
            dice_metric.reset()
            
            ## Save DETECT
            val_outputs = val_outputs[0].detach().cpu().numpy()
            
            pilImg = Image.fromarray(val_outputs[0]*255).convert(mode='L')
            pilImg.save(os.path.join(directory,'RESULTS','PYTHON',
                                     experiment,'DETECT-OUTPUT-TEST',
                                     filename))
                
            val_outputs = val_outputs.squeeze()[:-1,:-1]
            pilImg = Image.fromarray(val_outputs*255).convert(mode='L')
    
            ## Make and Save RESIZED
            
            temp_bbox = Image.open(os.path.join(bboxDir, filename))
            temp_resized_img = Image.open(os.path.join(resizedImgDir, filename))
            temp_resized_msk = Image.open(os.path.join(resizedMskDir, filename))
            
            np_bbox = np.uint8(np.asarray(temp_bbox)/255)
            np_resized_msk = np.asarray(temp_resized_msk)
                    
            stats = regionprops(np_bbox)[0].bbox
            
            tempOut = Image.new(mode="L", size=(480, 480))
            Image.Image.paste(tempOut,pilImg,box=(stats[1],stats[0]))
            
            tempOut.save(os.path.join(directory,'RESULTS','PYTHON',
                                     experiment,'Test_output',
                                     filename))
    
                
        mean_dice_val = np.mean(dice_vals)
    
    
    
    
    
    
    


