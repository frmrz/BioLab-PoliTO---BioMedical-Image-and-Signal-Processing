#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 12:03:35 2022

@author: francesco
"""

import os
import tempfile

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from monai.transforms import (
    AddChanneld,
    Compose,
    LoadImaged,
    RandCropByPosNegLabeld,
    SpatialPadd,
    ToTensord,
)

from monai.utils import set_determinism

from monai.data import (
    DataLoader,
    CacheDataset,
    load_decathlon_datalist,
)

from PIL import Image

SEED = 46  # can be anything
set_determinism(SEED)

directory = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/DATA/DEVELOPMENT'

root_dir = tempfile.mkdtemp() if directory is None else directory
print(root_dir)

train_transforms = Compose(
    [
        LoadImaged(keys=["image", "label-a1"]),
        AddChanneld(keys=["image", "label-a1"]),

        RandCropByPosNegLabeld(
            keys=["image", "label-a1"],
            image_key="image",
            label_key="label-a1",
            spatial_size=(96, 96),
            pos=1,
            neg=3,
            num_samples=6,
            image_threshold=100,
            allow_smaller=True,
        ),
        SpatialPadd(
            keys=["image", "label-a1"],
            spatial_size=(96, 96),
            ),
        ToTensord(keys=["image", "label-a1"]),
    ]
)

val_transforms = Compose(
    [
        LoadImaged(keys=["image", "label-a1"]),
        AddChanneld(keys=["image", "label-a1"]),

        RandCropByPosNegLabeld(
            keys=["image", "label-a1"],
            image_key="image",
            label_key="label-a1",
            spatial_size=(96, 96),
            pos=1,
            neg=3,
            num_samples=6,
            image_threshold=100,
            allow_smaller=True,

        ),
        SpatialPadd(
            keys=["image", "label-a1", "label-a1s", "label-tum", "label-hybrid"],
            spatial_size=(96, 96),
            ),
        ToTensord(keys=["image", "label-a1", "label-a1s", "label-tum", "label-hybrid"]),
    ]
)


split_JSON = "/IMT-patchExtraction.json"
datasets = directory + split_JSON

datalist = load_decathlon_datalist(datasets, True, "training")
val_files = load_decathlon_datalist(datasets, True, "validation")

wks = 16

train_ds = CacheDataset(
    data=datalist, transform=train_transforms, cache_num=16, cache_rate=1.0, num_workers=wks,
)

train_loader = DataLoader(
    train_ds, batch_size=128, shuffle=True, num_workers=wks, pin_memory=True
)

val_ds = CacheDataset(
    data=val_files, transform=val_transforms, cache_num=16, cache_rate=1.0, num_workers=wks
)

val_loader = DataLoader(
    val_ds, batch_size=128, shuffle=False, num_workers=wks, pin_memory=True
)


# save patches
outImgTrain = os.path.join(directory, 'TRAINING', 'PATCH-96x96-IMG')
outMskTrain = os.path.join(directory, 'TRAINING', 'PATCH-96x96-A1')

for ii in tqdm(range(len(train_ds))):
    try:
        data = train_ds[ii]
        q = 0
        for batch_item in data:
            filename = os.path.split(batch_item["image_meta_dict"]["filename_or_obj"])[
                1].split('.')[0]
    
            img = batch_item["image"].numpy().squeeze()
            lab = batch_item["label-a1"].numpy().squeeze()
    
            pilImg = Image.fromarray(img).convert(mode='L').rotate(270)
            pilLab = Image.fromarray(lab).convert(mode='L').rotate(270)
    
            pilImg.save(os.path.join(
                outImgTrain, filename + '_{:04d}.png'.format(q)))
            pilLab.save(os.path.join(
                outMskTrain, filename + '_{:04d}.png'.format(q)))

            q += 1
    except:
        print('Error in file at index {}'.format(ii))

outImgVal = os.path.join(directory, 'VALIDATION', 'PATCH-96x96-IMG')
outMskVal = os.path.join(directory, 'VALIDATION', 'PATCH-96x96-A1')

for ii in tqdm(range(len(val_ds))):
    try:
        data = val_ds[ii]
        q = 0
        for batch_item in data:
            filename = os.path.split(batch_item["image_meta_dict"]["filename_or_obj"])[
                1].split('.')[0]
    
            img = batch_item["image"].numpy().squeeze()
            lab = batch_item["label-a1"].numpy().squeeze()
    
            pilImg = Image.fromarray(img).convert(mode='L').rotate(270)
            pilLab = Image.fromarray(lab).convert(mode='L').rotate(270)

            pilImg.save(os.path.join(
                outImgVal, filename + '_{:04d}.png'.format(q)))
            pilLab.save(os.path.join(
                outMskVal, filename + '_{:04d}.png'.format(q)))

            q += 1
    except:
        print('Error in file at index {}'.format(ii))
