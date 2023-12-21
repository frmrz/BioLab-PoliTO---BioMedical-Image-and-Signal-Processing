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
# from Models import R2AttU_Net, R2U_Net, AttU_Net
# from network_models import AttU_Net_resnet
import albumentations as albu

def get_training_augmentation(p_aug=1):
    
    p_transform = random.random()
    
    if p_transform < p_aug:
        train_transform = [
    
            albu.HorizontalFlip(p=0.5),
    
            albu.ShiftScaleRotate(scale_limit=0.1, rotate_limit=10, shift_limit=0.1, p=1, border_mode=cv2.BORDER_WRAP),
    
            albu.augmentations.transforms.GaussNoise(var_limit=(10.0, 50.0), mean=0, per_channel=False, always_apply=False, p=0.3),
            albu.augmentations.transforms.CLAHE (clip_limit=4.0, tile_grid_size=(8, 8), always_apply=False, p=0.3),
            albu.augmentations.transforms.MultiplicativeNoise (multiplier=(0.9, 1.1), per_channel=False, elementwise=False, always_apply=False, p=0.2),
            albu.augmentations.transforms.OpticalDistortion (distort_limit=0.05, shift_limit=0.05, interpolation=1, border_mode=4, value=None, mask_value=None, always_apply=False, p=0.2),                   
            albu.OneOf(
                [
                    albu.augmentations.transforms.Sharpen (alpha=(0.2, 0.5), lightness=(0.5, 1.0), always_apply=False, p=0.5),
                    albu.Blur(blur_limit=3, p=1),
                ],
                p=0.9,
            ),
            # albu.Normalize(mean=(0.2767,0.2767,0.2767),
            #        std=(0.0703,0.0703,0.0703), p=1),
        ]
    else:
        train_transform = []
            
    return albu.Compose(train_transform)

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
    # import matplotlib.pyplot as plt
    
    ROOT = '/media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/'
    DATA_DIR = os.path.join(ROOT, 'DATA', 'DEVELOPMENT')
    RESULTS_DIR =  os.path.join(ROOT, 'RESULTS', 'PYTHON')
    
    Experiment_name = 'UNET-480x480-v4'
    # Experiment_folder = os.makedirs(os.path.join(RESULTS_DIR,Experiment_name), exist_ok = True)
    
        
    x_train_dir = os.path.join(DATA_DIR, 'TRAINING', 'IMAGES-RESIZED')
    y_train_dir = os.path.join(DATA_DIR, 'TRAINING', 'MASKS-HYBRID-RESIZED-v1')
    
    x_valid_dir = os.path.join(DATA_DIR, 'VALIDATION', 'IMAGES-RESIZED')
    y_valid_dir = os.path.join(DATA_DIR, 'VALIDATION', 'MASKS-HYBRID-RESIZED-v1')
    
    # x_test_dir = os.path.join(DATA_DIR, 'DATA', 'DL_FSHD', 'Test', 'images')
    # y_test_dir = os.path.join(DATA_DIR, 'DATA', 'DL_FSHD', 'Test', 'masks')
        
    # Lets look at data we have
    
    dataset = Dataset(x_train_dir, y_train_dir, classes=['imt'])
    
    image, mask = dataset[4] # get some sample
    # image=image.transpose(1,2,0)
    
    visualize(
        # image=image.transpose(1,2,0), 
        image=image.squeeze(), 
        muscle_mask=mask.squeeze(),
    )
    
    
    #### Visualize resulted augmented images and masks
    
    augmented_dataset = Dataset(
        x_train_dir, 
        y_train_dir, 
        augmentation=get_training_augmentation(), 
        classes=['imt'],
    )
    
    # same image with different random transforms
    for i in range(3):
        image, mask = augmented_dataset[1]
        # visualize(image=image.transpose(1,2,0), mask=mask.squeeze(-1))
        visualize(image=image, mask=mask)
        
    from torch.utils.tensorboard import SummaryWriter
    import segmentation_models_pytorch as smp
    # import torchvision
    # from torchsummary import summary
    
    # def freeze_encoder(model):
    #     for child in model.encoder.children():
    #         for param in child.parameters():
    #             param.requires_grad = False
    #     return
    
    # def unfreeze(model):
    #     for child in model.children():
    #         for param in child.parameters():
    #             param.requires_grad = True
    #     return

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
    train_dataset = Dataset(
        x_train_dir, 
        y_train_dir, 
        augmentation=get_training_augmentation(), 
        preprocessing=get_preprocessing(preprocessing_fn),
        classes=CLASSES,
    )
    
    valid_dataset = Dataset(
        x_valid_dir, 
        y_valid_dir, 
        augmentation=get_validation_augmentation(), 
        preprocessing=get_preprocessing(preprocessing_fn),
        classes=CLASSES,
    )
    
    train_loader = DataLoader(train_dataset, batch_size=8, shuffle=True, num_workers=16, pin_memory=True, drop_last=True)
    valid_loader = DataLoader(valid_dataset, batch_size=8, shuffle=True, num_workers=16, pin_memory=True, drop_last=True)
    
    # cd /media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/CODE/Python/segmentation_models
    # tensorboard --logdir /media/francesco/FMZ_archive/PROJECT-CUBS-SPIE/CODE/Python/segmentation_models/runs
    
    # Dice/F1 score - https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient
    # IoU/Jaccard score - https://en.wikipedia.org/wiki/Jaccard_index

        
    # loss = smp.utils.losses.DiceLoss()
    loss = smp.utils.losses.DiceLoss()
    
    metrics = [
        smp.utils.metrics.IoU(threshold=0.5)
    ]
    
    optimizer = torch.optim.Adam([ 
        dict(params=model.parameters(), lr=0.0002),
    ])
    
    
    # create epoch runners 
    # it is a simple loop of iterating over dataloader`s samples
    train_epoch = smp.utils.train.TrainEpoch(
        model, 
        loss=loss, 
        metrics=metrics, 
        optimizer=optimizer,
        device=DEVICE,
        verbose=True,
    )
    
    valid_epoch = smp.utils.train.ValidEpoch(
        model, 
        loss=loss, 
        metrics=metrics, 
        device=DEVICE,
        verbose=True,
    )
    
    # train model for 40 epochs
    
    max_score = 0
    min_loss = 10^9
    
    writer = SummaryWriter()
    
    patience_lr = 0
    patience_train = 0
    num_epochs_decay = 5
    # a=next(iter(train_loader))
    
    for i in range(0, 100):
        
        print('\nEpoch: {}'.format(i))
        train_logs = train_epoch.run(train_loader)
        valid_logs = valid_epoch.run(valid_loader)
        
        # ...log the running loss
        writer.add_scalars('training and validation loss',
                        {'Train': train_logs[loss.__name__],
                          'Val':valid_logs[loss.__name__]},
                        i)

        if valid_logs[loss.__name__] < (min_loss - 0.05*min_loss):
            min_loss = valid_logs[loss.__name__]
            
            max_score = valid_logs['iou_score']
            checkpoints_path = os.path.join(RESULTS_DIR,Experiment_name,'checkpoints')
            
            if not os.path.isdir(checkpoints_path):
                os.mkdir(checkpoints_path)
                
            torch.save(model, os.path.join(checkpoints_path,'best_model.pth'))
            print('Model saved!')
            
            # patience_lr = 0
            patience_train = 0
        else:
            # patience_lr += 1
            patience_train += 1
            print('Patience {}\n'.format(patience_train))
        
        # if patience_lr == 10:
        #     optimizer.param_groups[0]['lr'] = optimizer.param_groups[0]['lr']*0.1
        #     patience_lr = 0
        #     print('Decrease decoder learning rate by 10^-1!')
            
        # Decay learning rate                    
        if i + 1 == num_epochs_decay*1:
            optimizer.param_groups[0]['lr'] = 0.0001
            print ('Decay learning rate to lr: {}.'.format(optimizer.param_groups[0]['lr']))
        elif i + 1 == num_epochs_decay*2:
            optimizer.param_groups[0]['lr'] = 0.00005
            print ('Decay learning rate to lr: {}.'.format(optimizer.param_groups[0]['lr']))
        elif i + 1 == num_epochs_decay*3:
            optimizer.param_groups[0]['lr'] = 0.00001
            print ('Decay learning rate to lr: {}.'.format(optimizer.param_groups[0]['lr']))                 
        elif i + 1 == num_epochs_decay*4:
            optimizer.param_groups[0]['lr'] = 0.000001
            print ('Decay learning rate to lr: {}.'.format(optimizer.param_groups[0]['lr']))  
            
        if patience_train == 10:
            print('Early Stopping!')
            break
            
        # Make dir for epoch valid res
        folder_name = 'Epoch_{}_Loss_{:.2f}_Dice_{:.2f}'.format(i+1, valid_logs[loss.__name__], valid_logs['iou_score'])
        folder_path = os.path.join(RESULTS_DIR,Experiment_name,'Validation_snapshot', folder_name)
        os.makedirs(folder_path)
        
        random_list = []
        
        for i in range(0, 10):
            # any random numbers from 0 to 1000
            random_list.append(random.randint(0, len(valid_dataset)))

        for j in tqdm(random_list):
            
            image, gt_mask = valid_dataset[j]
            # image, gt_mask = train_dataset[j]
            
            gt_mask = gt_mask.squeeze()
            
            x_tensor = torch.from_numpy(image).to(DEVICE).unsqueeze(0)
            pr_mask = model.predict(x_tensor)
            # pr_mask = torch.sigmoid(pr_mask)
            pr_mask = (pr_mask.squeeze().cpu().numpy().round())
                
            img = image[0,]
            stack = np.concatenate((img, pr_mask), axis=1)
            # plt.imshow(img)
            # plt.imshow(pr_mask)
            # plt.imshow(stack)
            
            cv2.imwrite(os.path.join(folder_path,os.listdir(x_valid_dir)[j]),stack*255)
                        # pr_mask*255)
    
    # load best saved checkpoint
    best_model = torch.load(os.path.join(checkpoints_path,'best_model.pth'))
        
    #### TRAIN SET
    
    # create test dataset
    # train_dataset = Dataset(
    #     x_train_dir, 
    #     y_train_dir, 
    #     augmentation=get_validation_augmentation(), 
    #     preprocessing=get_preprocessing(preprocessing_fn),
    #     classes=CLASSES,
    # )
    
    # train_dataloader = DataLoader(train_dataset)
    
    # evaluate model on test set
    # train_epoch = smp.utils.train.ValidEpoch(
    #     model=best_model,
    #     loss=loss,
    #     metrics=metrics,
    #     device=DEVICE,
    # )
    
    # logs = train_epoch.run(train_dataloader)

        
    # Train_output_folder = os.makedirs(os.path.join('./',Experiment_name,'Train_output'), exist_ok = True)

    # for i in tqdm(range(len(train_dataset))):
        
    #     image_vis = train_dataset[i][0].astype('uint8')
    #     image, gt_mask = train_dataset[i]
        
    #     gt_mask = gt_mask.squeeze()
        
    #     x_tensor = torch.from_numpy(image).to(DEVICE).unsqueeze(0)
    #     pr_mask = best_model.predict(x_tensor)
    #     pr_mask = (pr_mask.squeeze().cpu().numpy().round())
            
    #     cv2.imwrite(os.path.join('./',Experiment_name,'Train_output',os.listdir(x_train_dir)[i]),
    #                 pr_mask*255)
        
    #     #     visualize(
    #     #         image=image_vis, 
    #     #         ground_truth_mask=gt_mask, 
    #     #         predicted_mask=pr_mask
    #     #     )
        
    #### VAL SET
    
    # create test dataset
    val_dataset = Dataset(
        x_valid_dir, 
        y_valid_dir, 
        augmentation=get_validation_augmentation(), 
        preprocessing=get_preprocessing(preprocessing_fn),
        classes=CLASSES,
    )
    
    val_dataloader = DataLoader(val_dataset)
    
    # evaluate model on test set
    val_epoch = smp.utils.train.ValidEpoch(
        model=best_model,
        loss=loss,
        metrics=metrics,
        device=DEVICE,
    )
    
    logs = val_epoch.run(val_dataloader)

        
    Val_output_folder = os.makedirs(os.path.join(RESULTS_DIR,Experiment_name,'Val_output'), exist_ok = True)

    for i in tqdm(range(len(val_dataset))):
        
        image_vis = val_dataset[i][0].astype('uint8')
        image, gt_mask = val_dataset[i]
        
        gt_mask = gt_mask.squeeze()
        
        x_tensor = torch.from_numpy(image).to(DEVICE).unsqueeze(0)
        pr_mask = best_model.predict(x_tensor)
        pr_mask = (pr_mask.squeeze().cpu().numpy().round())
            
        cv2.imwrite(os.path.join(RESULTS_DIR,Experiment_name,'Val_output',os.listdir(x_valid_dir)[i]),
                    pr_mask*255)
        
        #     visualize(
        #         image=image_vis, 
        #         ground_truth_mask=gt_mask, 
        #         predicted_mask=pr_mask
        #     )        
        
        
        
        
#### inference parrt


        
    ### TEST SET
    
    # # create test dataset
    # test_dataset = Dataset(
    #     x_test_dir, 
    #     y_test_dir, 
    #     augmentation=get_validation_augmentation(), 
    #     preprocessing=get_preprocessing(preprocessing_fn),
    #     classes=CLASSES,
    # )
    
    # test_dataloader = DataLoader(test_dataset)
    
    # # evaluate model on test set
    # test_epoch = smp.utils.train.ValidEpoch(
    #     model=best_model,
    #     loss=loss,
    #     metrics=metrics,
    #     device=DEVICE,
    # )
    
    # logs_test = test_epoch.run(test_dataloader)

        
    # Test_output_folder = os.makedirs(os.path.join(RESULTS_DIR,Experiment_name,'Test_output'), exist_ok = True)
    
    
    # for i in tqdm(range(len(test_dataset))):
        
    #     image_vis = test_dataset[i][0].astype('uint8')
    #     image, gt_mask = test_dataset[i]
        
    #     gt_mask = gt_mask.squeeze()
        
    #     x_tensor = torch.from_numpy(image).to(DEVICE).unsqueeze(0)
    #     pr_mask = best_model.predict(x_tensor)
    #     pr_mask = (pr_mask.squeeze().cpu().numpy().round())
    #     # pr_mask = (pr_mask.squeeze().cpu().numpy())
            
    #     cv2.imwrite(os.path.join(RESULTS_DIR,Experiment_name,'Test_output',os.listdir(x_test_dir)[i]),
    #                 pr_mask*255)
        

        
        