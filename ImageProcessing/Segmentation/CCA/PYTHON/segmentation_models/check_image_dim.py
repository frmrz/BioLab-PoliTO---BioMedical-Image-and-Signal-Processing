# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 10:01:33 2021

@author: franc
"""

import os
import cv2
from tqdm import tqdm
import matplotlib as plt
import numpy as np

DATA_PATH = os.path.dirname(os.path.realpath(__file__))

img_path = os.path.join(DATA_PATH, 'data', 'MuscleCNN', 'train')


imlist = os.listdir(img_path)

to_mean = np.zeros((480,512))
mean_to = np.zeros((len(imlist),1))
std_to = np.zeros((len(imlist),1))
all_pixels = []
q=0

for i in imlist:
    img = cv2.imread(os.path.join(img_path,i),cv2.IMREAD_GRAYSCALE)   
    img_flat = img.flatten().tolist()
    to_mean = to_mean + img
    all_pixels.append(img_flat)
    mean_to[q,0] = np.mean(img)
    std_to[q,0] = np.std(img)
    q+=1
    
all_mean = np.mean(to_mean/len(imlist))/255
all_std =  np.std(to_mean/len(imlist))/255

mm = np.mean(mean_to)/255
ss = np.mean(std_to)/255

from itertools import chain

alpix = list(chain(*all_pixels))
alpixN = np.array(alpix)

mall = np.mean(alpixN)/255
sall = np.var(alpixN)/255

# res_path = os.path.join(DATA_PATH, 'MuscleCNN', 'grayset', 'testannot')

# # os.makedirs(res_path)

# imlist = os.listdir(img_path)

# for n, id_ in tqdm(enumerate(imlist), total=len(imlist)):
#     # Load images
#     img = cv2.imread(os.path.join(img_path,id_),cv2.IMREAD_GRAYSCALE)
#     cv2.imwrite(os.path.join(res_path,id_),img)

# folder = 'UnetPlusPlus_Resnet101_DiceLoss'

# group = 'Test_ouput'

# img_path = os.path.join(DATA_PATH, folder, group)
# res_path = os.path.join(DATA_PATH, folder, group)

# # os.makedirs(res_path,exist_ok=True)
# imlist = os.listdir(img_path)

# for n, id_ in tqdm(enumerate(imlist), total=len(imlist)):
#     # Load images
#     img = cv2.imread(os.path.join(img_path,id_),cv2.IMREAD_GRAYSCALE)
    
#     imagem = cv2.bitwise_not(img)
    
#     # img_big = cv2.resize(img, (512,480), interpolation=cv2.INTER_CUBIC)
#     # img_th  = cv2.threshold(img_big,127.5,255,cv2.THRESH_BINARY)
    
#     # plt.pyplot.imshow(imagem)
    
#     cv2.imwrite(os.path.join(res_path,id_),imagem)
    
# group = 'Train_ouput'

# img_path = os.path.join(DATA_PATH, folder, group)
# res_path = os.path.join(DATA_PATH, folder, group)

# # os.makedirs(res_path,exist_ok=True)
# imlist = os.listdir(img_path)

# for n, id_ in tqdm(enumerate(imlist), total=len(imlist)):
#     # Load images
#     img = cv2.imread(os.path.join(img_path,id_),cv2.IMREAD_GRAYSCALE)
    
#     imagem = cv2.bitwise_not(img)
    
#     # img_big = cv2.resize(img, (512,480), interpolation=cv2.INTER_CUBIC)
#     # img_th  = cv2.threshold(img_big,127.5,255,cv2.THRESH_BINARY)
    
#     # plt.pyplot.imshow(imagem)
    
#     cv2.imwrite(os.path.join(res_path,id_),imagem)
    
# group = 'Val_ouput'

# img_path = os.path.join(DATA_PATH, folder, group)
# res_path = os.path.join(DATA_PATH, folder, group)

# # os.makedirs(res_path,exist_ok=True)
# imlist = os.listdir(img_path)

# for n, id_ in tqdm(enumerate(imlist), total=len(imlist)):
#     # Load images
#     img = cv2.imread(os.path.join(img_path,id_),cv2.IMREAD_GRAYSCALE)
    
#     imagem = cv2.bitwise_not(img)
    
#     # img_big = cv2.resize(img, (512,480), interpolation=cv2.INTER_CUBIC)
#     # img_th  = cv2.threshold(img_big,127.5,255,cv2.THRESH_BINARY)
    
#     # plt.pyplot.imshow(imagem)
    
    # cv2.imwrite(os.path.join(res_path,id_),imagem)