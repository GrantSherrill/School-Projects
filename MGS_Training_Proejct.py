# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 15:25:35 2023

@author: GrantS
"""
#%%
#This file serves as a record of all of the code I used in my pproject, there is an updated one without the failed models and more comments for future use
#%%
import pandas as pd
import os 
import numpy as np
import random
import cv2

from PIL import Image
from skimage import data, io, filters
from skimage.filters import threshold_otsu
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg
import scipy.misc
from scipy import ndimage

import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras import backend as K

from tensorflow.python.ops import array_ops
from tensorflow.python.ops import math_ops

from tensorflow.keras import layers
from tensorflow.keras.models import Model
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.layers import Input, Conv2D, MaxPooling2D, Dropout, UpSampling2D
from tensorflow.keras.layers import Concatenate
from tensorflow.keras.optimizers import Adam
#%%Importing data

DataDir=r'C:\Users\GrantS\Documents\Data\jm161mix'
model_file = DataDir+'/MGSProejctBGSubResnetUpdated.hdf5'

FileBase='jm161mixglutre_'  #phase contrast files to be analyzed
Filter='P'         #Filter to make binary masks
Targetf='P'

FrameMin=1   #First frame to be analyzed
FrameMax=400    #Last frame to be analyzed
ChannelMin=1
ChannelMax=5
PosMin=1
PosMax=35
nImage=0
flip=17
for iPos in range(PosMin,PosMax+1):
    for iChannel in range(ChannelMin,ChannelMax+1):
        ChannelDir=str(DataDir)+'/Pos'+str(iPos)+'/Ch'+str(iChannel)+'/'
        LabelDir=str(DataDir)+'/Pos'+str(iPos)+'/ChBGSubT'+str(iChannel)+'/'
        for iFrame in range(FrameMin,FrameMax+1):
            LabelFileName=FileBase+'ChBGSub'+str(iChannel)+'_'+Filter+'_'+str(iFrame)+'.tif'
            TargetFileName=FileBase+'Ch'+str(iChannel)+'_'+Targetf+'_'+str(iFrame)+'.tif'
            TestFileName=FileBase+'Ch'+str(iChannel)+'_P_'+str(iFrame)+'.tif'
            #print(ChannelDir+BinFileName)
            if os.path.exists(LabelDir+LabelFileName) and os.path.exists(ChannelDir+TargetFileName):
                nImage+=1
                if nImage==1:
                    TestIm=mpimg.imread(LabelDir+LabelFileName)#Use the 1st image to determine the dimensions of the image
                    TestIm2=mpimg.imread(ChannelDir+TargetFileName)
                    ImSize=np.ndim(TestIm)
                    ImSize2=np.ndim(TestIm2)

#Load training images
c=44
LabelArray=np.zeros([nImage,236,32,1],dtype=float)
TargetArray=np.zeros([nImage,236,32,1],dtype=float)
LabelCrop=np.zeros([nImage,236-c,32,1],dtype=float)
TargetCrop=np.zeros([nImage,236-c,32,1],dtype=float)
#TestArray=np.zeros([nImage,TestIm.shape[0],TestIm.shape[1],1])

iImage=0
for iPos in range(PosMin,PosMax+1):
    for iChannel in range(ChannelMin,ChannelMax+1):
        ChannelDir=str(DataDir)+'/Pos'+str(iPos)+'/Ch'+str(iChannel)+'/'
        LabelDir=str(DataDir)+'/Pos'+str(iPos)+'/ChBGSubT'+str(iChannel)+'/'

        for iFrame in range(FrameMin,FrameMax+1):
            #BinFileName=FileBase+'Ch'+str(iChannel)+'_'+Filter+'Bin_'+str(iFrame)+'.tif'
            TargetFileName=FileBase+'Ch'+str(iChannel)+'_'+Targetf+'_'+str(iFrame)+'.tif'
            TestFileName=FileBase+'Ch'+str(iChannel)+'_P_'+str(iFrame)+'.tif'
            LabelFileName=FileBase+'ChBGSub'+str(iChannel)+'_'+Filter+'_'+str(iFrame)+'.tif'
            if os.path.exists(LabelDir+LabelFileName) and os.path.exists(ChannelDir+TargetFileName):
                LabelArray[iImage,:TestIm.shape[0],:TestIm.shape[1],0]=mpimg.imread(LabelDir+LabelFileName)

                TargetArray[iImage,:TestIm2.shape[0],:TestIm2.shape[1],0]=mpimg.imread(ChannelDir+TargetFileName)
                if iPos<=flip:
                    LabelCrop[iImage,:,:,:]=LabelArray[iImage,c:,:,:]
                    TargetCrop[iImage,:,:,:]=TargetArray[iImage,c:,:,:]
                else:
                    LabelCrop[iImage,:,:,:]=LabelArray[iImage,:-c,:,:]
                    TargetCrop[iImage,:,:,:]=TargetArray[iImage,:-c,:,:]
                #LabelArray=mpimg.imread(LabelDir+LabelFileName)

                #TargetArray=mpimg.imread(ChannelDir+TargetFileName)
                iImage+=1      
      
  #%%              
plt.figure()
plt.imshow(LabelCrop[1000,:,:,0],cmap='Greys')
plt.colorbar()
Testarray1=LabelArray[0,:,:,0]  
TargetDims=TargetArray.shape
#%%Renormalizing

nimages=LabelCrop.shape[0]
newimage=np.zeros([nImage,TargetCrop.shape[1],TargetCrop.shape[2],1],dtype='float64')
labels=np.zeros([nImage,TargetCrop.shape[1],TargetCrop.shape[2],1],dtype='float64')
for n in range(nimages):
        
    image=TargetCrop[n,:,:,0]
    lab=LabelCrop[n,:,:,0]

    m=np.min(image)
    mx=np.max(image)
        
    #m=np.int32(m)
    #mx=np.int32(mx)
    newimage[n,:,:,0]=(TargetCrop[n,:,:,0]-m)/mx
    
    
    
    lab=LabelCrop[n,:,:,0]
    m=np.min(lab)
    mx=np.max(lab)
    labels[n,:,:,0]=(LabelCrop[n,:,:,0]-m)/mx
        #TargetArray[n,:,:,0]= newimage[n,:,:,0]-newimage[0,:,:,0].min()
label=np.zeros([nImage,TargetCrop.shape[1],TargetCrop.shape[2],1],dtype='float64')
for i in range(iImage):
    label[i,:,:,0]=1-(labels[i,:,:,0])
                

          
test=label[800,:,:,0]
plt.figure()
plt.imshow(test,cmap='Greys')
plt.colorbar()
#%%
from sklearn.model_selection import train_test_split
x_train, x_te, y_train, y_te = train_test_split(newimage, label, test_size=0.2)
x_test, x_eval, y_test, y_eval = train_test_split(x_te, y_te, test_size=0.5)
#%%DOES NOT SEEM TO LEARN WELL AT ALL
#With this architecture I can see some features, but it definitely is not great
from tensorflow.keras import layers
model = keras.Sequential(
[
    layers.Dense(100,input_shape=(x_train.shape[1],x_train.shape[2],1),activation='relu'),
    layers.Flatten(),
    layers.Dense(200,activation='relu'),
    layers.Dense(400,activation='relu'),
    layers.Dense(800,activation='relu'),
    layers.Dense(400,activation='relu'),
    layers.Dense(200,activation='relu'),
    layers.Dense((x_train.shape[1]*x_train.shape[2]),activation='sigmoid'),
    layers.Reshape((x_train.shape[1],x_train.shape[2])),
])
model.summary()
model.compile(loss=keras.losses.mean_squared_error, optimizer='adam', metrics=['mse']) 
#%%
model = keras.Sequential(
    [
        layers.Conv2D(input_shape=(x_train.shape[1],x_train.shape[2],1), kernel_size=(3, 3), filters=2,
                      activation='relu', padding='same'),
        layers.MaxPooling2D(pool_size=(2, 2)),
        layers.Conv2D(20, (3, 3), activation='relu'),
        layers.Dropout(0.5),
        layers.MaxPooling2D(pool_size=(2, 2)),
        layers.Flatten(),
        layers.Dense(20*4*4, activation='relu'),
        layers.Dropout(0.5),
        layers.Dense(x_train.shape[1]*x_train.shape[2], activation='sigmoid'),
        layers.Reshape((x_train.shape[1],x_train.shape[2]))
    ])

model.summary()
model.compile(loss=keras.losses.mean_squared_error, optimizer='adam', metrics=['mse']) 
#%%
model = keras.Sequential(
    [
        layers.Conv2D(input_shape=(x_train.shape[1],x_train.shape[2],1), kernel_size=(3, 3), filters=64,
                      activation='relu', padding='same'),
        layers.Conv2D(64, (3, 3), activation='relu',padding='same'),
        layers.MaxPooling2D(pool_size=(2, 2)),
        layers.Conv2D(128, (3, 3), activation='relu',padding='same'),
        layers.Conv2D(128, (3, 3), activation='relu',padding='same'),
        layers.Dropout(0.1),
        layers.MaxPooling2D(pool_size=(2, 2)),
        layers.Conv2D(256, (3, 3), activation='relu',padding='same'),
        layers.Conv2D(256, (3, 3), activation='relu',padding='same'),
        layers.Dropout(0.1),
        layers.Conv2D(128, (3, 3), activation='relu',padding='same'),
        layers.Conv2D(128, (3, 3), activation='relu',padding='same'),
        layers.Conv2D(64, (3, 3), activation='relu',padding='same'),
        layers.Conv2D(64, (3, 3), activation='relu',padding='same'),
        layers.Flatten(),
        layers.Dense(1000,activation='relu'),
        layers.Dropout(0.1),
        layers.Dense(x_train.shape[1]*x_train.shape[2], activation='sigmoid'),
        layers.Reshape((x_train.shape[1],x_train.shape[2]))
    ])

model.summary()
model.compile(loss=keras.losses.mean_squared_error, optimizer='adam', metrics=['mse']) 


#%%smaller unet MODEL

def smallnet(input_size = (x_train.shape[0],x_train.shape[1],1), final_activation = 'sigmoid', output_classes = 1):
    '''
    Generic U-Net declaration. Only the input size, final activation layer and 
    dimensionality of the output vary between the different uses.
    '''
    
    inputs = Input(input_size,  name='true_input')	
	
    conv1 = Conv2D(16, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(inputs)
    conv1 = Conv2D(16, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv1)
    pool1 = MaxPooling2D(pool_size=(2, 2))(conv1)
    conv2 = Conv2D(32, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(pool1)
    conv2 = Conv2D(32, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv2)
    drop2 = Dropout(0.5)(conv2)
    pool2 = MaxPooling2D(pool_size=(2, 2))(drop2)

    conv3 = Conv2D(64, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(pool2)
    conv3 = Conv2D(64, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv3)
    drop3 = Dropout(0.5)(conv3)

    up4 = Conv2D(32, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(UpSampling2D(size = (2,2))(drop3))
    merge4 = Concatenate(axis = 3)([drop2,up4]) 
    conv4 = Conv2D(32, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(merge4)
    conv4 = Conv2D(32, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv4)


    up5 = Conv2D(16, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(UpSampling2D(size = (2,2))(conv4))
    merge5 = Concatenate(axis = 3)([conv1,up5])
    conv5 = Conv2D(16, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(merge5)
    conv5 = Conv2D(16, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv5)
    conv6 = Conv2D(output_classes, 1, activation = final_activation, name = 'true_output')(conv5)
    
    model = Model(inputs = inputs, outputs = conv6)
    
    
    return model
 #%%      
def unet(input_size = (x_train.shape[0],x_train.shape[1],1), final_activation = 'sigmoid', output_classes = 1):
    '''
    Generic U-Net declaration. Only the input size, final activation layer and 
    dimensionality of the output vary between the different uses.
    '''
    
    inputs = Input(input_size,  name='true_input')	
	
    conv1 = Conv2D(64, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(inputs)
    conv1 = Conv2D(64, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv1)
    pool1 = MaxPooling2D(pool_size=(2, 2))(conv1)
    conv2 = Conv2D(128, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(pool1)
    conv2 = Conv2D(128, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv2)
    pool2 = MaxPooling2D(pool_size=(2, 2))(conv2)
    conv3 = Conv2D(256, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(pool2)
    conv3 = Conv2D(256, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv3)
    pool3 = MaxPooling2D(pool_size=(2, 2))(conv3)
    conv4 = Conv2D(512, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(pool3)
    conv4 = Conv2D(512, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv4)
    drop4 = Dropout(0.5)(conv4)
    pool4 = MaxPooling2D(pool_size=(2, 2))(drop4)

    conv5 = Conv2D(1024, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(pool4)
    conv5 = Conv2D(1024, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv5)
    drop5 = Dropout(0.5)(conv5)

    up6 = Conv2D(512, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(UpSampling2D(size = (2,2))(drop5))
    merge6 = Concatenate(axis = 3)([drop4,up6]) 
    conv6 = Conv2D(512, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(merge6)
    conv6 = Conv2D(512, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv6)

    up7 = Conv2D(256, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(UpSampling2D(size = (2,2))(conv6))
    merge7 = Concatenate(axis = 3)([conv3,up7])
    conv7 = Conv2D(256, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(merge7)
    conv7 = Conv2D(256, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv7)

    up8 = Conv2D(128, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(UpSampling2D(size = (2,2))(conv7))
    merge8 = Concatenate(axis = 3)([conv2,up8])
    conv8 = Conv2D(128, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(merge8)
    conv8 = Conv2D(128, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv8)

    up9 = Conv2D(64, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(UpSampling2D(size = (2,2))(conv8))
    merge9 = Concatenate(axis = 3)([conv1,up9])
    conv9 = Conv2D(64, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(merge9)
    conv9 = Conv2D(64, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv9)
    conv10 = Conv2D(output_classes, 1, activation = final_activation, name = 'true_output')(conv9)
    
    model = Model(inputs = inputs, outputs = conv10)
    
    return model

#%%      
def resnet(input_size = (x_train.shape[0],x_train.shape[1],1), final_activation = 'sigmoid', output_classes = 1):
    '''
    Generic U-Net declaration. Only the input size, final activation layer and 
    dimensionality of the output vary between the different uses.
    '''
    
    inputs = Input(input_size,  name='true_input')	
	
    temp1=inputs
    conv1 = Conv2D(64, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(inputs)
    conv1 = Conv2D(64, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv1)
    out1= layers.Add()([temp1,conv1])
    out1=layers.ReLU()(out1)
    pool1 = MaxPooling2D(pool_size=(2, 2))(out1)
    
    conv2 = Conv2D(128, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(pool1)
    conv2 = Conv2D(128, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv2)
    temp2=Conv2D(128, 3, padding = 'same', kernel_initializer = 'he_normal')(pool1)
    out2= layers.Add()([temp2,conv2])
    out2=layers.ReLU()(out2)
    pool2 = MaxPooling2D(pool_size=(2, 2))(out2)
    conv3 = Conv2D(256, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(pool2)
    conv3 = Conv2D(256, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv3)
    temp3=Conv2D(256, 3, padding = 'same', kernel_initializer = 'he_normal')(pool2)
    out3= layers.Add()([temp3,conv3])
    out3=layers.ReLU()(out3)
    pool3 = MaxPooling2D(pool_size=(2, 2))(out3)
    conv4 = Conv2D(512, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(pool3)
    conv4 = Conv2D(512, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv4)
    temp4=Conv2D(512, 3, padding = 'same', kernel_initializer = 'he_normal')(pool3)
    out4= layers.Add()([temp4,conv4])
    out4=layers.ReLU()(out4)
    drop4 = Dropout(0.5)(out4)
    pool4 = MaxPooling2D(pool_size=(2, 2))(drop4)

    conv5 = Conv2D(1024, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(pool4)
    conv5 = Conv2D(1024, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv5)
    temp5=Conv2D(1024, 3, padding = 'same', kernel_initializer = 'he_normal')(pool4)
    out5= layers.Add()([temp5,conv5])
    out5=layers.ReLU()(out5)
    drop5 = Dropout(0.5)(out5)

    up6 = Conv2D(512, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(UpSampling2D(size = (2,2))(drop5))
    merge6 = Concatenate(axis = 3)([drop4,up6]) 
    conv6 = Conv2D(512, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(merge6)
    conv6 = Conv2D(512, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv6)
    temp6=Conv2D(512, 3, padding = 'same', kernel_initializer = 'he_normal')(merge6)
    out6= layers.Add()([temp6,conv6])
    out6=layers.ReLU()(out6)

    up7 = Conv2D(256, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(UpSampling2D(size = (2,2))(out6))
    merge7 = Concatenate(axis = 3)([out3,up7])
    conv7 = Conv2D(256, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(merge7)
    conv7 = Conv2D(256, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv7)
    temp7=Conv2D(256, 3, padding = 'same', kernel_initializer = 'he_normal')(merge7)
    out7= layers.Add()([temp7,conv7])
    out7=layers.ReLU()(out7)

    up8 = Conv2D(128, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(UpSampling2D(size = (2,2))(out7))
    merge8 = Concatenate(axis = 3)([out2,up8])
    conv8 = Conv2D(128, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(merge8)
    conv8 = Conv2D(128, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv8)
    temp8=Conv2D(128, 3, padding = 'same', kernel_initializer = 'he_normal')(merge8)
    out8= layers.Add()([temp8,conv8])
    out8=layers.ReLU()(out8)

    up9 = Conv2D(64, 2, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(UpSampling2D(size = (2,2))(out8))
    merge9 = Concatenate(axis = 3)([out1,up9])
    conv9 = Conv2D(64, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(merge9)
    conv9 = Conv2D(64, 3, activation = 'relu', padding = 'same', kernel_initializer = 'he_normal')(conv9)
    temp9=Conv2D(64, 3, padding = 'same', kernel_initializer = 'he_normal')(merge9)
    out9= layers.Add()([temp9,conv9])
    out9=layers.ReLU()(out9)
    conv10 = Conv2D(output_classes, 1, activation = 'sigmoid', name = 'true_output')(out9)
    #out=layers.Dense(512, activation='relu')(conv10)
    #out=layers.Dense(1,activation='sigmoid')(out)
    model = Model(inputs = inputs, outputs = conv10)
    
    return model
                                   
#%%
input_size=TargetCrop.shape[1:]
model = resnet(input_size)
model.compile(optimizer = Adam(learning_rate = 1e-4), loss = 'mse', metrics = ['mse'])
model.summary()
model_checkpoint = ModelCheckpoint(model_file, monitor='loss',verbose=1, save_best_only=True)
history=model.fit(x_train,y_train,validation_data=(x_test,y_test),batch_size=32, epochs=50 )

#%% Save Local Copy of Model
model.save(model_file)
#%%
#Make the test data and replace imagearray

PredictArray=model.predict(x_eval)
#score=model.evaluate(x_eval,y_eval)
#%%
# look into training history
fig,ax = plt.subplots(1,1, sharex=True, figsize=(5,5))

#score=model.evaluate(x_eval,y_eval)

# accuracy
ax.plot(history.history['mse'], color='b')
ax.plot(history.history['val_mse'], ls='--', color='r')
ax.set_ylabel('Model Accuracy')
ax.set_xlabel('Epoch')
ax.legend(['train', 'test'], loc='best')
#ax.text(0.5,0.95,f'{score[1]:.2f}',horizontalalignment='center',verticalalignment='top', 
 #                        transform=ax.transAxes)
#ax[0].set_ylim(top=1)

plt.savefig(DataDir+'/TrainingUpdated.pdf', format="pdf", bbox_inches="tight")
#%%
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

predictions=model.predict(x_eval)
#%%Return a number for how many cells are "close"
BinTest=[]
differences=[]
avg=[]
stds=[]
mse=[]
nmse=[]
for i in range(len(PredictArray)):
    x_temp=PredictArray[i,:,:].flatten()
    y_temp=y_eval[i,:,:].flatten()
    mistakes= np.where( (x_temp> y_temp+0.1)| (x_temp<y_temp-0.1))[0]
    num_mistakes=len(mistakes)
    print(num_mistakes)
    diffarray=y_eval[i,:,:]-PredictArray[i,:,:]
    #plt.imshow(diffarray,cmap='Greys')
    differences.append(mistakes)
    meanerr=np.sum((x_temp-y_temp)**2)/len(x_temp)
    nmse.append(np.sqrt(np.sum((x_temp-y_temp)**2)/len(x_temp))/(np.average(y_temp)))
    avdiff=np.average(diffarray)
    mxdiff=np.max(diffarray)
    stddiff=np.std(diffarray)
    avg.append(avdiff)
    stds.append(stddiff)
    mse.append(meanerr)
#%%histogram of differences and stuff
histbins=np.linspace(np.min(avg),np.max(avg),15)
hbins=np.linspace(np.min(mse),np.max(mse),15)
fig,ax=plt.subplots(3,1, figsize=(12,12))
ax[0].axvline(np.average(avg),lw=2,c='red')
ax[0].hist(avg,histbins,edgecolor='black',density='true')
ax[0].set(xlim=(np.min(avg),np.max(avg)))
#ax[0].text(np.average(avg)+0.008,len(avg)/15+6,np.round(np.average(avg),4))
ax[0].set_title('Histogram of Average Pixel Difference')
ax[1].set_title('Histogram of Standard deviations')
ax[1].hist(stds,edgecolor='black')
ax[1].set(xlim=(0,np.max(stds)))
ax[1].axvline(np.average(stds),lw=2,c='red')
#ax[1].text(np.average(stds)+0.003,len(avg)/9+6,np.round(np.average(stds),4))
ax[2].set_title('Histogram of Mean Square Errors')
ax[2].hist(mse,hbins,edgecolor='black')
ax[2].set(xlim=(0,np.max(mse)))
ax[2].axvline(np.average(mse),lw=2,c='red')
#ax[2].text(np.average(mse)+0.003,len(avg)/9+6,np.round(np.average(mse),4))
#plt.savefig(DataDir+'/ErrorAnalysis.pdf', format="pdf", bbox_inches="tight")

#%%
TestData=model.predict(x_train)
#%%
Testnmse=[]
for i in range(len(TestData)):
    x_temp=TestData[i,:,:].flatten()
    y_temp=y_train[i,:,:].flatten()    
    Testnmse.append(np.sqrt(np.sum((x_temp-y_temp)**2)/len(x_temp))/(np.average(y_temp)))


#%%
histbins=np.linspace(np.min(Testnmse),np.max(Testnmse),15)
hbins=np.linspace(np.min(nmse),np.max(nmse),15)
fig,ax=plt.subplots(2,1, figsize=(12,12))
ax[0].set_title('Histogram of Normalized Mean Square Errors for Training')
ax[0].hist(Testnmse,histbins,edgecolor='black')
ax[0].set(xlim=(0,np.max(Testnmse)))
ax[0].axvline(np.average(Testnmse),lw=2,c='red')
#ax[1].text(np.average(stds)+0.003,len(avg)/9+6,np.round(np.average(stds),4))
ax[1].set_title('Histogram of Normalized Mean Square Errors for Evaluation')
ax[1].hist(nmse,hbins,edgecolor='black')
ax[1].set(xlim=(0,np.max(nmse)))
ax[1].axvline(np.average(nmse),lw=2,c='red')
#ax[2].text(np.average(mse)+0.003,len(avg)/9+6,np.round(np.average(mse),4))
plt.savefig(DataDir+'/ErrorAnalysis2.pdf', format="pdf", bbox_inches="tight")

#%%Plotting predicted and real data
i=0
n=48
for i in range(5):
    fig,ax=plt.subplots(1,4, figsize=(6,10))
    ax[0].imshow(x_eval[i+n,:,:], cmap="Greys")
    ax[0].set_title('Input')
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[1].imshow(y_eval[i+n,:,:], cmap="Greys")
    ax[1].set_title('Target')
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[2].imshow(PredictArray[i+n,:,:,0], cmap="Greys")
    ax[2].set_title('Predicition')
    ax[2].set_xticks([])
    ax[2].set_yticks([])
    ax[3].imshow(y_eval[i+n,:,:,0]-PredictArray[i+n,:,:,0], cmap="Greys")
    ax[3].set_title('Subtraction')
    ax[3].set_xticks([])
    ax[3].set_yticks([])
    #plt.figure()
    #plt.imshow(BinArray[i+50*n,:,:,0], cmap="Greys")
    i+=1
plt.subplots_adjust(wspace=1, hspace=0)
vartest=y_train[0,:,:,0]
#plt.savefig(DataDir+'/Predictions2.pdf', format="pdf", bbox_inches="tight")
#%%Reading in Real Data
    
    
    
DataDir=r'C:\Users\GrantS\Documents\Data\jm219'
ModelDir=r'C:\Users\GrantS\Documents\Data\jm161mix'
model_file = ModelDir+'/MGSProejctBGSubResnetUpdated.hdf5'

FileBase='jm219-glythy_'  #phase contrast files to be analyzed
Filter='P'         #Filter to make binary masks
Targetf='P'

FrameMin=1   #First frame to be analyzed
FrameMax=200    #Last frame to be analyzed
ChannelMin=1
ChannelMax=5
PosMin=21
PosMax=22
nImage=0
flip=17
#NOTE FOR J57 IT starts flipping from 219
for iPos in range(PosMin,PosMax+1):
    for iChannel in range(ChannelMin,ChannelMax+1):
        ChannelDir=str(DataDir)+'/Pos'+str(iPos)+'/Ch'+str(iChannel)+'/'
        for iFrame in range(FrameMin,FrameMax+1):
            
            TargetFileName=FileBase+'Ch'+str(iChannel)+'_'+Targetf+'_'+str(iFrame)+'.tif'
            TestFileName=FileBase+'Ch'+str(iChannel)+'_P_'+str(iFrame)+'.tif'
           
            if os.path.exists(ChannelDir+TargetFileName):
                nImage+=1
                if nImage==1:
                    
                    TestIm2=mpimg.imread(ChannelDir+TargetFileName)
                    
                    ImSize2=np.ndim(TestIm2)

#Load training images

c=44
TargetArray=np.zeros([nImage,236,32,1],dtype=float)
TargetCrop=np.zeros([nImage,236-c,32,1],dtype=float)
#TestArray=np.zeros([nImage,TestIm.shape[0],TestIm.shape[1],1])

iImage=0
for iPos in range(PosMin,PosMax+1):
    for iChannel in range(ChannelMin,ChannelMax+1):
        ChannelDir=str(DataDir)+'/Pos'+str(iPos)+'/Ch'+str(iChannel)+'/'

        for iFrame in range(FrameMin,FrameMax+1):
            #BinFileName=FileBase+'Ch'+str(iChannel)+'_'+Filter+'Bin_'+str(iFrame)+'.tif'
            TargetFileName=FileBase+'Ch'+str(iChannel)+'_'+Targetf+'_'+str(iFrame)+'.tif'
            TestFileName=FileBase+'Ch'+str(iChannel)+'_P_'+str(iFrame)+'.tif'
            if os.path.exists(ChannelDir+TargetFileName):

                TargetArray[iImage,:TestIm2.shape[0],:TestIm2.shape[1],0]=mpimg.imread(ChannelDir+TargetFileName)
                if iPos<=flip:
                    TargetCrop[iImage,:,:,:]=TargetArray[iImage,c:,:,:]
                else:
                    TargetCrop[iImage,:,:,:]=TargetArray[iImage,:-c,:,:]
                #LabelArray=mpimg.imread(LabelDir+LabelFileName)

                #TargetArray=mpimg.imread(ChannelDir+TargetFileName)
                iImage+=1      
        

test=TargetCrop[0,:,:,0]
#%%

nimages=TargetCrop.shape[0]
newimage=np.zeros([nImage,TargetCrop.shape[1],TargetCrop.shape[2],1],dtype='float64')
for n in range(nimages):
        
    image=TargetCrop[n,:,:,0]


    m=np.min(image)
    mx=np.max(image)
        
    #m=np.int32(m)
    #mx=np.int32(mx)
    newimage[n,:,:,0]=(TargetCrop[n,:,:,0]-m)/mx


                
#%%
model=tf.keras.models.load_model(model_file)
Prediction=model.predict(newimage)
#%%
i=0
n=401
for i in range(1):
    fig,ax=plt.subplots(1,2, figsize=(4,8))
    ax[0].imshow(newimage[i+n,:,:], cmap="Greys")
    ax[0].set_title('Input')
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[1].imshow(Prediction[i+n,:,:], cmap="Greys")
    ax[1].set_title('Prediction')
    ax[1].set_xticks([])
    ax[1].set_yticks([])

    #plt.imshow(BinArray[i+50*n,:,:,0], cmap="Greys")
    i+=1
#plt.savefig(DataDir+'/UnseenData2.pdf', format="pdf", bbox_inches="tight")

plt.subplots_adjust(wspace=1, hspace=0)
vartest=y_train[0,:,:,0]