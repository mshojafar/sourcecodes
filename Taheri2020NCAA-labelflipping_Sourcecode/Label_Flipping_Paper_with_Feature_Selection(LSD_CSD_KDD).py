# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 14:25:17 2019

@author: Rahim
"""
#*****************************************************************import Library*****************************************************************************
from __future__ import print_function
from sklearn.feature_selection import SelectFromModel
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split 
from sklearn.metrics import confusion_matrix
from sklearn import model_selection
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestRegressor
from scipy.sparse import csr_matrix, vstack, hstack
from scipy.sparse import coo_matrix
from keras.preprocessing.text import one_hot
from sklearn import metrics
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.semi_supervised import LabelPropagation
from sklearn.semi_supervised import LabelSpreading
from sklearn.semi_supervised import label_propagation
from sklearn.metrics import roc_auc_score
from sklearn.metrics import f1_score
from sklearn.cluster import KMeans
import math
#import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation , Flatten
from sklearn.metrics import log_loss
from keras.optimizers import SGD
from keras.layers.normalization import BatchNormalization
from keras.layers.convolutional import UpSampling2D
from keras.layers.convolutional import Conv2D, MaxPooling2D, MaxPooling1D
from keras.layers.embeddings import Embedding
from scipy import sparse
import pandas as pd
import numpy as np
#import random
import sklearn
from sklearn.metrics.pairwise import manhattan_distances
from keras.models import Model
from keras.layers import  Conv1D, multiply, GlobalMaxPool1D, Input , Lambda
import time
import argparse
#import math
from numpy import *
import os.path as osp
import scipy.sparse as sp
import pickle
from sklearn.metrics import accuracy_score
from warnings import simplefilter
#*********************************************************************************************************************************
CLASS = 'class'
CLASS_BEN = 'B'
CLASS_MAL = 'M'
DATA = 'data'
#********************************************Functions that will be used in this program*****************************************
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-tables', nargs='*', dest='input_tables')

    args = parser.parse_args()

    return args
#*********************************************************************************************************************************
def read_table(table_file):
    
        table = dict()
        
        with open(table_file, 'rb') as handle:
            while True:
                   try:
                           table = pickle.load(handle)
                   except EOFError:
                           break
        
        f_set=set()
        
        for k,v in table.items():
             for feature in v[DATA]:
                f_set.add(feature)
               
        return table , f_set
#******************************************************************************
def relevant_features(data, response_vector, features):
    rel_features = list()
    ranked_index=list()
    
    model =RandomForestRegressor()
    rfe = RFE(model, 1)
    fit = rfe.fit(data, response_vector)
    old_features=features

    for i in fit.ranking_:
        if i<len(features):
              rel_features.append(features[i])
    ranked_index=[old_features.index(x) for x in rel_features if x in old_features]
       
    return rel_features ,ranked_index
#*********************************************************************************************************************************
def build_table(tables):
    full_table = dict()

    file_set = set()
    
    for table in tables:
        file_set.update(table.keys())
        for key, val in table.items():
            full_table[key] = val
              
    files = list(file_set)
    return full_table, files
#*********************************************************************************************************************************
def convert_to_matrix(table, features, files):
    mat = sp.lil.lil_matrix((len(files), len(features)), dtype=np.int8)

    print("Input Data Size =  ", mat.get_shape())
    # the response vector
   
    cl = [0]*len(files)
    
    for key, val in table.items():
        k = files.index(key)
    
        if val[CLASS] is CLASS_BEN:
            cl[k] = 1
       
        for v in val[DATA]:
            try:
                idx = features.index(v)
                mat[k, idx] = 1
            except Exception as e:
                print(e)
                pass              
        
    return mat, cl
#******************************************************************************
def delete_row_lil(mat, i):
    if not isinstance(mat, sp.lil.lil_matrix):
        raise ValueError("works only for LIL format -- use .tolil() first")
    mat.rows = np.delete(mat.rows, i)
    mat.data = np.delete(mat.data, i)
    mat._shape = (mat._shape[0] - 1, mat._shape[1])
#*****************************************************************Main Function*******************************************************
def main():
    simplefilter(action='ignore', category=FutureWarning)
    args = parse_args()
    tables = []
    f_set = set()
    #read the data
    for t_files in args.input_tables:
        table, features = read_table(t_files)
        f_set = f_set.union(features)
        tables.append(table)    
    #************************************build table from data and convert to matrix***************************************************
    full_table, files = build_table(tables)
    files.sort()
    features = list(f_set)
    features.sort()
    mat, cl = convert_to_matrix(full_table, features, files)       
    print("************************Doing feature Ranking on all of the Data*************************")
    r_features,ranked_index = relevant_features(mat, cl, features)
    original_selected=ranked_index[1:301]
    data = sparse.lil_matrix(sparse.csr_matrix(mat)[:,original_selected])
    
    #******************************************Split data to train , test and validation**********************************************
    seed = 10
    test_size = 0.2
    X_train, X_test, Y_train, Y_test= train_test_split(data, cl, test_size= test_size, random_state=seed)
    test_size = 0.25
    X_train, X_val, Y_train, Y_val= train_test_split(X_train, Y_train, test_size= test_size, random_state=seed) 
    #***********************************************************************************************************************************
    print("                                                                                                                               ")
    print("                                                                                                                               ")
    print("*********Semi-Supervised Deep Learning Based Approach Against Label Flipping Attack in Malware Detection System*****************")
    print("                                                                                                                               ")   
    
    X_train=sparse.csr_matrix(X_train)
    print("row_train,column_train=", X_train.get_shape())
    print("                                                                   ")
    X_val=sparse.csr_matrix(X_val)
    row_val,column_val=X_val.get_shape()
    print("row_val,column_val=",X_val.get_shape())  
    print("                                                                   ")
    X_test=sparse.csr_matrix(X_test)
    row_test,column_test=X_test.get_shape()
    print("row_test,column_test=",X_test.get_shape()) 
    print("                                                                   ")
    print("********************************************************************")
    #**************************************************Model Definition*****************************************************************
    X_train_NoAttack=X_train.copy()
    Y_train_NoAttack=Y_train[:]
    
    X_val_NoAttack=X_val.copy()
    Y_val_NoAttack=Y_val[:]
    
    row_train_NoAttack,column_train_NoAttack=X_train_NoAttack.get_shape()
    model_main = Sequential()
    model_main.add(Embedding(row_train_NoAttack, 8, input_length=column_train_NoAttack))
    model_main.add(Conv1D(16,2, strides=2, padding='same'))
    model_main.add(MaxPooling1D(pool_size = (4), strides=(2)))
    model_main.add(Conv1D(32,2, strides=2, padding='same'))
    model_main.add(MaxPooling1D(pool_size = (4), strides=(2)))
    model_main.add(Conv1D(64,2, strides=2, padding='same'))
    model_main.add(Flatten())
    model_main.add(Dense(1, activation='sigmoid'))
    model_main.compile(optimizer='adam', loss='binary_crossentropy', metrics=['acc'])
    model_main.fit(X_train_NoAttack, Y_train_NoAttack, epochs=200, verbose=0)
   
    Y_CNN_NoAttack=model_main.predict(X_test, verbose=0)
    Y_predict_NoAttack=[0]*len(Y_CNN_NoAttack)
    
    for i in range(len(Y_CNN_NoAttack)):
        if Y_CNN_NoAttack[i]<0.5:
              Y_CNN_NoAttack[i]=0
        else:
              Y_CNN_NoAttack[i]=1

    for i in range(len(Y_CNN_NoAttack)):
        Y_predict_NoAttack[i]= int(Y_CNN_NoAttack[i])    
    #*****************************************************Result of Model without attack on X_test*****************************************
    print("********************************Result of Model without attack******************************************************************")
    loss, accuracy = model_main.evaluate(X_train_NoAttack, Y_train_NoAttack, verbose=2)
    print('Accuracy for Train set: %f' % (accuracy*100))
    print('Loss for Train set: %f' % (loss))
    print("                                                                   ")
    
    loss, accuracy = model_main.evaluate(X_val_NoAttack, Y_val_NoAttack, verbose=2)
    print('Accuracy for Validation set: %f' % (accuracy*100))
    print('Loss for Train Validation set: %f' % (loss))
    print("                                                                   ")
    
    loss, accuracy = model_main.evaluate(X_test, Y_test, verbose=2)
    print('Accuracy for Test set: %f' % (accuracy*100))
    print('Loss for Test set:: %f' % (loss))
    print("                                                                   ")

    TN_NoAttack, FP_NoAttack, FN_NoAttack, TP_NoAttack = confusion_matrix(Y_test,  Y_predict_NoAttack).ravel()
    print("TN_NoAttack=",TN_NoAttack)
    print("FP_NoAttack=",FP_NoAttack)
    print("FN_NoAttack=",FN_NoAttack)
    print("TP_NoAttack=",TP_NoAttack)
    print("                                                                   ")

    if (FP_NoAttack+TN_NoAttack)>0:
         FPR_NoAttack=FP_NoAttack/(FP_NoAttack+TN_NoAttack)
         print("The FPR_NoAttack result=", FPR_NoAttack)
                
    if (FP_NoAttack+TN_NoAttack)>0:
         TPR_NoAttack=TP_NoAttack/(TP_NoAttack+FN_NoAttack)
         print("The TPR_NoAttack result=", TPR_NoAttack)
                
    if (TN_NoAttack+FP_NoAttack)>0:
        TNR_NoAttack=TN_NoAttack/(TN_NoAttack+FP_NoAttack)
        print("The TNR_NoAttack result=", TNR_NoAttack)
                
    if (FN_NoAttack+TP_NoAttack)>0:
        FNR_NoAttack=FN_NoAttack/(FN_NoAttack+TP_NoAttack)
        print("The FNR_NoAttack result=", FNR_NoAttack)
                
    if ((TN_NoAttack/(TN_NoAttack+FP_NoAttack))+(TP_NoAttack/(TP_NoAttack+FP_NoAttack)))>0:
        AUC_NoAttack=1/(2*((TN_NoAttack/(TN_NoAttack+FP_NoAttack))+(TP_NoAttack/(TP_NoAttack+FP_NoAttack))))
        print("The AUC_NoAttack result=", AUC_NoAttack)
                
    if  (TP_NoAttack+TN_NoAttack+FP_NoAttack+FN_NoAttack)>0:   
        ACC_NoAttack=(TP_NoAttack+TN_NoAttack)/(TP_NoAttack+TN_NoAttack+FP_NoAttack+FN_NoAttack)
        print("The ACC_NoAttack result=", ACC_NoAttack)
                
    if  ((TP_NoAttack+FP_NoAttack)*(TP_NoAttack+FN_NoAttack)*(TN_NoAttack+FP_NoAttack)*(TN_NoAttack+FN_NoAttack))>0:
        MCC_NoAttack=(TP_NoAttack*TN_NoAttack-FP_NoAttack*FN_NoAttack)/math.sqrt((TP_NoAttack+FP_NoAttack)*(TP_NoAttack+FN_NoAttack)*(TN_NoAttack+FP_NoAttack)*(TN_NoAttack+FN_NoAttack))
        print("The Matthews correlation coefficient result=", MCC_NoAttack)
    print("                                                                                                                               ")
    print("*****************************************************End of Without Attack part************************************************")
    print("                                                                                                                               ")
    print("                                                                                                                               ")
    print("                                                                                                                               ")
    print("*****************************************************Label Flipping Attack*****************************************************")
    print("                                                                                                                               ")
    #************************** 
    # finding Malware of Train data
    malware_train= sparse.lil_matrix(X_train)
    cl_malware=list()
    z_m=0    
    count_m=0
    for i, j in enumerate(Y_train):
         if j == 1:
            delete_row_lil(malware_train, i-count_m)
            count_m=count_m+1
         else:
            cl_malware.insert(z_m, 1)
            z_m=z_m+1 
    #***************************
    #Finding Benign of Train data
    cl_X_train=list(Y_train) 
    benign_train=sparse.lil_matrix(X_train)
    z_b=0    
    count_b=0
    cl_benign=list()
    for i, j in enumerate(cl_X_train):
        if j == 0:
            delete_row_lil(benign_train, i-count_b)
            count_b=count_b+1
        else:
            cl_benign.insert(z_b, 1)
            z_b=z_b+1
    print("***********Size of Each Data Part:**********")        
    print("malware_train=", malware_train.get_shape())
    print("benign_train=", benign_train.get_shape())
    #***************************************************
    row_malware_train,column_malware_train=malware_train.get_shape()
    #Number_of_flipped_label=int(row_malware_train)
    
    X_train_LFA=X_train.copy()
    Y_train_LFA=Y_train[:]
    
    row_train_LFA,column_train_LFA=X_train_LFA.get_shape()
    clusterer = KMeans(n_clusters=2, random_state=10)
    X=X_train_LFA.toarray()
    t0=time.time()
    cluster_labels = clusterer.fit_predict(X)
    sample_silhouette_values = silhouette_samples(X, cluster_labels)
    #print("sample_silhouette_values=",sample_silhouette_values)
    
    flipped_Y_train=list(Y_train_LFA)
    counter=0
    for new_index in range(row_train_LFA): 
        if (sample_silhouette_values[new_index]<0.1):                           #and (flipped_Y_train[new_index]==0)
             flipped_Y_train[new_index]=abs(flipped_Y_train[new_index]-1)     #flipped_Y_train[new_index]=1
             counter=counter+1

    print("Flipped  counter=", counter)         
    t1=time.time()
    print("Time for Label Flipping Attack =",t1-t0)
    print("                                                                   ")
    
     #**************************************************************************
    model_main_LFA_Final = Sequential()
    model_main_LFA_Final.add(Embedding(row_train_LFA, 8, input_length=column_train_LFA))
    model_main_LFA_Final.add(Conv1D(16,2, strides=2, padding='same'))
    model_main_LFA_Final.add(MaxPooling1D(pool_size = (4), strides=(2)))
    model_main_LFA_Final.add(Conv1D(32,2, strides=2, padding='same'))
    model_main_LFA_Final.add(MaxPooling1D(pool_size = (4), strides=(2)))
    model_main_LFA_Final.add(Conv1D(64,2, strides=2, padding='same'))
    model_main_LFA_Final.add(Flatten())
    model_main_LFA_Final.add(Dense(1, activation='sigmoid'))
    model_main_LFA_Final.compile(optimizer='adam', loss='binary_crossentropy', metrics=['acc'])
    model_main_LFA_Final.fit(X_train_LFA, flipped_Y_train, epochs=200, verbose=0)
    
   
    Y_predict_LFA=model_main_LFA_Final.predict(X_test, verbose=0)
    Y_predict_LFA_Final=[0]*len(Y_predict_LFA)
    
    for i in range(len(Y_predict_LFA)):
        if Y_predict_LFA[i]<0.5:
              Y_predict_LFA[i]=0
        else:
              Y_predict_LFA[i]=1

    for i in range(len(Y_predict_LFA)):
        Y_predict_LFA_Final[i]= int(Y_predict_LFA[i])    
    #*****************************************************Result of Model with  LFA ******************************************************
    print("********************************Result of Model with LFA attack **************************************************************")
    print("                                                                                                                              ")
    loss, accuracy = model_main_LFA_Final.evaluate(X_train_LFA, flipped_Y_train, verbose=2)
    print('Accuracy for Train set: %f' % (accuracy*100))
    print('Loss for Train set: %f' % (loss))
    print("                                                                   ")
   
    loss, accuracy = model_main_LFA_Final.evaluate(X_test, Y_test, verbose=2)
    print('Accuracy for Test set: %f' % (accuracy*100))
    print('Loss for Test set:: %f' % (loss))
    print("                                                                   ")

    TN_LFA, FP_LFA, FN_LFA, TP_LFA = confusion_matrix(Y_test,  Y_predict_LFA_Final).ravel()
    print("TN_LFA=",TN_LFA)
    print("FP_LFA=",FP_LFA)
    print("FN_LFA=",FN_LFA)
    print("TP_LFA=",TP_LFA)
    print("                                                                   ")

    if (FP_LFA+TN_LFA)>0:
         FPR_LFA=FP_LFA/(FP_LFA+TN_LFA)
         print("The FPR_LFA result=", FPR_LFA)
                
    if (FP_LFA+TN_LFA)>0:
         TPR_LFA=TP_LFA/(TP_LFA+FN_LFA)
         print("The TPR_LFA result=", TPR_LFA)
                
    if (TN_LFA+FP_LFA)>0:
        TNR_LFA=TN_LFA/(TN_LFA+FP_LFA)
        print("The TNR_LFA result=", TNR_LFA)
                
    if (FN_LFA+TP_LFA)>0:
        FNR_LFA=FN_LFA/(FN_LFA+TP_LFA)
        print("The FNR_LFA result=", FNR_LFA)
                
    if ((TN_LFA/(TN_LFA+FP_LFA))+(TP_LFA/(TP_LFA+FP_LFA)))>0:
        AUC_LFA=1/(2*((TN_LFA/(TN_LFA+FP_LFA))+(TP_LFA/(TP_LFA+FP_LFA))))
        print("The AUC_LFA result=", AUC_LFA)
                
    if  (TP_LFA+TN_LFA+FP_LFA+FN_LFA)>0:   
        ACC_LFA=(TP_LFA+TN_LFA)/(TP_LFA+TN_LFA+FP_LFA+FN_LFA)
        print("The ACC_LFAk result=", ACC_LFA)
                
    if  ((TP_LFA+FP_LFA)*(TP_LFA+FN_LFA)*(TN_LFA+FP_LFA)*(TN_LFA+FN_LFA))>0:
        MCC_LFA=(TP_LFA*TN_LFA-FP_LFA*FN_LFA)/math.sqrt((TP_LFA+FP_LFA)*(TP_LFA+FN_LFA)*(TN_LFA+FP_LFA)*(TN_LFA+FN_LFA))
        print("The Matthews correlation coefficient result=", MCC_LFA)   
    print("                                                                                                                               ")
    print("************************************************End of Label Flipping Attack part**********************************************")
    print("                                                                                                                               ")
    print("                                                                                                                               ")
    print("                                                                                                                               ")
    print("*****************************************************KNN Based Semi-Supervised Defense(KSD)************************************")
    print("                                                                                                                               ")

    X_train_KNN=X_train.copy()
    Y_train_KNN=flipped_Y_train[:]
    
    X_val_KNN=X_val.copy()
    Y_val_KNN=Y_val[:]
    
    row_train_KNN,column_train_KNN=X_train_KNN.get_shape()
        
    Number_of_flipped_label=int(row_train_KNN/50)
    Y_train_corrected_By_KNN=list(Y_train_KNN)

    c=0
    m=0
    t2=time.time()
    
    for i in list(range(Number_of_flipped_label)):
         row_KNN=X_train_KNN.getrow(i)
         distances = sklearn.metrics.pairwise.manhattan_distances(row_KNN,X_val_KNN)
         indices = distances.argsort()[:10]
         d=indices[0]
         a=d[0:10]
        
         F=0
         for j in range(len(a)):
                 t=a[j]
                 F=F+Y_val_KNN[t]
         fraction=F/10
         if fraction>=0.5:
             Y_train_corrected_By_KNN[i]=1
             m=m+1
         else: 
             Y_train_corrected_By_KNN[i]=0
             c=c+1
    Y_train_corrected_By_KNN_Final=np.array(Y_train_corrected_By_KNN)  
    t3=time.time()
    print("Time for KNN Based Semi-Supervised Defense(KSD) =",t3-t2)
    print("                                                                   ")

    model_main_KNN = Sequential()
    model_main_KNN.add(Embedding(row_train_NoAttack, 8, input_length=column_train_NoAttack))
    model_main_KNN.add(Conv1D(16,2, strides=2, padding='same'))
    model_main_KNN.add(MaxPooling1D(pool_size = (4), strides=(2)))
    model_main_KNN.add(Conv1D(32,2, strides=2, padding='same'))
    model_main_KNN.add(MaxPooling1D(pool_size = (4), strides=(2)))
    model_main_KNN.add(Conv1D(64,2, strides=2, padding='same'))
    model_main_KNN.add(Flatten())
    model_main_KNN.add(Dense(1, activation='sigmoid'))
    model_main_KNN.compile(optimizer='adam', loss='binary_crossentropy', metrics=['acc'])
    model_main_KNN.fit(X_train_KNN,Y_train_corrected_By_KNN_Final, epochs=20, batch_size=32, verbose=0)
    Y_predict_KNN=model_main_KNN.predict(X_test, verbose=0)
    
    Y_predict_KNN_Final=[0]*len(Y_predict_KNN)
    for i in range(len(Y_predict_KNN)):
        if Y_predict_KNN[i]<0.5:
              Y_predict_KNN[i]=0
        else:
              Y_predict_KNN[i]=1
    
    for i in range(len(Y_predict_KNN)):
        Y_predict_KNN_Final[i]= int(Y_predict_KNN[i])
    #*****************************************************Result of Model After KNN Based Defense*****************************************
    print("************************Result After KNN_Based Defense************************************************************************")
    print("                                                                                                                               ")

    loss, accuracy = model_main_KNN.evaluate(X_train_KNN, Y_train_KNN, verbose=0)
    print('Accuracy for Train set: %f' % (accuracy*100))
    print('Loss for Train set: %f' % (loss))
    print("                                                                   ")

    loss, accuracy = model_main_KNN.evaluate(X_test, Y_test, batch_size=32, verbose=0)
    print('Accuracy After KNN-Based Defense: %f' % (accuracy*100))
    print('Loss After KNN-Based Defense: %f' % (loss))
    print("                                                                   ")

    TN_KNN, FP_KNN, FN_KNN, TP_KNN = confusion_matrix(Y_test,   Y_predict_KNN_Final).ravel()
    print("TN_KNN=",TN_KNN)
    print("FP_KNN=",FP_KNN)
    print("FN_KNN=",FN_KNN)
    print("TP_KNN=",TP_KNN)
    print("                                                                   ")

    if (FP_KNN+TN_KNN)>0:
         FPR_KNN=FP_KNN/(FP_KNN+TN_KNN)
         print("The FPR_KNN result=", FPR_KNN)
                
    if (FP_KNN+TN_KNN)>0:
         TPR_KNN=TP_KNN/(TP_KNN+FN_KNN)
         print("The TPR_KNN result=", TPR_KNN)
                
    if (TN_KNN+FP_KNN)>0:
        TNR_KNN=TN_KNN/(TN_KNN+FP_KNN)
        print("The TNR_KNN result=", TNR_KNN)
                
    if (FN_KNN+TP_KNN)>0:
        FNR_KNN=FN_KNN/(FN_KNN+TP_KNN)
        print("The FNR_KNN result=", FNR_KNN)
                
    if ((TN_KNN/(TN_KNN+FP_KNN))+(TP_KNN/(TP_KNN+FP_KNN)))>0:
        AUC_KNN=1/(2*((TN_KNN/(TN_KNN+FP_KNN))+(TP_KNN/(TP_KNN+FP_KNN))))
        print("The AUC_KNN result=", AUC_KNN)
                
    if  (TP_KNN+TN_KNN+FP_KNN+FN_KNN)>0:   
        ACC_KNN=(TP_KNN+TN_KNN)/(TP_KNN+TN_KNN+FP_KNN+FN_KNN)
        print("The ACC_KNN result=", ACC_KNN)
                
    if  ((TP_KNN+FP_KNN)*(TP_KNN+FN_KNN)*(TN_KNN+FP_KNN)*(TN_KNN+FN_KNN))>0:
        MCC_KNN=(TP_KNN*TN_KNN-FP_KNN*FN_KNN)/math.sqrt((TP_KNN+FP_KNN)*(TP_KNN+FN_KNN)*(TN_KNN+FP_KNN)*(TN_KNN+FN_KNN))
        print("The Matthews correlation coefficient result=", MCC_KNN)
    print("                                                                                                                               ")
    print("************************************************End of KNN Based Semi-Supervised Defense(KSD) part*****************************")
    print("                                                                                                                               ")
    print("                                                                                                                               ")
    print("                                                                                                                               ")    
    print("*****************************************************Label Based Semi-supervised Defense(LSD)**********************************")
    print("                                                                                                                               ")
    #***********************label Propagation and Label Spreading for Using in Label Based Semi-supervised Defense(LSD) *******************
    X_train_LSD=X_train.copy()
    Y_train_LSD=flipped_Y_train[:]
    
    X_val_LSD=X_val.copy()
    Y_val_LSD=Y_val[:]
    row_val_LSD,column_val_LSD=X_val_LSD.get_shape()
    row_train_LSD,column_train_LSD=X_train_LSD.get_shape()
    
    t4=time.time()
    
    labels = np.full(row_train_LSD, -1)
    for i in range(row_val_LSD):
        labels[i] = Y_val_LSD[i]

    X=X_train_LSD.toarray()
    label_spread = label_propagation.LabelSpreading(kernel='knn', alpha=0.8)
    label_propa=label_propagation.LabelPropagation(kernel='knn', gamma=20, n_neighbors=7, max_iter=1000, tol=0.001, n_jobs=None)
    label_spread.fit(X, labels)
    label_propa.fit(X, labels)
    output_labels_spread = label_spread.transduction_
    output_labels_propa = label_propa.transduction_
    #*******************Convolutional Neural Network for Using in Label Based Semi-supervised Defense(LSD) ******************************
    CNN_model_for_LSD = Sequential()
    CNN_model_for_LSD.add(Embedding(row_train_LSD, 8, input_length=column_train_LSD))
    CNN_model_for_LSD.add(Conv1D(16,2, strides=2, padding='same'))
    CNN_model_for_LSD.add(MaxPooling1D(pool_size = (4), strides=(2)))
    CNN_model_for_LSD.add(Conv1D(32,2, strides=2, padding='same'))
    CNN_model_for_LSD.add(MaxPooling1D(pool_size = (4), strides=(2)))
    CNN_model_for_LSD.add(Conv1D(64,2, strides=2, padding='same'))
    CNN_model_for_LSD.add(Flatten())
    
    CNN_model_for_LSD.add(Dense(1, activation='sigmoid'))
    CNN_model_for_LSD.compile(optimizer='adam', loss='binary_crossentropy', metrics=['acc'])
    CNN_model_for_LSD.fit(X_train_LSD, Y_train_LSD, epochs=200, verbose=0)
    
    Y_predict_CNN_for_LSD=CNN_model_for_LSD.predict(X_train_LSD, verbose=0)
    
    Y_predict_CNN_LSD_Final=[0]*len(Y_predict_CNN_for_LSD)
    for i in range(len(Y_predict_CNN_for_LSD)):
        if Y_predict_CNN_for_LSD[i]<0.5:
              Y_predict_CNN_for_LSD[i]=0
        else:
               Y_predict_CNN_for_LSD[i]=1
    
    for i in range(len(Y_predict_CNN_for_LSD)):
        Y_predict_CNN_LSD_Final[i]= int(Y_predict_CNN_for_LSD[i])
    #*******************************************Voting Between CNN , label Propagation and Label Spreading**************************     
    Y_predict_LSD_Final=[0]*len(Y_train)
    for i in range(len(Y_train)):
        c=Y_train_LSD[i]+Y_predict_CNN_LSD_Final[i]+output_labels_propa[i]+output_labels_spread[i]
        if 2<=c:
            Y_predict_LSD_Final[i]=1
        else:
            Y_predict_LSD_Final[i]=0
    t5=time.time()
    print("Time for Label Based Semi-supervised Defense =",t5-t4)
    print("                                                                                                                               ")
    #*********************************************************************************************************************************
    model_main_LSD = Sequential()
    model_main_LSD.add(Embedding(row_train_LSD, 8, input_length=column_train_LSD))
    model_main_LSD.add(Conv1D(16,2, strides=2, padding='same'))
    model_main_LSD.add(MaxPooling1D(pool_size = (4), strides=(2)))
    model_main_LSD.add(Conv1D(32,2, strides=2, padding='same'))
    model_main_LSD.add(MaxPooling1D(pool_size = (4), strides=(2)))
    model_main_LSD.add(Conv1D(64,2, strides=2, padding='same'))
    model_main_LSD.add(Flatten())
    model_main_LSD.add(Dense(1, activation='sigmoid'))
    model_main_LSD.compile(optimizer='adam', loss='binary_crossentropy', metrics=['acc'])
    model_main_LSD.fit(X_train_LSD, Y_predict_LSD_Final, epochs=200, verbose=0)
        
    Y_predict_LSD_Defense=model_main_LSD.predict(X_test, verbose=0)
    Y_predict_LSD_Defense_Final=[0]*len(Y_predict_LSD_Defense)
    
    for i in range(len(Y_predict_LSD_Defense)):
        if Y_predict_LSD_Defense[i]<0.5:
              Y_predict_LSD_Defense[i]=0
        else:
               Y_predict_LSD_Defense[i]=1
    
    for i in range(len(Y_predict_LSD_Defense)):
        Y_predict_LSD_Defense_Final[i]= int(Y_predict_LSD_Defense[i])  
    #**************************************Result of Model after Label Based Semi-supervised Defense(LSD)**********************************
    print("************************Result of Model after Label Based Semi-supervised Defense(LSD)*****************************************")
    print("                                                                                                                               ")
    loss, accuracy = model_main.evaluate(X_train, Y_predict_LSD_Final, verbose=2)
    print('Accuracy for Train set: %f' % (accuracy*100))
    print('Loss for Train set: %f' % (loss))
    print("                                                                   ")

    loss, accuracy = model_main.evaluate(X_test, Y_test, verbose=2)
    print('Accuracy for Test set: %f' % (accuracy*100))
    print('Loss for Test set:: %f' % (loss))
    print("                                                                   ")

    TN_LSD, FP_LSD, FN_LSD, TP_LSD = confusion_matrix(Y_test,  Y_predict_LSD_Defense_Final).ravel()
    print("TN_LSD=",TN_LSD)
    print("FP_LSD=",FP_LSD)
    print("FN_LSD=",FN_LSD)
    print("TP_LSD=",TP_LSD)
    print("                                                                   ")

    if (FP_LSD+TN_LSD)>0:
         FPR_LSD=FP_LSD/(FP_LSD+TN_LSD)
         print("The FPR_LSD result=", FPR_LSD)
                
    if (FP_LSD+TN_LSD)>0:
         TPR_LSD=TP_LSD/(TP_LSD+FN_LSD)
         print("The TPR_LSD result=", TPR_LSD)
                
    if (TN_LSD+FP_LSD)>0:
        TNR_LSD=TN_LSD/(TN_LSD+FP_LSD)
        print("The TNR_LSD result=", TNR_LSD)
                
    if (FN_LSD+TP_LSD)>0:
        FNR_LSD=FN_LSD/(FN_LSD+TP_LSD)
        print("The FNR_LSD result=", FNR_LSD)
                
    if ((TN_LSD/(TN_LSD+FP_LSD))+(TP_LSD/(TP_LSD+FP_LSD)))>0:
        AUC_LSD=1/(2*((TN_LSD/(TN_LSD+FP_LSD))+(TP_LSD/(TP_LSD+FP_LSD))))
        print("The AUC result=", AUC_LSD)
                
    if  (TP_LSD+TN_LSD+FP_LSD+FN_LSD)>0:   
        ACC_LSD=(TP_LSD+TN_LSD)/(TP_LSD+TN_LSD+FP_LSD+FN_LSD)
        print("The ACC result=", ACC_LSD)
                
    if  ((TP_LSD+FP_LSD)*(TP_LSD+FN_LSD)*(TN_LSD+FP_LSD)*(TN_LSD+FN_LSD))>0:
        MCC_LSD=(TP_LSD*TN_LSD-FP_LSD*FN_LSD)/math.sqrt((TP_LSD+FP_LSD)*(TP_LSD+FN_LSD)*(TN_LSD+FP_LSD)*(TN_LSD+FN_LSD))
        print("The Matthews correlation coefficient result=", MCC_LSD) 
    print("                                                                                                                               ")
    print("*****************************************************End of Label Based Semi-supervised Defense(LSD)***************************")   
    print("                                                                                                                               ")
    print("                                                                                                                               ")
    print("                                                                                                                               ")
    print("*****************************************************Clustering Based Semi-supervised Defense(CSD)*****************************")
    print("                                                                                                                               ")

    X_train_CSD=X_train.copy()
    Y_train_CSD=flipped_Y_train[:]
    
    X_val_CSD=X_val.copy()
    Y_val_CSD=Y_val[:]
    row_train_CSD,column_train_CSD=X_train_CSD.get_shape()
    
    t6=time.time()

    Y_predict_val_from_CNN_Model=model_main.predict(X_val_CSD, verbose=0)
    
    Y_predict_val_from_CNN_Model_Final=[0]*len(Y_predict_val_from_CNN_Model)
    for i in range(len(Y_predict_val_from_CNN_Model)):
        if Y_predict_val_from_CNN_Model[i]<0.5:
              Y_predict_val_from_CNN_Model[i]=0
        else:
              Y_predict_val_from_CNN_Model[i]=1
    for i in range(len(Y_predict_val_from_CNN_Model)):
       Y_predict_val_from_CNN_Model_Final[i]= int(Y_predict_val_from_CNN_Model[i])
        
    adjusted_rand_score_val=metrics.adjusted_rand_score(Y_val_CSD, Y_predict_val_from_CNN_Model_Final)
    adjusted_mutual_info_score_val=metrics.adjusted_mutual_info_score(Y_val_CSD, Y_predict_val_from_CNN_Model_Final) 
    homogeneity_score_val=metrics.homogeneity_score(Y_val_CSD, Y_predict_val_from_CNN_Model_Final) 
    fowlkes_mallows_score_val=metrics.fowlkes_mallows_score(Y_val_CSD, Y_predict_val_from_CNN_Model_Final) 
  
    for i in range(20):        #row_train
        Y_temp=Y_val_CSD.copy()

        row=X_train_CSD.getrow(i)
        X_temp = sp.lil.lil_matrix(sparse.csr_matrix(sparse.vstack((X_val_CSD, row))))
        Y_temp.append(Y_train_CSD[i])
        
        Y_predict_CNN_compute_CSD=model_main.predict(X_temp, verbose=0)
        
        Y_predict_temp=[0]*len(Y_predict_CNN_compute_CSD)
        
        for n in range(len(Y_predict_CNN_compute_CSD)):
            if Y_predict_CNN_compute_CSD[n]<0.5:
                  Y_predict_CNN_compute_CSD[n]=0
            else:
                   Y_predict_CNN_compute_CSD[n]=1
         
        for m in range(len(Y_predict_CNN_compute_CSD)):
            Y_predict_temp[m]= int(Y_predict_CNN_compute_CSD[m])

        adjusted_rand_score_temp=metrics.adjusted_rand_score(Y_temp, Y_predict_temp)
        adjusted_mutual_info_score_temp=metrics.adjusted_mutual_info_score(Y_temp, Y_predict_temp) 
        homogeneity_score_temp=metrics.homogeneity_score(Y_temp, Y_predict_temp) 
        fowlkes_mallows_score_temp=metrics.fowlkes_mallows_score(Y_temp, Y_predict_temp)
        
        landa1=abs(adjusted_rand_score_temp-adjusted_rand_score_val)
        landa2=abs(adjusted_mutual_info_score_temp-adjusted_mutual_info_score_val)
        landa3=abs(homogeneity_score_temp-homogeneity_score_val)
        landa4=abs(fowlkes_mallows_score_temp-fowlkes_mallows_score_val)
        
        sum_of_diffrences=landa1+landa2+landa3+landa4
        
        if sum_of_diffrences<0.1:
            X_val_CSD = sp.lil.lil_matrix(sparse.csr_matrix(sparse.vstack((X_val_CSD, row))))
            Y_val_CSD.append(Y_train_CSD[i])          
            Y_predict_CNN_inside_CSD=model_main.predict(X_val_CSD, verbose=0)
            
            Y_predict_CNN_inside_CSD_Final=[0]*len(Y_predict_CNN_inside_CSD)                   #Y_predict_CNN_inside
            for j in range(len(Y_predict_CNN_inside_CSD)):                                     #Y_predict_CNN_inside
                if Y_predict_CNN_inside_CSD[j]<0.5:
                      Y_predict_CNN_inside_CSD[j]=0
                else:
                      Y_predict_CNN_inside_CSD[j]=1
                       
            for k in range(len(Y_predict_CNN_inside_CSD)):                              #Y_predict_CNN_inside
                Y_predict_CNN_inside_CSD_Final[k]= int(Y_predict_CNN_inside_CSD[k])

            adjusted_rand_score_val=metrics.adjusted_rand_score(Y_val_CSD, Y_predict_CNN_inside_CSD_Final)
            adjusted_mutual_info_score_val=metrics.adjusted_mutual_info_score(Y_val_CSD, Y_predict_CNN_inside_CSD_Final) 
            homogeneity_score_val=metrics.homogeneity_score(Y_val_CSD, Y_predict_CNN_inside_CSD_Final) 
            fowlkes_mallows_score_val=metrics.fowlkes_mallows_score(Y_val_CSD, Y_predict_CNN_inside_CSD_Final) 
    t7=time.time()
    print("Time for Clustering Based Semi-supervised Defense =",t7-t6)
    print("                                                                   ")
    #**************************************************************************************** 
    X_train_Final_CSD= X_val_CSD.copy()  
    Y_train_Final_CSD=Y_val_CSD.copy()
    row_train_CSD_Final,col_train_CSD_Final=X_train_Final_CSD.get_shape()    
    
    model_main_CSD = Sequential()
    model_main_CSD.add(Embedding(row_train_CSD_Final, 8, input_length=col_train_CSD_Final))
    model_main_CSD.add(Conv1D(16,2, strides=2, padding='same'))
    model_main_CSD.add(MaxPooling1D(pool_size = (4), strides=(2)))
    model_main_CSD.add(Conv1D(32,2, strides=2, padding='same'))
    model_main_CSD.add(MaxPooling1D(pool_size = (4), strides=(2)))
    model_main_CSD.add(Conv1D(64,2, strides=2, padding='same'))
    model_main_CSD.add(Flatten())
    model_main_CSD.add(Dense(1, activation='sigmoid'))
    model_main_CSD.compile(optimizer='adam', loss='binary_crossentropy', metrics=['acc'])
    model_main_CSD.fit(X_train_Final_CSD, Y_train_Final_CSD, epochs=200, verbose=0)
       
    Y_test_predict_CSD=model_main_CSD.predict(X_test, verbose=0)
    
    Y_test_predict_CSD_Final=[0]*len(Y_test_predict_CSD)
    for i in range(len(Y_test_predict_CSD)):
        if Y_test_predict_CSD[i]<0.5:
              Y_test_predict_CSD[i]=0
        else:
              Y_test_predict_CSD[i]=1
    
    for i in range(len(Y_test_predict_CSD)):
        Y_test_predict_CSD_Final[i]= int(Y_test_predict_CSD[i])
        
    #*****************************************************Result of Model after Clustering Based Semi-supervised Defense(CSD)**************
    print("***********************Result of Model after Clustering Based Semi-supervised Defense(CSD)*************************************")  
    print("                                                                                                                               ")

    loss, accuracy = model_main_CSD.evaluate(X_train_Final_CSD, Y_train_Final_CSD, verbose=2)
    print('Accuracy for New Train set: %f' % (accuracy*100))
    print('Loss for New Train set: %f' % (loss))
    print("                                                                   ")

    loss, accuracy = model_main_CSD.evaluate(X_test, Y_test, verbose=2)
    print('Accuracy for Test set: %f' % (accuracy*100))
    print('Loss for Test set:: %f' % (loss))
    print("                                                                   ")

    TN_CSD, FP_CSD, FN_CSD, TP_CSD = confusion_matrix(Y_test,  Y_test_predict_CSD_Final).ravel()
    print("TN_CSD=",TN_CSD)
    print("FP_CSD=",FP_CSD)
    print("FN_CSD=",FN_CSD)
    print("TP_CSD=",TP_CSD)
    print("                                                                   ")

    if (FP_CSD+TN_CSD)>0:
         FPR_CSD=FP_CSD/(FP_CSD+TN_CSD)
         print("The FPR_CSD result=", FPR_CSD)
                
    if (FP_CSD+TN_CSD)>0:
         TPR_CSD=TP_CSD/(TP_CSD+FN_CSD)
         print("The TPR_CSD result=", TPR_CSD)
                
    if (TN_CSD+FP_CSD)>0:
        TNR_CSD=TN_CSD/(TN_CSD+FP_CSD)
        print("The TNR_CSD result=", TNR_CSD)
                
    if (FN_CSD+TP_CSD)>0:
        FNR_CSD=FN_CSD/(FN_CSD+TP_CSD)
        print("The FNR_CSD result=", FNR_CSD)
                
    if ((TN_CSD/(TN_CSD+FP_CSD))+(TP_CSD/(TP_CSD+FP_CSD)))>0:
        AUC_CSD=1/(2*((TN_CSD/(TN_CSD+FP_CSD))+(TP_CSD/(TP_CSD+FP_CSD))))
        print("The AUC_CSD result=", AUC_CSD)
                
    if  (TP_CSD+TN_CSD+FP_CSD+FN_CSD)>0:   
        ACC_CSD=(TP_CSD+TN_CSD)/(TP_CSD+TN_CSD+FP_CSD+FN_CSD)
        print("The ACC_CSD result=", ACC_CSD)
                
    if  ((TP_CSD+FP_CSD)*(TP_CSD+FN_CSD)*(TN_CSD+FP_CSD)*(TN_CSD+FN_CSD))>0:
        MCC_CSD=(TP_CSD*TN_CSD-FP_CSD*FN_CSD)/math.sqrt((TP_CSD+FP_CSD)*(TP_CSD+FN_CSD)*(TN_CSD+FP_CSD)*(TN_CSD+FN_CSD))
        print("The Matthews correlation coefficient result=", MCC_CSD)
    print("                                                                   ")
    print("************************************************End of Clustering Based Semi-supervised Defense(LSD)***************************")   
    print("                                                                                                                               ")
    print("                                                                                                                               ")
    print("                                                                                                                               ")
#******************************************************************************************************************************************
if __name__ == "__main__":
    main()
#******************************************************************************  