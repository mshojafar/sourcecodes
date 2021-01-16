# -*- coding: utf-8 -*-
 
"""
Created on Fri May 25 12:03:10 2018

@author: Rahim
#this approch use the distribution of Benign data to poison thet test data
"""
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
from scipy import sparse
import pandas as pd
import numpy as np
import random
import time
import argparse
import math
from numpy import *
import os.path as osp
import scipy.sparse as sp
import pickle
from sklearn import metrics
from sklearn.metrics import accuracy_score
#******************************************************************************
CLASS = 'class'
CLASS_BEN = 'B'
CLASS_MAL = 'M'
DATA = 'data'
#********************************************Functions that will be used in this program*******************************************************************************************
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-tables', nargs='*', dest='input_tables')

    args = parser.parse_args()

    return args
#******************************************************************************
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
def build_table(tables):
    full_table = dict()

    file_set = set()
    
    for table in tables:
        file_set.update(table.keys())
        for key, val in table.items():
            full_table[key] = val
              
    files = list(file_set)
    return full_table, files
#******************************************************************************
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
#*****************************************************************Main Function*********************************************************************************************************
def main():
    args = parse_args()

    tables = []
    f_set = set()
    
    #read the data
    for t_files in args.input_tables:
        table, features = read_table(t_files)
        f_set = f_set.union(features)
        tables.append(table)
    print("                                                                                         ")
    print("                                                                                         ")
    print("*****************************************************************************************")
    print("********Using Benign Distribution + Random Forest Classifier + GAN countermeasure********") 
    print("*****************************************************************************************")

    #*build table from data and convert to matrix 
    full_table, files = build_table(tables)
    files.sort()
    features = list(f_set)
    features.sort()
    mat, cl = convert_to_matrix(full_table, features, files) 
   
    #Doing feature Ranking on all of the Data
    print("************************Doing feature Ranking on all of the Data*************************")
    t0=time.time()
    r_features,ranked_index = relevant_features(mat, cl, features)
    t1=time.time()
    print("Time of Feature Ranking=",t1-t0)
    print("******************************************************************************************")
 
    original_selected=ranked_index[1:301]
    data = sparse.lil_matrix(sparse.csr_matrix(mat)[:,original_selected])
    seed = 10
    test_size = 0.2
    X_train, X_test, Y_train, Y_test= train_test_split(data, cl, test_size= test_size, random_state=seed)
    test_size = 0.25
    X_train, X_val, Y_train, Y_val= train_test_split(X_train, Y_train, test_size= test_size, random_state=seed)  
    #**************************************************************************
    num_trees = 100
    max_features = 3
    t0=time.time()
    kfold = KFold(n_splits=10, random_state=10)
    model = RandomForestClassifier(n_estimators=num_trees, max_features=max_features)
    model.fit(X_train, Y_train)
    t1=time.time()
    print("Time for Clssification Algorithm is runing on 300 high-ranked features =",t1-t0)
    print("************************************Result without attack *******************************************************************************************")
    # compute Classification Accuracy in train and test and Validation
    scoring = 'accuracy'
    results = model_selection.cross_val_score(model, X_train,Y_train, cv=kfold, scoring=scoring)
    print(("The accuracy of Classification in train: %.3f (%.3f)") % (results.mean(), results.std()))
    #********************* compute Classification Accuracy in validation***************************
    scoring = 'accuracy'
    results = model_selection.cross_val_score(model, X_val,Y_val, cv=kfold, scoring=scoring)
    print(("The accuracy of Classification in validation: %.3f (%.3f)") % (results.mean(), results.std()))
    #********************* compute Classification Accuracy in test*********************************
    scoring = 'accuracy'
    results = model_selection.cross_val_score(model, X_test,Y_test, cv=kfold, scoring=scoring)
    print(("The accuracy of Classification in test: %.3f (%.3f)") % (results.mean(), results.std()))
    #********************* compute Classification Accuracy in Validation***************************
    predictions_val = model.predict(X_val)
    print("classification_report by validation:")
    print(classification_report(Y_val, predictions_val))
    #********************* compute Classification Accuracy in train********************************
    predictions = model.predict(X_test)
    print("classification_report by test:")
    print(classification_report(Y_test, predictions))
    #********************* compute Logarithmic Loss in Train*********************************
    scoring = 'neg_log_loss'
    results = model_selection.cross_val_score(model, X_train,Y_train, cv=kfold, scoring=scoring)
    print(("The Loss of Classification in train data: %.3f (%.3f)") % (results.mean(), results.std()))
    #********************* compute Logarithmic Loss in validation**************************** 
    scoring = 'neg_log_loss'
    results = model_selection.cross_val_score(model, X_val,Y_val, cv=kfold, scoring=scoring)
    print(("The Loss of Classification in validation data:: %.3f (%.3f)") % (results.mean(), results.std()))
    #********************* compute Logarithmic Loss in Test***********************************
    scoring = 'neg_log_loss'
    results = model_selection.cross_val_score(model, X_test,Y_test, cv=kfold, scoring=scoring)
    print(("The Loss of Classification in test data:: %.3f (%.3f)") % (results.mean(), results.std()))
    #********************* compute Area Under ROC Curve in Train******************************
    scoring = 'roc_auc'
    results = model_selection.cross_val_score(model, X_train,Y_train, cv=kfold, scoring=scoring)
    print(("The Area Under ROC Curve in Train: %.3f (%.3f)") % (results.mean(), results.std()))
    #********************* compute Area Under ROC Curve in Validation*************************
    scoring = 'roc_auc'
    results = model_selection.cross_val_score(model, X_val,Y_val, cv=kfold, scoring=scoring)
    print(("The Area Under ROC Curve in Validation: %.3f (%.3f)") % (results.mean(), results.std()))
    #********************* compute Area Under ROC Curve in Test*******************************
    scoring = 'roc_auc'
    results = model_selection.cross_val_score(model, X_test,Y_test, cv=kfold, scoring=scoring)
    print(("The Area Under ROC Curve in test: %.3f (%.3f)") % (results.mean(), results.std()))
    #*****************************Compute FPR and TPR in Validation**************************
    cm=confusion_matrix(Y_test, predictions)
    print("confusion_matrix=")
    print(cm)
    TP=cm[0][0]
    print("TP=",TP)
    FP=cm[0][1]
    print("FP=",FP)
    FN=cm[1][0]
    print("FN=",FN)
    TN=cm[1][1]
    print("TN=",TN)
    FPR=FP/(FP+TN)
    print("The FPR result=", FPR)
    TPR=TP/(TP+FN)
    print("The TPR result=", TPR)
    
    TNR=TN/(TN+FP)
    print("The TNR result=", TNR)
    
    FNR=FN/(FN+TP)
    print("The FNR result=", FNR)
    
    AUC=1/(2*((TN/(TN+FP))+(TP/(TP+FP))))
    print("The AUC result=", AUC)
    
    ACC=(TP+TN)/(TP+TN+FP+FN)
    print("The ACC result=", ACC)
    
    MCC=(TP*TN-FP*FN)/math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    print("The Matthews correlation coefficient result=", MCC)
    
    print("*******************************End of Result without attack:*****************************************************************************************")
    #*************************************************************************************************************************************************************
    # finding Malware of test data
    malware_test= sparse.lil_matrix(X_test)
    cl_malware=list()
    z_m=0    
    count_m=0
    for i, j in enumerate(Y_test):
        if j == 1:
            delete_row_lil(malware_test, i-count_m)
            count_m=count_m+1
        else:
            cl_malware.insert(z_m, 1)
            z_m=z_m+1

    #**************************
    #finding Benign of test data
    benign_test = sparse.lil_matrix(X_test)
    cl_benign=list()
    z_b=0    
    count_b=0
    for i, j in enumerate(Y_test):
        if j == 0:
            delete_row_lil(benign_test.tolil(), i-count_b)
            count_b=count_b+1
        else:
            cl_benign.insert(z_b, 1)
            z_b=z_b+1
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
    cl_benign_train=list()
    for i, j in enumerate(cl_X_train):
        if j == 0:
            delete_row_lil(benign_train, i-count_b)
            count_b=count_b+1
        else:
            cl_benign_train.insert(z_b, 1)
            z_b=z_b+1
    print("***********Size of Each Data Part:**********")        
    print("malware_train=", malware_train.get_shape())
    print("benign_train=", benign_train.get_shape())
    print("malware_test=", malware_test.get_shape())
    print("benign_test=", benign_test.get_shape())
    #***************************************************
    t0=time.time()
    ranked_features_in_benign,ranked_index_of_benign = relevant_features(benign_train,cl_benign_train, features)  
    t1=time.time()
    print("Time for Ranking benign_train to find important features =",t1-t0)
    #***************************************************************************************************************************************************************
    numbers=list()
    numbers=[3,6,9,12,15,18,21,24,27,30,60]
    X_test = sp.lil.lil_matrix(X_test)
    
    for loop in range(10):
        print("************************************************************************************************************************************************************************************")
        print("Result related to loop number : ",loop)
       
        Malware_Test=sparse.lil_matrix(malware_test.copy()) 
        row_of_Malware,column_of_Malware=Malware_Test.get_shape()
        index_of_row=list(range(row_of_Malware))
        random.shuffle(index_of_row)

        number_of_row_to_change=int(row_of_Malware/10)
        selected_row=index_of_row[0:number_of_row_to_change]
        
        for i, v in enumerate(numbers):
            print("*****************************************************************************************************************************************************")
            print("*********************selected features :",int(v) )
            print("************************************Result after attack *************************")
            max_index_of_column=int(v)+1
            t0=time.time()
            rw_test,cl_test=X_test.get_shape()
            poison_data=sp.lil.lil_matrix((0,cl_test),dtype=np.int8)
            Malware_Test=sparse.lil_matrix(malware_test.copy())

            counter_of_poisoned_point=0
            
            for m,value in enumerate(selected_row):
                flag=0
                for i, j in enumerate(ranked_index_of_benign[1:max_index_of_column]):
                    for k,l in enumerate(original_selected):
                      if j==l:     
                          if Malware_Test[value,l]==0:
                                Malware_Test[value,l]=1
                                flag=1
                if flag==1:
                    counter_of_poisoned_point=counter_of_poisoned_point+1
           
            
            Benign_Test=sparse.lil_matrix(benign_test.copy()) 
            poison_data = sp.lil.lil_matrix(sparse.csr_matrix(sparse.vstack((Benign_Test, Malware_Test))))
            r,w=poison_data.get_shape()
            Y_test=Y_test[0:r]
            
            t1=time.time()                   
            print("Time related to applying attack in this number of Features= ",t1-t0)
            
            print("Number of poisoned Malware= ",counter_of_poisoned_point)
            #********************* compute Classification Accuracy in test*********************************
            scoring = 'accuracy'
            results = model_selection.cross_val_score(model, poison_data,Y_test, cv=kfold, scoring=scoring)
            print(("The accuracy of Classification in test: %.3f (%.3f)") % (results.mean(), results.std()))
           
            #********************* compute Classification Accuracy in train********************************
            predictions = model.predict(poison_data)
            print("classification_report by test:")
            print(classification_report(Y_test, predictions))
          
            #********************* compute Logarithmic Loss in Test***********************************
            scoring = 'neg_log_loss'
            results = model_selection.cross_val_score(model, poison_data,Y_test, cv=kfold, scoring=scoring)
            print(("The Loss of Classification in test data:: %.3f (%.3f)") % (results.mean(), results.std()))
           
            #********************* compute Area Under ROC Curve in Test*******************************
            scoring = 'roc_auc'
            results = model_selection.cross_val_score(model, poison_data,Y_test, cv=kfold, scoring=scoring)
            print(("The Area Under ROC Curve in test: %.3f (%.3f)") % (results.mean(), results.std()))
            #*****************************Compute FPR and TPR in Validation**************************
            cm=confusion_matrix(Y_test, predictions)
            print("confusion_matrix=")
            print(cm)
            TP=cm[0][0]
            print("TP=",TP)
            FP=cm[0][1]
            print("FP=",FP)
            FN=cm[1][0]
            print("FN=",FN)
            TN=cm[1][1]
            print("TN=",TN)
            FPR=FP/(FP+TN)
            print("The FPR result=", FPR)
            
            TPR=TP/(TP+FN)
            print("The TPR result=", TPR)
            
            TNR=TN/(TN+FP)
            print("The TNR result=", TNR)
            
            FNR=FN/(FN+TP)
            print("The FNR result=", FNR)
            
            AUC=1/(2*((TN/(TN+FP))+(TP/(TP+FP))))
            print("The AUC result=", AUC)
            
            ACC=(TP+TN)/(TP+TN+FP+FN)
            print("The ACC result=", ACC)
            
            MCC=(TP*TN-FP*FN)/math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
            print("The Matthews correlation coefficient result=", MCC)
            
            print("********************************************************************************")
            print("*******************Result after applying GAN countermeasure**************")
            t0=time.time()


            model2 = ExtraTreesClassifier(n_estimators=250,random_state=0)
            model2.fit(benign_train, cl_benign_train)
            importances = model2.feature_importances_
            indices = np.argsort(importances)[::-1]
            

            importance_of_Features_in_benign_train=list()
            for f in range(60):
                importance_of_Features_in_benign_train.append(indices[f])

            #******************************Runing the Logistic Regression and  finding Some Sampels Near to Hyperplain*****************************
            poison_model = LogisticRegression() 
            poison_model.fit(X_train,Y_train)
            print("Result related to Logistic Regression:")
            scoring = 'accuracy'
            poison_results = model_selection.cross_val_score(poison_model, X_train,Y_train, cv=kfold, scoring=scoring)
            print(("The accuracy of Classification in train: %.3f (%.3f)") % (poison_results.mean(), poison_results.std()))
            #********************* compute Logistic Regression Accuracy in validation without change ***************************
            scoring = 'accuracy'
            results = model_selection.cross_val_score(poison_model, X_val,Y_val, cv=kfold, scoring=scoring)
            print(("The accuracy of Classification in validation: %.3f (%.3f)") % (poison_results.mean(), poison_results.std()))
            #********************* compute Logistic Regression Accuracy in test without change *********************************
            scoring = 'accuracy'
            results = model_selection.cross_val_score(poison_model, X_test,Y_test, cv=kfold, scoring=scoring)
            print(("The accuracy of Classification in test: %.3f (%.3f)") % (poison_results.mean(), poison_results.std()))
            #**********************Declration of Variables for finding desision value *************
            print("**************************************************************************************************")           
            temp=sparse.lil_matrix(X_train)
            a,b=temp.get_shape()
            decision_value=np.array([])
            selected_cl_malware_train=list()
            selected_malware_train = sparse.lil_matrix(X_train)       
            #**********************finding malware_train and related desision value **********************************
            counter_of_malware_train=0
            count_deleted=0
            for j in range(a):
                row=temp.getrow(j)
                if Y_train[j]==0:
                    decision_value=np.append(decision_value,poison_model.decision_function(row))
                    selected_cl_malware_train.insert(counter_of_malware_train, 0)
                    counter_of_malware_train=counter_of_malware_train+1
                else:
                    delete_row_lil(selected_malware_train.tolil(), j-count_deleted)
                    count_deleted=count_deleted+1 
            #**********************sort the absolute value of decision value for malware_train*************************      
            decision_value=np.absolute(decision_value)
            indices=decision_value.argsort()
              
            #************** Declration of Variables for selecting data*************************************************
            number_of_row_malware_train,number_of_column_malware_train=malware_train.get_shape()
               
            number_of_row_selected_malware_train=int(number_of_row_malware_train/10)
            
            #****************Selecting index related to 10 percent of malware_train with minimum decision value*******
            Selected_rows_as_less_likely=list()
            Selected_rows_as_less_likely=indices[:number_of_row_selected_malware_train]
                                 
            Malware_less_likely=sp.lil.lil_matrix((0, number_of_column_malware_train), dtype=np.int8)
            cl_less_likely=list()
            counter_for_cl_less_likely=0            
            for i,row_number in enumerate(Selected_rows_as_less_likely):
                selected_row=malware_train.getrow(row_number)
                Malware_less_likely= sp.lil.lil_matrix(sparse.csr_matrix(sparse.vstack((selected_row, Malware_less_likely))))
                cl_less_likely.insert(counter_for_cl_less_likely,0)
                counter_for_cl_less_likely=counter_for_cl_less_likely+1
                
            number_of_row_in_Malware_less_likely,number_of_column_in_Malware_less_likely=Malware_less_likely.get_shape()
            #****************finding Benign like samples********************************************************************************
            poisoned_data=sp.lil.lil_matrix((0, number_of_column_malware_train), dtype=np.int8)
            c=0
            for counter_of_Malware_less_likely in range(number_of_row_in_Malware_less_likely):
                selected_sample=Malware_less_likely.getrow(counter_of_Malware_less_likely)

                
                c=0
                for S in range(number_of_column_in_Malware_less_likely):     
                        index_for_change=random.randint(0,number_of_column_in_Malware_less_likely-1)
                        if selected_sample[0,index_for_change]==0:
                            selected_sample[0,index_for_change]=1 
                        label=model.predict(selected_sample)
                        if label==int(1):
                            poisoned_data= sp.lil.lil_matrix(sparse.csr_matrix(sparse.vstack((selected_sample, poisoned_data))))
                            c=c+1
                            break       
             
            Number_of_row_in_poisoned_data,Number_of_column_in_poison_data=poisoned_data.get_shape()
            Y_poisoin=list()
            for index in range(Number_of_row_in_poisoned_data):
                Y_poisoin.append(0)
            #***************************************************************************************************************************
            
            poisoned_data_X=poisoned_data.copy()
            poisoned_data_Y=Y_poisoin[:]
            second_test_set=0.2
            X_poisoned_train, X_poisoned_test, Y_poisoned_train, Y_poisoned_test= train_test_split(poisoned_data_X, poisoned_data_Y, test_size= second_test_set, random_state=seed)
            poison_data_for_retraining = sp.lil.lil_matrix(sparse.csr_matrix(sparse.vstack((X_train, X_poisoned_train))))
            poison_Class_for_retraining = Y_train + Y_poisoned_train
            
            num_trees = 100
            max_features = 3
            kfold = KFold(n_splits=10, random_state=10)
            model_for_counter_measure = RandomForestClassifier(n_estimators=num_trees, max_features=max_features)
            model_for_counter_measure.fit(poison_data_for_retraining,poison_Class_for_retraining)
    
    
            poison_data_for_test_after_retraining = sp.lil.lil_matrix(sparse.csr_matrix(sparse.vstack((X_test, X_poisoned_test))))
            poison_Class_for_test_after_retraining= Y_test + Y_poisoned_test
    
            t1=time.time()
            print("Time related to applying  GAN countermeasure in this number of Features: ",t1-t0)
            #********************* compute Classification Accuracy in test*********************************
            scoring = 'accuracy'
            results = model_selection.cross_val_score(model_for_counter_measure, poison_data_for_test_after_retraining,poison_Class_for_test_after_retraining, cv=kfold, scoring=scoring)
            print(("The accuracy of Classification in test: %.3f (%.3f)") % (results.mean(), results.std()))
           
            #********************* compute Classification Accuracy in train********************************
            predictions = model.predict(poison_data_for_test_after_retraining)
            print("classification_report by test:")
            print(classification_report(poison_Class_for_test_after_retraining, predictions))
          
            #********************* compute Logarithmic Loss in Test***********************************
            scoring = 'neg_log_loss'
            results = model_selection.cross_val_score(model_for_counter_measure, poison_data_for_test_after_retraining , poison_Class_for_test_after_retraining, cv=kfold, scoring=scoring)
            print(("The Loss of Classification in test data:: %.3f (%.3f)") % (results.mean(), results.std()))
           
            #********************* compute Area Under ROC Curve in Test*******************************
            scoring = 'roc_auc'
            results = model_selection.cross_val_score(model_for_counter_measure, poison_data_for_test_after_retraining , poison_Class_for_test_after_retraining, cv=kfold, scoring=scoring)
            print(("The Area Under ROC Curve in test: %.3f (%.3f)") % (results.mean(), results.std()))
            #*****************************Compute FPR and TPR in Validation**************************
            cm=confusion_matrix(poison_Class_for_test_after_retraining, predictions)
            print("confusion_matrix=")
            print(cm)
            TP=cm[0][0]
            print("TP=",TP)
            FP=cm[0][1]
            print("FP=",FP)
            FN=cm[1][0]
            print("FN=",FN)
            TN=cm[1][1]
            print("TN=",TN)
            FPR=FP/(FP+TN)
            print("The FPR result=", FPR)
            
            TPR=TP/(TP+FN)
            print("The TPR result=", TPR)
            
            TNR=TN/(TN+FP)
            print("The TNR result=", TNR)
            
            FNR=FN/(FN+TP)
            print("The FNR result=", FNR)
            
            AUC=1/(2*((TN/(TN+FP))+(TP/(TP+FP))))
            print("The AUC result=", AUC)
            
            ACC=(TP+TN)/(TP+TN+FP+FN)
            print("The ACC result=", ACC)
            
            MCC=(TP*TN-FP*FN)/math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
            print("The Matthews correlation coefficient result=", MCC)

            print("Result Related to this Numbers of features is finished:",int(v))
            Malware_Test=sparse.lil_matrix(malware_test.copy()) 
            selected_row=index_of_row[0:number_of_row_to_change]
            original_selected=ranked_index[1:301]
        print("End of loop number : ",loop)
        print("************************************************************************************************************************************************************************************")
#********************************************************************************************************************************************************************
if __name__ == "__main__":
    main()
#******************************************************************************
