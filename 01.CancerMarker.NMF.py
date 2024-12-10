#!/usr/bin/env python
# coding: utf-8
# For Colon as example
# In[1]:


import os,copy,sys
from collections import Counter
import random
import multiprocessing
from glob import glob
import pickle

import pandas
import numpy

import scipy
from scipy import stats
from scipy.stats import rankdata
from scipy.stats import ttest_ind
from scipy.stats import pearsonr


from sklearn.linear_model import SGDRegressor
from sklearn.linear_model import PassiveAggressiveRegressor
from sklearn.linear_model import *
from sklearn.ensemble import  *
from sklearn.gaussian_process import *
from sklearn.neighbors import *
from sklearn.impute import KNNImputer
from sklearn.svm import *

from sklearn.model_selection import KFold,StratifiedKFold,train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV

from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import *

from sklearn.tree import *

from sklearn import metrics
from sklearn.metrics import confusion_matrix,roc_curve
from sklearn.metrics import roc_auc_score

from sklearn.metrics import silhouette_score    #### New

from sklearn.decomposition import NMF    #### New

from sklearn.neighbors import *

from sklearn.discriminant_analysis import *
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

from sklearn.kernel_ridge import KernelRidge

from sklearn.naive_bayes import *
from sklearn.naive_bayes import GaussianNB
from sklearn.naive_bayes import CategoricalNB
from sklearn.naive_bayes import ComplementNB



from statannot import add_stat_annotation

from statsmodels.stats.proportion import proportion_confint
from statsmodels.stats.contingency_tables import mcnemar

from plot_metric.functions import BinaryClassification
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
get_ipython().magic(u'matplotlib inline')

import warnings
warnings.filterwarnings("ignore")


import lightgbm
#import xgboost
#import catboost

import seaborn



 


# In[2]:


PARAMETERS = {
    'TOO':        ( 'Breast' , 'Colon' , 'Eso' , 'Gastric' , 'Lung' , 'Liver' , 'Pancreas' ),
    'Pathology':  ( 'Normal' , 'Cancer' , 'Benign' ),

    'Breast': {
        'Input_Matrix':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Breast____.10X_CpG_matrix.tsv',
        'Input_Sample':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Breast____.sample.tsv',
        'Input_CpG':     '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/CpG_Statistics_More/____Breast____.CpG.V2.tsv',
        'Output_Filter': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Breast____.10X_CpG_matrix.filtered.tsv',
        'Output_Impute': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Breast____.10X_CpG_matrix.imputed.tsv',
        'Output_Refer':  '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Breast____.10X_CpG_matrix.Imputed_with_InternalReference.tsv',
        'OutputPath':    '/Path/to/Analysis/31_01____NMF_with_InternalReference/Breast/',
        'K_range':       [ 2 , 3 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ] ,
        'Cutoff_Nonmissing':  0.7 ,
        'Cutoff_Sup':         0.65 , 
        'Cutoff_Inf':         0.35 ,
    } , 

    'Colon': {
        'Input_Matrix':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Colon____.10X_CpG_matrix.tsv',
        'Input_Sample':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Colon____.sample.tsv',
        'Input_CpG':     '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/CpG_Statistics_More/____Colon____.CpG.V2.tsv',
        'Output_Filter': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Colon____.10X_CpG_matrix.filtered.tsv',
        'Output_Impute': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Colon____.10X_CpG_matrix.imputed.tsv',
        'Output_Refer':  '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Colon____.10X_CpG_matrix.Imputed_with_InternalReference.tsv',
        'OutputPath':    '/Path/to/Analysis/31_01____NMF_with_InternalReference/Colon/',
        'K_range':       [ 2 , 3 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ] ,
        'Cutoff_Nonmissing':  0.7 ,
        'Cutoff_Sup':         0.65 , 
        'Cutoff_Inf':         0.35 ,
    } , 

    'Eso': {
        'Input_Matrix':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Eso____.10X_CpG_matrix.tsv',
        'Input_Sample':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Eso____.sample.tsv',
        'Input_CpG':     '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/CpG_Statistics_More/____Eso____.CpG.V2.tsv',
        'Output_Filter': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Eso____.10X_CpG_matrix.filtered.tsv',
        'Output_Impute': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Eso____.10X_CpG_matrix.imputed.tsv',
        'Output_Refer':  '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Eso____.10X_CpG_matrix.Imputed_with_InternalReference.tsv',
        'OutputPath':    '/Path/to/Analysis/31_01____NMF_with_InternalReference/Eso/',
        'K_range':       [ 2 , 3 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ] ,
        'Cutoff_Nonmissing':  0.7 ,
        'Cutoff_Sup':         0.65 , 
        'Cutoff_Inf':         0.35 ,
    } , 

    'Gastric': {
        'Input_Matrix':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Gastric____.10X_CpG_matrix.tsv',
        'Input_Sample':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Gastric____.sample.tsv',
        'Input_CpG':     '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/CpG_Statistics_More/____Gastric____.CpG.V2.tsv',
        'Output_Filter': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Gastric____.10X_CpG_matrix.filtered.tsv',
        'Output_Impute': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Gastric____.10X_CpG_matrix.imputed.tsv',
        'Output_Refer':  '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Gastric____.10X_CpG_matrix.Imputed_with_InternalReference.tsv',
        'OutputPath':    '/Path/to/Analysis/31_01____NMF_with_InternalReference/Gastric/',
        'K_range':       [ 2 , 3 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ] ,
        'Cutoff_Nonmissing':  0.7 ,
        'Cutoff_Sup':         0.65 , 
        'Cutoff_Inf':         0.35 ,
    } , 

    'Lung': {
        'Input_Matrix':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Lung____.10X_CpG_matrix.tsv',
        'Input_Sample':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Lung____.sample.tsv',
        'Input_CpG':     '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/CpG_Statistics_More/____Lung____.CpG.V2.tsv',
        'Output_Filter': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Lung____.10X_CpG_matrix.filtered.tsv',
        'Output_Impute': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Lung____.10X_CpG_matrix.imputed.tsv',
        'Output_Refer':  '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Lung____.10X_CpG_matrix.Imputed_with_InternalReference.tsv',
        'OutputPath':    '/Path/to/Analysis/31_01____NMF_with_InternalReference/Lung/',
        'K_range':       [ 2 , 3 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ] ,
        'Cutoff_Nonmissing':  0.7 ,
        'Cutoff_Sup':         0.65 , 
        'Cutoff_Inf':         0.35 ,
    } , 

    'Liver': {
        'Input_Matrix':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Liver____.10X_CpG_matrix.tsv',
        'Input_Sample':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Liver____.sample.tsv',
        'Input_CpG':     '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/CpG_Statistics_More/____Liver____.CpG.V2.tsv',
        'Output_Filter': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Liver____.10X_CpG_matrix.filtered.tsv',
        'Output_Impute': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Liver____.10X_CpG_matrix.imputed.tsv',
        'Output_Refer':  '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Liver____.10X_CpG_matrix.Imputed_with_InternalReference.tsv',
        'OutputPath':    '/Path/to/Analysis/31_01____NMF_with_InternalReference/Liver/',
        'K_range':       [ 2 , 3 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ] ,
        'Cutoff_Nonmissing':  0.7 ,
        'Cutoff_Sup':         0.65 , 
        'Cutoff_Inf':         0.35 ,
    } , 

    'Pancreas': {
        'Input_Matrix':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Pancreas____.10X_CpG_matrix.tsv',
        'Input_Sample':  '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Pancreas____.sample.tsv',
        'Input_CpG':     '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/CpG_Statistics_More/____Pancreas____.CpG.V2.tsv',
        'Output_Filter': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Pancreas____.10X_CpG_matrix.filtered.tsv',
        'Output_Impute': '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Pancreas____.10X_CpG_matrix.imputed.tsv',
        'Output_Refer':  '/Path/to/Analysis/31_01____NMF_with_InternalReference/____Pancreas____.10X_CpG_matrix.Imputed_with_InternalReference.tsv',
        'OutputPath':    '/Path/to/Analysis/31_01____NMF_with_InternalReference/Pancreas/',
        'K_range':       [ 2 , 3 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ] ,
        'Cutoff_Nonmissing':  0.7 ,
        'Cutoff_Sup':         0.65 , 
        'Cutoff_Inf':         0.35 ,
    },
}


# In[3]:


#### Input and output files
#### Parameters

#### ( 'Breast' , 'Colon' , 'Eso' , 'Gastric' , 'Lung' , 'Liver' , 'Pancreas' )

global parameters
#parameters = PARAMETERS['Breast']
parameters = PARAMETERS['Colon']
#parameters = PARAMETERS['Eso']
#parameters = PARAMETERS['Gastric']
#parameters = PARAMETERS['Lung']
#parameters = PARAMETERS['Liver']
#parameters = PARAMETERS['Pancreas']

for K in parameters.keys():
    print( K , '==>' , parameters[ K ] )
del K


# In[4]:


#### Global variable
global DATA
class DATA :
    pass


# In[5]:


def Read_Phenotype(  ):
    Phenotypes  =  pandas.read_csv(  parameters['Input_Sample']  ,  header=0  ,  index_col=0  ,  sep='\t'  )
    print(  Phenotypes.shape  )
    ##
    DATA.Phenotypes = Phenotypes
    ##
    DATA.Sample = list( Phenotypes.index.values )
    print(  len(  DATA.Sample  )  )
    ##
    ##


Read_Phenotype(  )


DATA.Phenotypes


# In[6]:


def Read_full_data():
    Available_CpG_list = [  ]
    ########
    count_line = 0         #### 计数
    count_available = 0    #### 计数
    print('Read from: ' , parameters['Input_CpG'] )
    inputTEXT = open( parameters['Input_CpG'] , 'r' )    #### 只读
    next( inputTEXT )    #### headline
    ## 0: CpG
    ## 1: Nonmissing/All
    ## 2: AUC
    while True:
        try:
            line2list = next( inputTEXT ).split('\n')[0].split('\r')[0].split('\t')
            CpG        =        str(line2list[0])
            Nonmissing = float( str(line2list[1]).split('=')[-1] )
            AUC        =      float(line2list[2])        if line2list[2]!='NA'    else 0.5
            count_line = count_line + 1    #### 计数
            ##
            if (
                ( Nonmissing >= parameters['Cutoff_Nonmissing'] )
                and
                ( ( AUC >= parameters['Cutoff_Sup'] )  or  ( AUC <= parameters['Cutoff_Inf'] ) )
            ):
                Available_CpG_list.append(  CpG  )
                count_available = count_available + 1    #### 计数
            else:
                pass
            ##
            if count_line % 100000 == 0:
                sys.stdout.write('\r        Line:'+str(count_line)+'        '  )
                sys.stdout.flush()
            else:
                pass
        except StopIteration:
            break
    inputTEXT.close()    #### 读毕
    print('\tChecking ==> Lines: ' , count_line  )
    print('\tChecking ==> Available: ' , count_available  )
    print('\tChecking ==> Available_CpG_list: ' , len(Available_CpG_list)  )
    print('\n')
    ##
    ##
    ##
    ##
    print('Read from: ' , parameters['Input_Matrix'] )
    Full_matrix_DF  =  pandas.read_csv(  parameters['Input_Matrix']  ,  header=0  ,  index_col=0  ,  sep='\t'  )
    print( Full_matrix_DF.shape )
    print('\n')
    ##
    ##
    ##
    ##
    print('Filtering ...... ')
    Available_matrix_DF  =  Full_matrix_DF.loc[ Available_CpG_list ]
    print( Available_matrix_DF.shape )
    print('\n')
    ##
    ##
    ##
    ##
    print('Write down: ' , parameters['Output_Filter'] )
    Available_matrix_DF.to_csv(  parameters['Output_Filter']  ,  sep='\t'  )
    ##
    ##
    ##
    ##



Read_full_data()
    


# In[7]:


def Read_bulked_data():
    print('Read from: ' , parameters['Output_Filter'] )
    X  =  pandas.read_csv( parameters['Output_Filter'] , header=0 , index_col=0  ,  sep='\t' ).T
    print(  X.shape  )
    ##
    ##
    X  =  X.loc[ DATA.Sample ]
    print(  X.shape  )
    ##
    ##
    ##X  =  X.dropna( axis=1 )
    print('Impute ......')
    Imputer = KNNImputer()
    X  =  pandas.DataFrame(  Imputer.fit_transform(  X  )  ,  columns=X.columns  ,  index=X.index  )
    del Imputer
    print(  X.shape  )
    ##
    print('Write down: ' , parameters['Output_Impute'] )
    X.T.to_csv(  parameters['Output_Impute']  ,  sep='\t'  )
    ##
    ##
    ##
    ##
    print('Add internal reference ......')
    Num_Sample  =  X.shape[0]
    print('Number of samples: ' , Num_Sample  )
    X['InternalReference:Meth=1']  =  [100]*Num_Sample
    X['InternalReference:Meth=0']  =  [ 0 ]*Num_Sample
    print(  X.shape  )
    ##
    print('Write down: ' , parameters['Output_Refer'] )
    X.T.to_csv(  parameters['Output_Refer']  ,  sep='\t'  )
    ##
    ##
    ##
    ##
    Marker  =  X.columns.values
    print(  len(Marker)  )
    DATA.X       =  X
    DATA.Marker  =  Marker
    



Read_bulked_data()


DATA.X



# In[8]:


def Calculate_silhouette( W , H ):
    K  =  W.shape[1]
    print( W.shape , H.shape , K )
    ##
    D = {  }
    for i in range( K ):
        D[i]  =  numpy.matrix(W[:,i]).T * numpy.matrix(H[i])
    ##
    D_stack = numpy.hstack( [ D[i]  for  i  in  range(K) ] ).T
    print( D_stack.shape )
    ##
    Label = [  ]
    for i in range( K ):
        for j in range( H.shape[1] ):
            Label.append( i )
    print( len(Label) , len(Label)==D_stack.shape[0] )
    ##
    silhouette  =  silhouette_score(  D_stack  ,  Label  )
    return  silhouette



def Calculate_Sparseness( X , W , H ):
    rss = numpy.linalg.norm( X-numpy.matrix(W)*numpy.matrix(H) )
    Sparseness  =  1 - rss/numpy.linalg.norm(X)
    return  Sparseness



# In[9]:


def Plot_by_page(  Hmatrix_DF  ,  Compo_Num  ,  OutputFigure_PDF  ):
    print('Plot figure: ' , OutputFigure_PDF )
    PDF  =  PdfPages( OutputFigure_PDF )
    ##
    for j in range( Compo_Num ):
        Component = 'Component:'+str(j)
        fig , ax  =  plt.subplots( 1 , 1 , figsize=(10,10) )
        ax.hist(  Hmatrix_DF.loc[ Component ]   ,   bins=50  )
        ax.set_title(  'Component '+str(j)+' of '+str(Compo_Num)  )
        PDF.savefig( fig )
        del Component
        ##
    ##
    PDF.close()


# In[10]:


def Decompose():
    ####DATA.silhouette = {  }
    DATA.rss = {  }
    ####DATA.Sparseness = {  }
    DATA.StatFrac = {  }
    ########
    for K in parameters['K_range']:
        print('K ==> ' , K )
        ##
        model = NMF( n_components=K , init='random' )
        W = model.fit_transform( DATA.X )
        H = model.components_
        R = DATA.X - numpy.matrix(W)*numpy.matrix(H)
        ##
        print('score')
        DATA.rss[ K ]         =  model.reconstruction_err_
        ####DATA.silhouette[ K ]  =  Calculate_silhouette( W , H )
        ####DATA.Sparseness[ K ]  =  Calculate_Sparseness( DATA.X , W , H )
        ##
        print('Write down')
        W_DF  =  pandas.DataFrame(  W  ,  columns=['Component:'+str(i) for i in range(K)]  ,  index=DATA.Sample                              )
        H_DF  =  pandas.DataFrame(  H  ,  columns=DATA.Marker                              ,  index=['Component:'+str(i) for i in range(K)]  )
        OutputFile_W    =  parameters['OutputPath']+str(K)+'____W____Sample_Components.tsv'
        OutputFile_W_T  =  parameters['OutputPath']+str(K)+'____W____Components_Sample.tsv'
        OutputFile_H    =  parameters['OutputPath']+str(K)+'____H____Components_Marker.tsv'
        OutputFile_H_T  =  parameters['OutputPath']+str(K)+'____H____Marker_Components.tsv'
        OutputFile_R    =  parameters['OutputPath']+str(K)+'____Residual____Sample_Marker.tsv'
        OutputFile_R_T  =  parameters['OutputPath']+str(K)+'____Residual____Marker_Sample.tsv'
        W_DF.to_csv(    OutputFile_W    ,  sep='\t'  )
        W_DF.T.to_csv(  OutputFile_W_T  ,  sep='\t'  )
        H_DF.to_csv(    OutputFile_H    ,  sep='\t'  )
        H_DF.T.to_csv(  OutputFile_H_T  ,  sep='\t'  )
        R.to_csv(       OutputFile_R    ,  sep='\t'  )
        R.T.to_csv(     OutputFile_R_T  ,  sep='\t'  )
        ##
        ##
        ##
        Output_IR   =  parameters['OutputPath']+str(K)+'____InRef____Components_Marker.tsv'
        Output_IR_T =  parameters['OutputPath']+str(K)+'____InRef____Marker_Components.tsv'
        H_DF[ ['InternalReference:Meth=1','InternalReference:Meth=0'] ].to_csv(    Output_IR    ,  sep='\t'  )
        H_DF[ ['InternalReference:Meth=1','InternalReference:Meth=0'] ].T.to_csv(  Output_IR_T  ,  sep='\t'  )
        ##
        ##
        ##
        print('Scaling ...... ')
        for i in range(K):
            Component = 'Component:'+str(i)
            ##multiplier  =  100 / float(H_DF.loc[ Component ]['InternalReference:Meth=1'])    
            multiplier  =  100 / float(H_DF.loc[ Component ]['InternalReference:Meth=1'])        if H_DF.loc[ Component ]['InternalReference:Meth=1'] > 0        else 100/max(H_DF.loc[Component])
            ####multiplier  =  100 / max( max(H_DF.loc[Component]) , float(H_DF.loc[ Component ]['InternalReference:Meth=1']) )
            H_DF.loc[ Component ]  =  H_DF.loc[ Component ] * multiplier
            W_DF[ Component ]      =  W_DF[ Component ]     / multiplier
            del Component , multiplier
        ##
        rss_check         =      numpy.linalg.norm( DATA.X  -  numpy.matrix(W_DF) * numpy.matrix(H_DF) )
        Sparseness_check  =  1 - numpy.linalg.norm( DATA.X  -  numpy.matrix(W_DF) * numpy.matrix(H_DF) ) / numpy.linalg.norm( DATA.X )
        Sparseness_check2 =  1 - numpy.linalg.norm( R )                                                  / numpy.linalg.norm( DATA.X )
        print( model.reconstruction_err_  ,  rss_check , abs(model.reconstruction_err_ - rss_check)<0.000001  , Sparseness_check2 , Sparseness_check ,  abs(Sparseness_check2 - Sparseness_check)<0.00000001 )
        del  rss_check , Sparseness_check , Sparseness_check2
        ##
        Output_Scaled_W    =  parameters['OutputPath']+str(K)+'____scaled_W____Sample_Components.tsv'
        Output_Scaled_W_T  =  parameters['OutputPath']+str(K)+'____scaled_W____Components_Sample.tsv'
        Output_Scaled_H    =  parameters['OutputPath']+str(K)+'____scaled_H____Components_Marker.tsv'
        Output_Scaled_H_T  =  parameters['OutputPath']+str(K)+'____scaled_H____Marker_Components.tsv'
        W_DF.to_csv(    Output_Scaled_W    ,  sep='\t'  )
        W_DF.T.to_csv(  Output_Scaled_W_T  ,  sep='\t'  )
        H_DF.to_csv(    Output_Scaled_H    ,  sep='\t'  )
        H_DF.T.to_csv(  Output_Scaled_H_T  ,  sep='\t'  )
        ##
        ##
        ##
        Output_PDF_H  =  parameters['OutputPath']+str(K)+'____scaled_H.pdf'
        Plot_by_page(  H_DF  ,  K  ,  Output_PDF_H  )
        ##
        ##
        ##
        DATA.StatFrac[ K ]  =  list( W_DF.sum( axis=1 ) )
        ##
        ##
        ##
        del model
        del W , H , R
        del W_DF , H_DF
        del OutputFile_W , OutputFile_W_T , OutputFile_H , OutputFile_H_T , OutputFile_R , OutputFile_R_T
        del Output_IR , Output_IR_T
        del Output_Scaled_W , Output_Scaled_W_T , Output_Scaled_H , Output_Scaled_H_T
        del Output_PDF_H
        ##
    ##






Decompose()





# In[11]:


def Output_Score():
    OutputFile_score  =  parameters['OutputPath']+'____Score_MoreComponents____.tsv'
    print('Write down: ' , OutputFile_score  )
    outputTEXT = open( OutputFile_score , 'w' )
    ####outputTEXT.write(  ('\t').join([ 'Num_Components' , 'rss' , 'Sparseness' , 'silhouette' ]) +'\n'  )
    outputTEXT.write(  ('\t').join([ 'Num_Components' , 'rss' ]) +'\n'  )
    for K in parameters['K_range']:
        ####outputTEXT.write(  ('\t').join([  str(K)  ,  str(DATA.rss[K])  ,  str(DATA.Sparseness[K])  ,  str(DATA.silhouette[K])  ])  +'\n'  )
        outputTEXT.write(  ('\t').join([  str(K)  ,  str(DATA.rss[K])  ])  +'\n'  )
    outputTEXT.close()
    



Output_Score()




# In[12]:


def Write_down_fraction_stastics():
    OutputFile_Frac  =  parameters['OutputPath']+'____Statistics_FractionSUM____.tsv'
    print('Write down: ' , OutputFile_Frac  )
    outputTEXT = open( OutputFile_Frac , 'w' )
    outputTEXT.write(  ('\t').join([ 'Num_Complement' , 'Mean' , 'STD' , 'Minimum' , '25_percentile' , 'Median' , '75_percentile' , 'Maximum' ])  +'\n'  )
    for K in parameters['K_range']:
        Mean     =  numpy.mean( DATA.StatFrac[ K ] )
        STD      =  numpy.std(  DATA.StatFrac[ K ] )
        ##
        Minimum  =  min( DATA.StatFrac[ K ] )
        Q25      =  numpy.percentile( DATA.StatFrac[ K ]  ,  25 )
        Median   =  numpy.median(     DATA.StatFrac[ K ]        )
        Q75      =  numpy.percentile( DATA.StatFrac[ K ]  ,  75 )
        Maximum  =  max( DATA.StatFrac[ K ] )
        ##
        outputTEXT.write(  ('\t').join([ str(K) , str(Mean) , str(STD) , str(Minimum) , str(Q25) , str(Median) , str(Q75) , str(Maximum) ])  +'\n'  )
        ##
    outputTEXT.close()



Write_down_fraction_stastics()



# In[ ]:





# In[ ]:





# In[ ]:




