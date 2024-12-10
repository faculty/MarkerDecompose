#!/usr/bin/env python
# coding: utf-8

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





from pandarallel import pandarallel
pandarallel.initialize(nb_workers=50)



 


# In[2]:


global parameters

parameters = {
    'TOO':        ( 'Breast' , 'Colon' , 'Eso' , 'Gastric' , 'Lung' , 'Liver' , 'Pancreas' ),
    ####'Pathology':  ( 'Normal' , 'Cancer' , 'Benign' ),
    'Pathology':  ( 'Normal' , 'Cancer' ),
    ####'Type':   ['Breast:Normal', 'Breast:Cancer', 'Colon:Normal', 'Colon:Cancer', 'Eso:Normal', 'Eso:Cancer', 'Gastric:Normal', 'Gastric:Cancer', 'Lung:Normal', 'Lung:Cancer', 'Lung:Benign', 'Liver:Normal', 'Liver:Cancer', 'Pancreas:Normal', 'Pancreas:Cancer', 'Plasma:Normal', 'Blood:Normal'],
    'Type':       ['Breast:Normal', 'Breast:Cancer', 'Colon:Normal', 'Colon:Cancer', 'Eso:Normal', 'Eso:Cancer', 'Gastric:Normal', 'Gastric:Cancer', 'Lung:Normal', 'Lung:Cancer', 'Liver:Normal', 'Liver:Cancer', 'Pancreas:Normal', 'Pancreas:Cancer', 'Plasma:Normal', 'Blood:Normal'],

    'Input_Matrix':  '/Path/to/Analysis/51_01____Bulked_RRBS_matrix____7cancers/7Cancers_with_HealthyBlood.10X_CpG_matrix.tsv',
    'Input_Stats':   '/Path/to/Analysis/51_01____Bulked_RRBS_matrix____7cancers/7Cancers_with_HealthyBlood.10X_CpG_statistics.tsv',

    'Output_Filter': '/Path/to/Analysis/51_02____NMF_with_InternalReference/7Cancers_with_HealthyBlood.10X_CpG_matrix.filtered.tsv',
    'Output_Impute': '/Path/to/Analysis/51_02____NMF_with_InternalReference/7Cancers_with_HealthyBlood.10X_CpG_matrix.imputed.tsv',
    'Output_Refer':  '/Path/to/Analysis/51_02____NMF_with_InternalReference/7Cancers_with_HealthyBlood.10X_CpG_matrix.Imputed_with_InternalReference.tsv',
    'OutputPath':    '/Path/to/Analysis/51_02____NMF_with_InternalReference/NMF/',
    'K_range':       [ 5 , 6 , 7 , 8 , 9 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ] ,

    'Cutoff_Nonmissing':  0.5  ,
    'Cutoff_Num_TOO':     7    , 
    'Cutoff_and_or':      'and',
    'Cutoff_Pvalue':      0.05 ,
}


# In[3]:


for K in parameters.keys():
    print( K , '==>' , parameters[ K ] )
del K


# In[4]:


#### Global variable
global DATA
class DATA :
    pass


# In[ ]:





# In[5]:


def Scan_line(  x  ):
    Nonmiss_N_Breast , Nonmiss_C_Breast = float(x['Nonmissing|Breast:Normal'  ].split('=')[-1])  ,  float(x['Nonmissing|Breast:Cancer'  ].split('=')[-1])
    Nonmiss_N_Colon  , Nonmiss_C_Colon  = float(x['Nonmissing|Colon:Normal'   ].split('=')[-1])  ,  float(x['Nonmissing|Colon:Cancer'   ].split('=')[-1])
    Nonmiss_N_Eso    , Nonmiss_C_Eso    = float(x['Nonmissing|Eso:Normal'     ].split('=')[-1])  ,  float(x['Nonmissing|Eso:Cancer'     ].split('=')[-1])
    Nonmiss_N_Gas    , Nonmiss_C_Gas    = float(x['Nonmissing|Gastric:Normal' ].split('=')[-1])  ,  float(x['Nonmissing|Gastric:Cancer' ].split('=')[-1])
    Nonmiss_N_Lung   , Nonmiss_C_Lung   = float(x['Nonmissing|Lung:Normal'    ].split('=')[-1])  ,  float(x['Nonmissing|Lung:Cancer'    ].split('=')[-1])
    Nonmiss_N_Liver  , Nonmiss_C_Liver  = float(x['Nonmissing|Liver:Normal'   ].split('=')[-1])  ,  float(x['Nonmissing|Liver:Cancer'   ].split('=')[-1])
    Nonmiss_N_Panc   , Nonmiss_C_Panc   = float(x['Nonmissing|Pancreas:Normal'].split('=')[-1])  ,  float(x['Nonmissing|Pancreas:Cancer'].split('=')[-1])
    ##
    Normal_Pvalue  =  float(x['ANOVA_P|Normal'])    if str(x['ANOVA_P|Normal']) != 'NA'    else 1.0
    Cancer_Pvalue  =  float(x['ANOVA_P|Cancer'])    if str(x['ANOVA_P|Cancer']) != 'NA'    else 1.0
    ##
    ##
    Count_Nonmiss_N = 0    #### 计数
    Count_Nonmiss_C = 0    #### 计数
    Switch_Nonmiss  = 'off'
    Switch_Pvalue   = 'off'
    FinalSwitch     = 0
    ##
    for Nonmiss_value in  ( Nonmiss_N_Breast , Nonmiss_N_Colon , Nonmiss_N_Eso , Nonmiss_N_Gas , Nonmiss_N_Lung , Nonmiss_N_Liver , Nonmiss_N_Panc ):
        if  Nonmiss_value  >= parameters['Cutoff_Nonmissing']:
            Count_Nonmiss_N = Count_Nonmiss_N + 1    #### 计数
        else:
            pass
    for Nonmiss_value in  ( Nonmiss_C_Breast , Nonmiss_C_Colon , Nonmiss_C_Eso , Nonmiss_C_Gas , Nonmiss_C_Lung , Nonmiss_C_Liver , Nonmiss_C_Panc ):
        if  Nonmiss_value  >= parameters['Cutoff_Nonmissing']:
            Count_Nonmiss_C = Count_Nonmiss_C + 1    #### 计数
        else:
            pass
    ##
    if  parameters['Cutoff_and_or'] == 'and':
        if  ( Count_Nonmiss_N >= parameters['Cutoff_Num_TOO'] )    and   ( Count_Nonmiss_C >= parameters['Cutoff_Num_TOO'] ):
            Switch_Nonmiss  = 'on'
        else:
            pass
    elif parameters['Cutoff_and_or'] == 'or':
        if  ( Count_Nonmiss_N >= parameters['Cutoff_Num_TOO'] )    or    ( Count_Nonmiss_C >= parameters['Cutoff_Num_TOO'] ):
            Switch_Nonmiss  = 'on'
        else:
            pass
    else:
        print('Error ==> and/or ' , parameters['Cutoff_Num_TOO']  )
    ##
    if  ( Normal_Pvalue < parameters['Cutoff_Pvalue'] )    or    ( Cancer_Pvalue < parameters['Cutoff_Pvalue'] ):
        Switch_Pvalue  = 'on'
    else:
        pass
    ##
    ##
    ##
    if  ( Switch_Nonmiss == 'on' )  and  ( Switch_Pvalue == 'on' ):
        FinalSwitch = 1
    else:
        FinalSwitch = 0
    ##
    return  FinalSwitch
 


# In[6]:


def Scan_Statistics(  ):
    Matrix_DF  =  pandas.read_csv(  parameters['Input_Stats']  ,  header=0  ,  index_col=0  ,  sep='\t'  )
    print(  Matrix_DF.shape  )
    ##
    Output_DF  =  Matrix_DF.parallel_apply(  Scan_line  ,  axis=1  )
    Available_CpG  =  Output_DF[ Output_DF==1 ].index.values
    print('Filtered CpG: ' , len(Available_CpG)  )
    ##
    ##
    DATA.Available_CpG  =  list(Available_CpG)


Scan_Statistics(  )




# In[7]:


def Read_full_data():
    print('Read from: ' , parameters['Input_Matrix'] )
    Full_matrix_DF  =  pandas.read_csv(  parameters['Input_Matrix']  ,  header=0  ,  index_col=0  ,  sep='\t'  )
    print( Full_matrix_DF.shape )
    print('\n')
    ##
    ##
    ##
    ##
    print('Filtering ...... ')
    Available_matrix_DF  =  Full_matrix_DF.loc[ DATA.Available_CpG ]
    print( Available_matrix_DF.shape )
    print('\n')
    ##
    ##
    ##
    print('Sorting column head ...... ')
    Type_SampleID____dict = {  }
    for Type in parameters['Type']:
        Type_SampleID____dict[ Type ]  =  [  ]
    ##
    for SampleID in Available_matrix_DF.columns.values:
        Type = SampleID.split(' ---- ')[0]
        if Type in Type_SampleID____dict.keys():
            Type_SampleID____dict[ Type ].append(  SampleID  )
    ##
    SampleID_sortedList = [  ]
    for Type in parameters['Type']:
        for SampleID in Type_SampleID____dict[ Type ]:
            SampleID_sortedList.append( SampleID )
    print( len(SampleID_sortedList) )
    ##
    DATA.SampleID_sorted  =  tuple(SampleID_sortedList)
    ##
    Available_matrix_DF  =  Available_matrix_DF[  SampleID_sortedList  ]
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


DATA.SampleID_sorted
    


# In[9]:



def Read_bulked_data():
    print('Read from: ' , parameters['Output_Filter'] )
    X  =  pandas.read_csv( parameters['Output_Filter'] , header=0 , index_col=0  ,  sep='\t' ).T
    print(  X.shape  )
    ##
    ##
    X  =  X.loc[ list(DATA.SampleID_sorted) ]
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



'''
def Read_bulked_data():
    print('Read again from: ' , parameters['Output_Refer'] )
    X  =  pandas.read_csv( parameters['Output_Refer'] , header=0 , index_col=0  ,  sep='\t' ).T
    print(  X.shape  )
    ##
    ##
    Marker  =  X.columns.values
    print(  len( Marker )  )
    DATA.X       =  X
    DATA.Marker  =  Marker
'''



Read_bulked_data()


DATA.X



# In[10]:


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



# In[11]:


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


# In[12]:


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
        W_DF  =  pandas.DataFrame(  W  ,  columns=['Component:'+str(i) for i in range(K)]  ,  index=list(DATA.SampleID_sorted)               )
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





# In[13]:


def Output_Score():
    OutputFile_score  =  parameters['OutputPath']+'____Score____.tsv'
    print('Write down: ' , OutputFile_score  )
    outputTEXT = open( OutputFile_score , 'w' )
    ####outputTEXT.write(  ('\t').join([ 'Num_Components' , 'rss' , 'Sparseness' , 'silhouette' ]) +'\n'  )
    outputTEXT.write(  ('\t').join([ 'Num_Components' , 'rss' ]) +'\n'  )
    for K in parameters['K_range']:
        ####outputTEXT.write(  ('\t').join([  str(K)  ,  str(DATA.rss[K])  ,  str(DATA.Sparseness[K])  ,  str(DATA.silhouette[K])  ])  +'\n'  )
        outputTEXT.write(  ('\t').join([  str(K)  ,  str(DATA.rss[K])  ])  +'\n'  )
    outputTEXT.close()
    



Output_Score()




# In[14]:


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




