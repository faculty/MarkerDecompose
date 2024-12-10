#!/usr/bin/python3
# -*- coding: utf-8 -*-



print('Import: sys , os , math')
import sys
from optparse import OptionParser
import os
import glob
import numpy
import pandas
import math

print('Import: scipy')
import scipy
from scipy import stats
from scipy.stats import rankdata
from scipy.stats import ttest_ind
from scipy.stats import pearsonr

print('Import: matplotlib')
import matplotlib
import matplotlib.pyplot
from matplotlib.backends.backend_pdf import PdfPages

print('Import: others')
import seaborn

import warnings
warnings.filterwarnings("ignore")









'''
parameters_Colon____colorectal_1X  =  {
    'Healthy_blood': ( 'NormalPlasma.1X' , 'HeaBloodCell.1X' ),
    '2nd_NMF':       ( 'decom' , 'merge' ),
    'High/Low':      ( 'HighMeth' , 'LowMeth' ),
    'TOO':             'Colon',
    ##
    'CFEA_data_type':  'CFEA_RRBS',
    'CFEA_X':          '1X',
    'CFEA_Stat':     ( 'AUC' , 'Distance' , 'logFC' , 'Utest' ),
    'Input_CFEA':   '/Path/to/Analysis/33_02____Bulked_RRBS_CpG_HighMeth_LowMeth/Add_CFEA/',
    'Input_Compo':  '/Path/to/Analysis/34_01____Component_VS_HealthyBlood____NMF_used_CpG/',
    'Input_Path':   '/Path/to/Analysis/34_02____Component_VS_HealthyBlood____NMF_used_CpG____Add_CFEA/',
    'Output_Path':  '/Path/to/Analysis/34_03____Distribution_cfDNA_CN/',
    ##
    'Top': 1000,
}
parameters_Colon____colorectal_10X  =  {
    'Healthy_blood': ( 'NormalPlasma.1X' , 'HeaBloodCell.1X' ),
    '2nd_NMF':       ( 'decom' , 'merge' ),
    'High/Low':      ( 'HighMeth' , 'LowMeth' ),
    'TOO':             'Colon',
    ##
    'CFEA_data_type':  'CFEA_RRBS',
    'CFEA_X':          '10X',
    'CFEA_Stat':     ( 'AUC' , 'Distance' , 'logFC' , 'Utest' ),
    'Input_CFEA':   '/Path/to/Analysis/33_02____Bulked_RRBS_CpG_HighMeth_LowMeth/Add_CFEA/',
    'Input_Compo':  '/Path/to/Analysis/34_01____Component_VS_HealthyBlood____NMF_used_CpG/',
    'Input_Path':   '/Path/to/Analysis/34_02____Component_VS_HealthyBlood____NMF_used_CpG____Add_CFEA/',
    'Output_Path':  '/Path/to/Analysis/34_03____Distribution_cfDNA_CN/',
    ##
    'Top': 1000,
}
'''

parameters_Colon____blood_10X____colorectal_1X  =  {
    'Healthy_blood': ( 'NormalPlasma.10X' , 'HeaBloodCell.10X' ),
    '2nd_NMF':       ( 'decom' , 'merge' ),
    'High/Low':      ( 'HighMeth' , 'LowMeth' ),
    'TOO':             'Colon',
    ##
    'CFEA_data_type':  'CFEA_RRBS',
    'CFEA_X':          '1X',
    'CFEA_Stat':     ( 'AUC' , 'Distance' , 'logFC' , 'Utest' ),
    'Input_CFEA':   '/Path/to/Analysis/33_02____Bulked_RRBS_CpG_HighMeth_LowMeth/Add_CFEA/',
    'Input_Compo':  '/Path/to/Analysis/34_01____Component_VS_HealthyBlood____NMF_used_CpG/',
    'Input_Path':   '/Path/to/Analysis/34_02____Component_VS_HealthyBlood____NMF_used_CpG____Add_CFEA/',
    'Output_Path':  '/Path/to/Analysis/34_03____Distribution_cfDNA_CN/',
    ##
    'Top': 1000,
    'Major_Stat':   'Ratio_mean:distance' , #### Ratio_mean:distance    Ratio_mean:logFC    Ratio_median:distance    Ratio_median:logFC    Ratio_across_samples:distance    Ratio_across_samples:logFC ####
}
parameters_Colon____blood_10X____colorectal_10X  =  {
    'Healthy_blood': ( 'NormalPlasma.10X' , 'HeaBloodCell.10X' ),
    '2nd_NMF':       ( 'decom' , 'merge' ),
    'High/Low':      ( 'HighMeth' , 'LowMeth' ),
    'TOO':             'Colon',
    ##
    'CFEA_data_type':  'CFEA_RRBS',
    'CFEA_X':          '10X',
    'CFEA_Stat':     ( 'AUC' , 'Distance' , 'logFC' , 'Utest' ),
    'Input_CFEA':   '/Path/to/Analysis/33_02____Bulked_RRBS_CpG_HighMeth_LowMeth/Add_CFEA/',
    'Input_Compo':  '/Path/to/Analysis/34_01____Component_VS_HealthyBlood____NMF_used_CpG/',
    'Input_Path':   '/Path/to/Analysis/34_02____Component_VS_HealthyBlood____NMF_used_CpG____Add_CFEA/',
    'Output_Path':  '/Path/to/Analysis/34_03____Distribution_cfDNA_CN/',
    ##
    'Top': 1000,
    'Major_Stat':   'Ratio_mean:distance' , #### Ratio_mean:distance    Ratio_mean:logFC    Ratio_median:distance    Ratio_median:logFC    Ratio_across_samples:distance    Ratio_across_samples:logFC ####
}








'''
parameters_Lung____lung_1X  =  {
    'Healthy_blood': ( 'NormalPlasma.1X' , 'HeaBloodCell.1X' ),
    '2nd_NMF':       ( 'decom' , 'merge' ),
    'High/Low':      ( 'HighMeth' , 'LowMeth' ),
    'TOO':             'Lung',
    ##
    'CFEA_data_type':  'CFEA_RRBS',
    'CFEA_X':          '1X',
    'CFEA_Stat':     ( 'AUC' , 'Distance' , 'logFC' , 'Utest' ),
    'Input_CFEA':   '/Path/to/Analysis/33_02____Bulked_RRBS_CpG_HighMeth_LowMeth/Add_CFEA/',
    'Input_Compo':  '/Path/to/Analysis/34_01____Component_VS_HealthyBlood____NMF_used_CpG/',
    'Input_Path':   '/Path/to/Analysis/34_02____Component_VS_HealthyBlood____NMF_used_CpG____Add_CFEA/',
    'Output_Path':  '/Path/to/Analysis/34_03____Distribution_cfDNA_CN/',
    ##
    'Top': 1000,
}
parameters_Lung____lung_10X  =  {
    'Healthy_blood': ( 'NormalPlasma.1X' , 'HeaBloodCell.1X' ),
    '2nd_NMF':       ( 'decom' , 'merge' ),
    'High/Low':      ( 'HighMeth' , 'LowMeth' ),
    'TOO':             'Lung',
    ##
    'CFEA_data_type':  'CFEA_RRBS',
    'CFEA_X':          '10X',
    'CFEA_Stat':     ( 'AUC' , 'Distance' , 'logFC' , 'Utest' ),
    'Input_CFEA':   '/Path/to/Analysis/33_02____Bulked_RRBS_CpG_HighMeth_LowMeth/Add_CFEA/',
    'Input_Compo':  '/Path/to/Analysis/34_01____Component_VS_HealthyBlood____NMF_used_CpG/',
    'Input_Path':   '/Path/to/Analysis/34_02____Component_VS_HealthyBlood____NMF_used_CpG____Add_CFEA/',
    'Output_Path':  '/Path/to/Analysis/34_03____Distribution_cfDNA_CN/',
    ##
    'Top': 1000,
}
'''

parameters_Lung____blood_10X____lung_1X  =  {
    'Healthy_blood': ( 'NormalPlasma.10X' , 'HeaBloodCell.10X' ),
    '2nd_NMF':       ( 'decom' , 'merge' ),
    'High/Low':      ( 'HighMeth' , 'LowMeth' ),
    'TOO':             'Lung',
    ##
    'CFEA_data_type':  'CFEA_RRBS',
    'CFEA_X':          '1X',
    'CFEA_Stat':     ( 'AUC' , 'Distance' , 'logFC' , 'Utest' ),
    'Input_CFEA':   '/Path/to/Analysis/33_02____Bulked_RRBS_CpG_HighMeth_LowMeth/Add_CFEA/',
    'Input_Compo':  '/Path/to/Analysis/34_01____Component_VS_HealthyBlood____NMF_used_CpG/',
    'Input_Path':   '/Path/to/Analysis/34_02____Component_VS_HealthyBlood____NMF_used_CpG____Add_CFEA/',
    'Output_Path':  '/Path/to/Analysis/34_03____Distribution_cfDNA_CN/',
    ##
    'Top': 1000,
    'Major_Stat':   'Ratio_mean:distance' , #### Ratio_mean:distance    Ratio_mean:logFC    Ratio_median:distance    Ratio_median:logFC    Ratio_across_samples:distance    Ratio_across_samples:logFC ####
}
parameters_Lung____blood_10X____lung_10X  =  {
    'Healthy_blood': ( 'NormalPlasma.10X' , 'HeaBloodCell.10X' ),
    '2nd_NMF':       ( 'decom' , 'merge' ),
    'High/Low':      ( 'HighMeth' , 'LowMeth' ),
    'TOO':             'Lung',
    ##
    'CFEA_data_type':  'CFEA_RRBS',
    'CFEA_X':          '10X',
    'CFEA_Stat':     ( 'AUC' , 'Distance' , 'logFC' , 'Utest' ),
    'Input_CFEA':   '/Path/to/Analysis/33_02____Bulked_RRBS_CpG_HighMeth_LowMeth/Add_CFEA/',
    'Input_Compo':  '/Path/to/Analysis/34_01____Component_VS_HealthyBlood____NMF_used_CpG/',
    'Input_Path':   '/Path/to/Analysis/34_02____Component_VS_HealthyBlood____NMF_used_CpG____Add_CFEA/',
    'Output_Path':  '/Path/to/Analysis/34_03____Distribution_cfDNA_CN/',
    ##
    'Top': 1000,
    'Major_Stat':   'Ratio_mean:distance' , #### Ratio_mean:distance    Ratio_mean:logFC    Ratio_median:distance    Ratio_median:logFC    Ratio_across_samples:distance    Ratio_across_samples:logFC ####
}





def main_loop( parameters ):
    HL_Stat____CFEA_DF = {  }
    ########
    for HL in parameters['High/Low']:
        HL_Stat____CFEA_DF[ HL ] = {  }
        for Stat in parameters['CFEA_Stat']:
            InputFile = parameters['Input_CFEA'] + str(parameters['TOO'])+'.'+str(HL)+'.Sorted_by_'+str(Stat)+'.Add_'+str(parameters['CFEA_data_type'])+'.'+str(parameters['CFEA_X'])+'.tsv'
            print('Read CFEA from: ' , InputFile )
            CFEA_DF  =  pandas.read_csv(  InputFile  ,  header=0  ,  index_col=0  ,  sep='\t'  ).iloc[  : parameters['Top']]
            print( CFEA_DF.shape )
            HL_Stat____CFEA_DF[ HL ][ Stat ]  =  CFEA_DF.rename( index={X:'Bulked_Data_sorted_by_CFEA_'+str(Stat)+'|'+X    for    X    in    CFEA_DF.index} )
            del InputFile , CFEA_DF
        ##
    ##
    ##
    ##
    ##
    for HeaBlood in parameters['Healthy_blood']:
        for NMF_2nd in parameters['2nd_NMF']:
            K____list = [ ]
            K_compo____dict = {  }
            ########
            Input_Components  =  parameters['Input_Compo'] + str(HeaBlood)+'____VS____Truncated_Meth.2nd_by_'+str(NMF_2nd)+'.'+str(parameters['TOO'])+'/____Component____.tsv'
            print('Read components from: ' , Input_Components )
            inputTEXT = open( Input_Components , 'r' )    #### 只读
            next( inputTEXT )    #### headline
            while True:
                try:
                    line2list = next(inputTEXT).split('\n')[0].split('\r')[0].split('\t')
                    K      =    int(line2list[0])
                    compo  =  tuple(line2list[1].split(','))
                    ##
                    K____list.append(  int(K)  )
                    K_compo____dict[ K ]  =  tuple( compo )
                    ##
                    del K , compo
                except StopIteration:
                    break
            inputTEXT.close()    #### 读毕
            del Input_Components
            ##
            for K in K____list:
                print( K , '  ==>  ' , K_compo____dict[K] )
            ##
            ##
            ##
            InputPath  = parameters['Input_Path']  + str(HeaBlood)+'____2nd_by_'+str(NMF_2nd)+'.'+str(parameters['TOO'])+'____'+str(parameters['CFEA_data_type'])+'.'+str(parameters['CFEA_X'])+'/'
            OutputPath = parameters['Output_Path'] + str(HeaBlood)+'____2nd_by_'+str(NMF_2nd)+'.'+str(parameters['TOO'])+'____'+str(parameters['CFEA_data_type'])+'.'+str(parameters['CFEA_X'])+'/'
            print('Input path:  ' , InputPath  )
            print('Output path: ' , OutputPath )
            ##
            if os.path.exists( OutputPath ) == True:
                pass
            else:
                os.mkdir(  OutputPath  )
                print( 'Make Output path: ' , OutputPath )
            ##
            ##
            ##
            for K in K____list:
                print( 'NMF: ' , K )
                ##
                CDT_DF  =  pandas.DataFrame(  index=['HighMeth:'+str(i+1) for i in range(parameters['Top'])]+['']*50+['LowMeth:'+str(i+1) for i in range(parameters['Top'])]  )
                ##
                Input_Sub_Path  = str(InputPath) +'NMF_'+str(K)+'/'
                Output_PDF      = str(OutputPath)+'NMF_'+str(K)+'.pdf'
                ##print('Output PDF: ' , Output_PDF )
                ##PDF  =  PdfPages( Output_PDF )
                ##
                Output_Stat     = str(OutputPath)+'NMF_'+str(K)+'.statistics.tsv'
                outputTEXT_stat = open( Output_Stat , 'w' )
                outputTEXT_stat.write(  ('\t').join([ 'Component' , 'HighMeth/LowMeth' , 'Component-Blood' , 'Bulked_data_sorted_by' , 'Bulked_data' , 't-test' , 'U-test' ])  +'\n'  )
                ##
                ##
                ##
                Output_MultiCom = str(OutputPath)+'NMF_'+str(K)+'.MultiComparison.tsv'
                outputTEXT_MCom = open( Output_MultiCom , 'w' )
                outputTEXT_MCom.write(  ('\t').join([ 
                    'Component' , 'HighMeth/LowMeth' , 'Component/Bulked' , 
                    'CFEA_RRBS|Ratio_mean:distance'           , 'CFEA_RRBS|Ratio_mean:logFC'           ,
                    'CFEA_RRBS|Ratio_median:distance'         , 'CFEA_RRBS|Ratio_median:logFC'         ,
                    'CFEA_RRBS|Ratio_across_samples:distance' , 'CFEA_RRBS|Ratio_across_samples:logFC' ,
                ])  +'\n'  )
                ##
                ##
                ##
                for Component in K_compo____dict[ K ]:
                    CDT_DF[ str(Component)+' ---- Compo-Blood'     ]  =  numpy.nan
                    CDT_DF[ str(Component)+' ---- Bulked_AUC'      ]  =  numpy.nan
                    CDT_DF[ str(Component)+' ---- Bulked_Distance' ]  =  numpy.nan
                    CDT_DF[ str(Component)+' ---- Bulked_logF'     ]  =  numpy.nan
                    CDT_DF[ str(Component)+' ---- Bulked_Utest'    ]  =  numpy.nan
                    CDT_DF[ str(Component)+' ---- Gap'             ]  =  numpy.nan
                    ##
                    for High_Low in parameters['High/Low']:
                        InputFile   =  str(Input_Sub_Path)  +str(Component).replace(':','_') +'.'+str(High_Low)+'.tsv'
                        print('Read from:  ' , InputFile  )
                        ##
                        Compo_DF  =  pandas.read_csv( InputFile  ,  header=0  ,  index_col=0  ,  sep='\t'  ).iloc[  : parameters['Top']]
                        ##
                        Compo_VS_CFEA   =                                      Compo_DF[ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|'+str(parameters['Major_Stat']) ]
                        Bulked_AUC      =  HL_Stat____CFEA_DF[ High_Low ][ 'AUC'      ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|'+str(parameters['Major_Stat']) ]
                        Bulked_Distance =  HL_Stat____CFEA_DF[ High_Low ][ 'Distance' ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|'+str(parameters['Major_Stat']) ]
                        Bulked_logFC    =  HL_Stat____CFEA_DF[ High_Low ][ 'logFC'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|'+str(parameters['Major_Stat']) ]
                        Bulked_Utest    =  HL_Stat____CFEA_DF[ High_Low ][ 'Utest'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|'+str(parameters['Major_Stat']) ]
                        ##
                        ##
                        ##
                        Compo_VS_CFEA____list    =  list( Compo_VS_CFEA.dropna()   )
                        Bulked_AUC____list       =  list( Bulked_AUC.dropna()      )
                        Bulked_Distance____list  =  list( Bulked_Distance.dropna() )
                        Bulked_logFC____list     =  list( Bulked_logFC.dropna()    )
                        Bulked_Utest____list     =  list( Bulked_Utest.dropna()    )
                        Compo_Mean  , Compo_STD   =  numpy.mean( Compo_VS_CFEA____list )   ,  numpy.std( Compo_VS_CFEA____list )
                        B_AUC_Mean  , B_AUC_STD   =  numpy.mean( Bulked_AUC____list )      ,  numpy.std( Bulked_AUC____list )
                        B_Dist_Mean , B_Dist_STD  =  numpy.mean( Bulked_Distance____list ) ,  numpy.std( Bulked_Distance____list )
                        B_FC_Mean   , B_FC_STD    =  numpy.mean( Bulked_logFC____list )    ,  numpy.std( Bulked_logFC____list )
                        B_U_Mean    , B_U_STD     =  numpy.mean( Bulked_Utest____list )    ,  numpy.std( Bulked_Utest____list )
                        Compo_VS_B_AUC____Ttest   =  str( scipy.stats.ttest_ind(    Compo_VS_CFEA____list , Bulked_AUC____list , equal_var=False )[1] )
                        Compo_VS_B_AUC____Utest   =  str( scipy.stats.mannwhitneyu( Compo_VS_CFEA____list , Bulked_AUC____list )[1] )
                        Compo_VS_B_Dist____Ttest  =  str( scipy.stats.ttest_ind(    Compo_VS_CFEA____list , Bulked_Distance____list , equal_var=False )[1] )
                        Compo_VS_B_Dist____Utest  =  str( scipy.stats.mannwhitneyu( Compo_VS_CFEA____list , Bulked_Distance____list )[1] )
                        Compo_VS_B_logFC____Ttest =  str( scipy.stats.ttest_ind(    Compo_VS_CFEA____list , Bulked_logFC____list , equal_var=False )[1] )
                        Compo_VS_B_logFC____Utest =  str( scipy.stats.mannwhitneyu( Compo_VS_CFEA____list , Bulked_logFC____list )[1] )
                        Compo_VS_B_Utest____Ttest =  str( scipy.stats.ttest_ind(    Compo_VS_CFEA____list , Bulked_Utest____list , equal_var=False )[1] )
                        Compo_VS_B_Utest____Utest =  str( scipy.stats.mannwhitneyu( Compo_VS_CFEA____list , Bulked_Utest____list )[1] )
                        ##
                        outputTEXT_stat.write(  ('\t').join([ str(Component) , str(High_Low) , str(Compo_Mean)+'+/-'+str(Compo_STD) , 'AUC'      ,  str(B_AUC_Mean)+'+/-'+str(B_AUC_STD)  , str(Compo_VS_B_AUC____Ttest)   , str(Compo_VS_B_AUC____Utest)   ])  +'\n'  )
                        outputTEXT_stat.write(  ('\t').join([ str(Component) , str(High_Low) , str(Compo_Mean)+'+/-'+str(Compo_STD) , 'Distance' , str(B_Dist_Mean)+'+/-'+str(B_Dist_STD) , str(Compo_VS_B_Dist____Ttest)  , str(Compo_VS_B_Dist____Utest)  ])  +'\n'  )
                        outputTEXT_stat.write(  ('\t').join([ str(Component) , str(High_Low) , str(Compo_Mean)+'+/-'+str(Compo_STD) , 'logFC'    ,   str(B_FC_Mean)+'+/-'+str(B_FC_STD)   , str(Compo_VS_B_logFC____Ttest) , str(Compo_VS_B_logFC____Utest) ])  +'\n'  )
                        outputTEXT_stat.write(  ('\t').join([ str(Component) , str(High_Low) , str(Compo_Mean)+'+/-'+str(Compo_STD) , 'U-test'   ,    str(B_U_Mean)+'+/-'+str(B_U_STD)    , str(Compo_VS_B_Utest____Ttest) , str(Compo_VS_B_Utest____Utest) ])  +'\n'  )
                        ##
                        ##
                        ##
                        Compo_VS_CFEA____dict = {
                            'mean:distance':    Compo_DF[ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_mean:distance' ].mean()           ,
                            'mean:logFC':       Compo_DF[ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_mean:logFC'    ].mean()           ,
                            'median:distance':  Compo_DF[ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_median:distance' ].mean()         ,
                            'median:logFC':     Compo_DF[ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_median:logFC'    ].mean()         ,
                            'across:distance':  Compo_DF[ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_across_samples:distance' ].mean() ,
                            'across:logFC':     Compo_DF[ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_across_samples:logFC'    ].mean() ,
                        }
                        Bulked_AUC____dict = {
                            'mean:distance':    HL_Stat____CFEA_DF[ High_Low ][ 'AUC'      ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_mean:distance' ].mean()           ,
                            'mean:logFC':       HL_Stat____CFEA_DF[ High_Low ][ 'AUC'      ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_mean:logFC'    ].mean()           ,
                            'median:distance':  HL_Stat____CFEA_DF[ High_Low ][ 'AUC'      ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_median:distance' ].mean()         ,
                            'median:logFC':     HL_Stat____CFEA_DF[ High_Low ][ 'AUC'      ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_median:logFC'    ].mean()         ,
                            'across:distance':  HL_Stat____CFEA_DF[ High_Low ][ 'AUC'      ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_across_samples:distance' ].mean() ,
                            'across:logFC':     HL_Stat____CFEA_DF[ High_Low ][ 'AUC'      ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_across_samples:logFC'    ].mean() ,
                        }
                        Bulked_Distance____dict = {
                            'mean:distance':    HL_Stat____CFEA_DF[ High_Low ][ 'Distance' ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_mean:distance' ].mean()           ,
                            'mean:logFC':       HL_Stat____CFEA_DF[ High_Low ][ 'Distance' ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_mean:logFC'    ].mean()           ,
                            'median:distance':  HL_Stat____CFEA_DF[ High_Low ][ 'Distance' ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_median:distance' ].mean()         ,
                            'median:logFC':     HL_Stat____CFEA_DF[ High_Low ][ 'Distance' ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_median:logFC'    ].mean()         ,
                            'across:distance':  HL_Stat____CFEA_DF[ High_Low ][ 'Distance' ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_across_samples:distance' ].mean() ,
                            'across:logFC':     HL_Stat____CFEA_DF[ High_Low ][ 'Distance' ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_across_samples:logFC'    ].mean() ,
                        }
                        Bulked_logFC____dict = {
                            'mean:distance':    HL_Stat____CFEA_DF[ High_Low ][ 'logFC'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_mean:distance' ].mean()           ,
                            'mean:logFC':       HL_Stat____CFEA_DF[ High_Low ][ 'logFC'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_mean:logFC'    ].mean()           ,
                            'median:distance':  HL_Stat____CFEA_DF[ High_Low ][ 'logFC'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_median:distance' ].mean()         ,
                            'median:logFC':     HL_Stat____CFEA_DF[ High_Low ][ 'logFC'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_median:logFC'    ].mean()         ,
                            'across:distance':  HL_Stat____CFEA_DF[ High_Low ][ 'logFC'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_across_samples:distance' ].mean() ,
                            'across:logFC':     HL_Stat____CFEA_DF[ High_Low ][ 'logFC'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_across_samples:logFC'    ].mean() ,
                        }
                        Bulked_Utest____dict = {
                            'mean:distance':    HL_Stat____CFEA_DF[ High_Low ][ 'Utest'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_mean:distance' ].mean()           ,
                            'mean:logFC':       HL_Stat____CFEA_DF[ High_Low ][ 'Utest'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_mean:logFC'    ].mean()           ,
                            'median:distance':  HL_Stat____CFEA_DF[ High_Low ][ 'Utest'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_median:distance' ].mean()         ,
                            'median:logFC':     HL_Stat____CFEA_DF[ High_Low ][ 'Utest'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_median:logFC'    ].mean()         ,
                            'across:distance':  HL_Stat____CFEA_DF[ High_Low ][ 'Utest'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_across_samples:distance' ].mean() ,
                            'across:logFC':     HL_Stat____CFEA_DF[ High_Low ][ 'Utest'    ][ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_across_samples:logFC'    ].mean() ,
                        }
                        ##
                        ##outputTEXT_MCom.write(('\t').join([ 'Component' , 'HighMeth/LowMeth' , 'Component/Bulked' , 'CFEA_RRBS|Ratio_mean:distance' , 'CFEA_RRBS|Ratio_mean:logFC' , 'CFEA_RRBS|Ratio_median:distance' , 'CFEA_RRBS|Ratio_median:logFC' , 'CFEA_RRBS|Ratio_across_samples:distance' , 'CFEA_RRBS|Ratio_across_samples:logFC' , ])  +'\n'  )
                        outputTEXT_MCom.write(  ('\t').join([ str(Component) , str(High_Low) , 'Component-Blood'           ,  str(  Compo_VS_CFEA____dict['mean:distance'])  ,  str(  Compo_VS_CFEA____dict['mean:logFC'])  ,  str(  Compo_VS_CFEA____dict['median:distance'])  ,  str(  Compo_VS_CFEA____dict['median:logFC'])  ,  str(  Compo_VS_CFEA____dict['across:distance'])  ,  str(  Compo_VS_CFEA____dict['across:logFC'])   ])  +'\n'  )
                        outputTEXT_MCom.write(  ('\t').join([ str(Component) , str(High_Low) , 'Bulked:Sorted_by_AUC'      ,  str(     Bulked_AUC____dict['mean:distance'])  ,  str(     Bulked_AUC____dict['mean:logFC'])  ,  str(     Bulked_AUC____dict['median:distance'])  ,  str(     Bulked_AUC____dict['median:logFC'])  ,  str(     Bulked_AUC____dict['across:distance'])  ,  str(     Bulked_AUC____dict['across:logFC'])   ])  +'\n'  )
                        outputTEXT_MCom.write(  ('\t').join([ str(Component) , str(High_Low) , 'Bulked:Sorted_by_Distance' ,  str(Bulked_Distance____dict['mean:distance'])  ,  str(Bulked_Distance____dict['mean:logFC'])  ,  str(Bulked_Distance____dict['median:distance'])  ,  str(Bulked_Distance____dict['median:logFC'])  ,  str(Bulked_Distance____dict['across:distance'])  ,  str(Bulked_Distance____dict['across:logFC'])   ])  +'\n'  )
                        outputTEXT_MCom.write(  ('\t').join([ str(Component) , str(High_Low) , 'Bulked:Sorted_by_logFC'    ,  str(   Bulked_logFC____dict['mean:distance'])  ,  str(   Bulked_logFC____dict['mean:logFC'])  ,  str(   Bulked_logFC____dict['median:distance'])  ,  str(   Bulked_logFC____dict['median:logFC'])  ,  str(   Bulked_logFC____dict['across:distance'])  ,  str(   Bulked_logFC____dict['across:logFC'])   ])  +'\n'  )
                        outputTEXT_MCom.write(  ('\t').join([ str(Component) , str(High_Low) , 'Bulked:Sorted_by_U-test'   ,  str(   Bulked_Utest____dict['mean:distance'])  ,  str(   Bulked_Utest____dict['mean:logFC'])  ,  str(   Bulked_Utest____dict['median:distance'])  ,  str(   Bulked_Utest____dict['median:logFC'])  ,  str(   Bulked_Utest____dict['across:distance'])  ,  str(   Bulked_Utest____dict['across:logFC'])   ])  +'\n'  )
                        ##
                        ##
                        ##
                        CDT_DF[ str(Component)+' ---- Compo-Blood'     ].loc[[str(High_Low)+':'+str(i+1) for i in range(parameters['Top'])]]  =  list(Compo_VS_CFEA)   
                        CDT_DF[ str(Component)+' ---- Bulked_AUC'      ].loc[[str(High_Low)+':'+str(i+1) for i in range(parameters['Top'])]]  =  list(Bulked_AUC)      
                        CDT_DF[ str(Component)+' ---- Bulked_Distance' ].loc[[str(High_Low)+':'+str(i+1) for i in range(parameters['Top'])]]  =  list(Bulked_Distance) 
                        CDT_DF[ str(Component)+' ---- Bulked_logF'     ].loc[[str(High_Low)+':'+str(i+1) for i in range(parameters['Top'])]]  =  list(Bulked_logFC)    
                        CDT_DF[ str(Component)+' ---- Bulked_Utest'    ].loc[[str(High_Low)+':'+str(i+1) for i in range(parameters['Top'])]]  =  list(Bulked_Utest)    
                        CDT_DF[ str(Component)+' ---- Gap'             ].loc[[str(High_Low)+':'+str(i+1) for i in range(parameters['Top'])]]  =  numpy.nan
                        ##print(  CDT_DF[ [ str(Component)+' ---- Compo-Blood' , str(Component)+' ---- Bulked_AUC' , str(Component)+' ---- Bulked_Distance' , str(Component)+' ---- Bulked_logF' , str(Component)+' ---- Bulked_Utest' , str(Component)+' ---- Gap'] ] )
                        ##
                        '''
                        Data_Series  =  pandas.DataFrame(
                            pandas.concat( [ Compo_VS_CFEA , Bulked_AUC , Bulked_Distance , Bulked_logFC , Bulked_Utest ] , axis=0  )  
                        ).rename( columns={ parameters['CFEA_data_type']+'|'+parameters['CFEA_X']+'|Ratio_across_samples:distance':  'Cancer/normal' })
                        ##
                        Group_Series  =  pandas.DataFrame(
                            ['Compo-Blood']*parameters['Top']  +  ['Bulked\nAUC']*parameters['Top']  +  ['Bulked\nDistance']*parameters['Top']  +  ['Bulked\nlogFC']*parameters['Top']  +  ['Bulked\nUtest']*parameters['Top'] ,
                            index=Data_Series.index , 
                            columns=['Component/bulked']   
                        )
                        ##
                        Top_DF  =  pandas.concat( [ Data_Series , Group_Series ] , axis=1 )
                        print( Top_DF )
                        '''
                        ##
                        ##
                        ##
                        '''
                        matplotlib.pyplot.figure( )
                        matplotlib.pyplot.clf( )
                        seaborn.swarmplot(        data=Top_DF  ,  x='Component/bulked'  ,  y='Cancer/normal'  ,  size=1   )
                        ax = seaborn.violinplot(  data=Top_DF  ,  x='Component/bulked'  ,  y='Cancer/normal'  ,  color='white'  ,  linewidth=0.5  ,  cut=0  )
                        ax.collections[0].set_edgecolor('black')
                        ax.collections[1].set_edgecolor('black')
                        ax.collections[2].set_edgecolor('black')
                        ax.collections[3].set_edgecolor('black')
                        ax.collections[4].set_edgecolor('black')
                        ax.collections[5].set_edgecolor('black')
                        matplotlib.pyplot.title(  str(Component)+'|'+str(High_Low)  )
                        PDF.savefig()
                        '''
                        ##
                        del InputFile , Compo_DF 
                        del Compo_VS_CFEA          , Bulked_AUC             , Bulked_Distance          , Bulked_logFC         , Bulked_Utest
                        del Compo_VS_CFEA____list  , Bulked_AUC____list     , Bulked_Distance____list  , Bulked_logFC____list , Bulked_Utest____list
                        del Compo_VS_CFEA____dict  , Bulked_AUC____dict     , Bulked_Distance____dict  , Bulked_logFC____dict , Bulked_Utest____dict
                        del Compo_Mean , Compo_STD , B_AUC_Mean , B_AUC_STD , B_Dist_Mean , B_Dist_STD , B_FC_Mean , B_FC_STD , B_U_Mean , B_U_STD  
                        del Compo_VS_B_AUC____Ttest , Compo_VS_B_AUC____Utest , Compo_VS_B_Dist____Ttest , Compo_VS_B_Dist____Utest , Compo_VS_B_logFC____Ttest , Compo_VS_B_logFC____Utest , Compo_VS_B_Utest____Ttest , Compo_VS_B_Utest____Utest
                        ####del Data_Series , Group_Series , Top_DF
                ##
                ###PDF.close()
                del Input_Sub_Path , Output_PDF
                outputTEXT_stat.close()
                outputTEXT_MCom.close()
                ##
                ##
                ##
                Output_TSV = str(OutputPath)+'NMF_'+str(K)+'.tsv'
                CDT_DF.to_csv( Output_TSV , sep='\t' )
                del CDT_DF 
                ##
                ##
                ##
                Output_CDT = str(OutputPath)+'NMF_'+str(K)+'.cdt'
                inputTEXT  = open( Output_TSV , 'r' )    #### 只读
                outputTEXT = open( Output_CDT , 'w' )
                headline_list = next( inputTEXT ).split('\n')[0].split('\r')[0].split('\t')[ 1 :   ]
                outputTEXT.write(  ('\t').join( ['ID' , 'NAME' , 'GWEIGHT'] +           headline_list    )  +'\n'  )
                outputTEXT.write(  ('\t').join( ['EWEIGHT'  ,  ''  ,  ''  ] + ['1']*len(headline_list)   )  +'\n'  )
                while True:
                    try:
                        line2list = next( inputTEXT ).split('\n')[0].split('\r')[0].split('\t')
                        outputTEXT.write( ('\t').join( [ str(line2list[0]) , '' , '1' ]  +  line2list[ 1 :   ] )    +'\n'  )
                    except StopIteration:
                        break
                inputTEXT.close()    #### 读毕
                outputTEXT.close()
                del Output_TSV , Output_CDT
            ##
            ##
            ##
            del InputPath , OutputPath
            del K____list , K_compo____dict
        ##
    ##
    print('\n')

















if __name__ == '__main__':
    print( '\n' )

    #### 开始执行函数 ####

    ##main_loop( parameters_Colon____colorectal_1X  )
    ##main_loop( parameters_Colon____colorectal_10X  )

    main_loop( parameters_Colon____blood_10X____colorectal_1X  )
    main_loop( parameters_Colon____blood_10X____colorectal_10X  )



    ##main_loop( parameters_Lung____lung_1X   )
    ##main_loop( parameters_Lung____lung_10X   )

    main_loop( parameters_Lung____blood_10X____lung_1X   )
    main_loop( parameters_Lung____blood_10X____lung_10X   )
    
    #### 执行函数完毕 ####
    
    print( '\nThe end.\n' )










