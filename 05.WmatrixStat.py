#!/usr/bin/python3
# -*- coding: utf-8 -*-




import sys
from optparse import OptionParser
import os
import glob
import math
import numpy
import pandas
import scipy.stats










parameters  =  {
    'Type': [
        'Breast:Normal', 'Breast:Cancer', 'Colon:Normal', 'Colon:Cancer', 'Eso:Normal'     , 'Eso:Cancer'     , 'Gastric:Normal', 'Gastric:Cancer', 
        'Lung:Normal'  , 'Lung:Cancer'  , 'Liver:Normal', 'Liver:Cancer', 'Pancreas:Normal', 'Pancreas:Cancer', 
        'Plasma:Normal', 'Blood:Normal' ,
    ],
    'AvailableType': [ 
        'Breast:Cancer', 'Colon:Cancer', 'Eso:Cancer', 'Gastric:Cancer', 'Lung:Cancer', 'Liver:Cancer', 'Pancreas:Cancer',
        'Plasma:Normal', 'Blood:Normal', 
    ],
    'CancerType': [ 
        'Breast:Cancer', 'Colon:Cancer', 'Eso:Cancer', 'Gastric:Cancer', 'Lung:Cancer', 'Liver:Cancer', 'Pancreas:Cancer',
    ],
    'K':          ( 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ),
    'WorkPath':   '/Path/to/Analysis/51_03____Components/',
}





class main_loop:
    def __init__( self , parameters ):
        for K in parameters['K']:
            InputFile  =  str(parameters['WorkPath']) +str(K)+'____scaled_W____Components_Sample.txt'
            Input_CDT  =  str(parameters['WorkPath']) +str(K)+'____scaled_W____Components_Sample.cdt'
            OutputFile =  str(parameters['WorkPath']) +str(K)+'____scaled_W____Components_Sample.Statistics____Fraction_of_Cancer.tsv'
            OutputRepo =  str(parameters['WorkPath']) +str(K)+'____scaled_W____Components_Sample.StatSummary____Fraction_of_Cancer.tsv'
            ##
            self.Scan_Fraction(  InputFile  ,  Input_CDT  ,  OutputFile  ,  OutputRepo  )
            ##
            del  InputFile  ,  Input_CDT  ,  OutputFile  ,  OutputRepo
    

    def Scan_Fraction(  self  ,  InputFile  ,  Input_CDT  ,  OutputFile  ,  OutputRepo  ):
        print('Read fraction matrix from: ' , InputFile  )
        Fraction_DF  =  pandas.read_csv(  InputFile  ,  header=0  ,  index_col=0  ,  sep='\t'  )
        print( Fraction_DF.shape )
        ##
        ##
        ##
        Sample_dict = {  }
        for AvailableType in parameters['AvailableType']:
            Sample_dict[ AvailableType ]  =  [  ]
        for Sample in Fraction_DF.columns.values:
            Type = Sample.split(' ---- ')[0]
            if Type in Sample_dict.keys():
                Sample_dict[ Type ].append( Sample )
            else:
                pass
        for AvailableType in parameters['AvailableType']:
            print( AvailableType , '  ==>  ' , len(Sample_dict[ AvailableType ])  )
        ##
        ##
        ##
        print('Do statistics ...... ')
        Stat_DF = pandas.DataFrame( columns=[
            'Average|Breast:Cancer', 'Average|Colon:Cancer', 'Average|Eso:Cancer'     , 'Average|Gastric:Cancer', 
            'Average|Lung:Cancer'  , 'Average|Liver:Cancer', 'Average|Pancreas:Cancer',
            'Average|Plasma:Normal', 'Average|Blood:Normal',
            'Breast----Colon,Eso,Gastric,Lung,Liver,Pancreas(Cancer);Plasma,Blood(Normal)',
            'Colon----Breast,Eso,Gastric,Lung,Liver,Pancreas(Cancer);Plasma,Blood(Normal)',
            'Eso----Breast,Colon,Gastric,Lung,Liver,Pancreas(Cancer);Plasma,Blood(Normal)',
            'Gastric----Breast,Colon,Eso,Lung,Liver,Pancreas(Cancer);Plasma,Blood(Normal)',
            'Lung----Breast,Colon,Eso,Gastric,Liver,Pancreas(Cancer);Plasma,Blood(Normal)',
            'Liver----Breast,Colon,Eso,Gastric,Lung,Pancreas(Cancer);Plasma,Blood(Normal)',
            'Pancreas----Breast,Colon,Eso,Gastric,Lung,Liver(Cancer);Plasma,Blood(Normal)',
            'Breast(max_P)', 'Colon(max_P)', 'Eso(max_P)', 'Gastric(max_P)', 'Lung(max_P)', 'Liver(max_P)', 'Pancreas(max_P)',
        ] )
        for Component in Fraction_DF.index.values:
            ComponentStat_dict = {  }
            ComponentStat_dict['Average|Breast:Cancer']    =  str( Fraction_DF.loc[ Component ][ Sample_dict['Breast:Cancer'  ] ].mean() )
            ComponentStat_dict['Average|Colon:Cancer']     =  str( Fraction_DF.loc[ Component ][ Sample_dict['Colon:Cancer'   ] ].mean() )
            ComponentStat_dict['Average|Eso:Cancer']       =  str( Fraction_DF.loc[ Component ][ Sample_dict['Eso:Cancer'     ] ].mean() )
            ComponentStat_dict['Average|Gastric:Cancer']   =  str( Fraction_DF.loc[ Component ][ Sample_dict['Gastric:Cancer' ] ].mean() )
            ComponentStat_dict['Average|Lung:Cancer']      =  str( Fraction_DF.loc[ Component ][ Sample_dict['Lung:Cancer'    ] ].mean() )
            ComponentStat_dict['Average|Liver:Cancer']     =  str( Fraction_DF.loc[ Component ][ Sample_dict['Liver:Cancer'   ] ].mean() )
            ComponentStat_dict['Average|Pancreas:Cancer']  =  str( Fraction_DF.loc[ Component ][ Sample_dict['Pancreas:Cancer'] ].mean() )
            ComponentStat_dict['Average|Plasma:Normal']    =  str( Fraction_DF.loc[ Component ][ Sample_dict['Plasma:Normal'  ] ].mean() )
            ComponentStat_dict['Average|Blood:Normal']     =  str( Fraction_DF.loc[ Component ][ Sample_dict['Blood:Normal'   ] ].mean() )
            ##
            for OneType in parameters['CancerType']:
                OneType_sample     =  list( Sample_dict[ OneType ] )
                OneType_data_list  =  list([  float(X)  for  X  in  Fraction_DF.loc[ Component ][ OneType_sample ]  ])
                ##
                Pvalue_list = [  ]
                ##
                for AnotherType in parameters['AvailableType']:
                    AnotherType_sample     =  list( Sample_dict[ AnotherType ] )
                    AnotherType_data_list  =  list([  float(X)  for  X  in  Fraction_DF.loc[ Component ][ AnotherType_sample ]  ])
                    ##
                    if OneType == AnotherType:
                        if (OneType_sample==AnotherType_sample)  or  (OneType_data_list==AnotherType_data_list):
                            pass
                        else:
                            print('Error  ==>  type and data ')
                    else:
                        Pvalue  =  scipy.stats.ttest_ind( OneType_data_list  , AnotherType_data_list , alternative='greater' )[1]
                        Pvalue_list.append(  float(Pvalue)  )
                    ##
                    del AnotherType_sample , AnotherType_data_list
                ##
                Pvalue_output  =  (';').join([ str(X)  for  X  in  Pvalue_list ])
                ##
                if   OneType=='Breast:Cancer':
                    ComponentStat_dict[ 'Breast----Colon,Eso,Gastric,Lung,Liver,Pancreas(Cancer);Plasma,Blood(Normal)' ] = str(Pvalue_output)
                    ComponentStat_dict[ 'Breast(max_P)' ]  =  max(Pvalue_list)
                elif OneType=='Colon:Cancer':
                    ComponentStat_dict[ 'Colon----Breast,Eso,Gastric,Lung,Liver,Pancreas(Cancer);Plasma,Blood(Normal)' ] = str(Pvalue_output)
                    ComponentStat_dict[ 'Colon(max_P)' ]  =  max(Pvalue_list)
                elif OneType=='Eso:Cancer':
                    ComponentStat_dict[ 'Eso----Breast,Colon,Gastric,Lung,Liver,Pancreas(Cancer);Plasma,Blood(Normal)' ] = str(Pvalue_output)
                    ComponentStat_dict[ 'Eso(max_P)' ]  =  max(Pvalue_list)
                elif OneType=='Gastric:Cancer':
                    ComponentStat_dict[ 'Gastric----Breast,Colon,Eso,Lung,Liver,Pancreas(Cancer);Plasma,Blood(Normal)' ] = str(Pvalue_output)
                    ComponentStat_dict[ 'Gastric(max_P)' ]  =  max(Pvalue_list)
                elif OneType=='Lung:Cancer':
                    ComponentStat_dict[ 'Lung----Breast,Colon,Eso,Gastric,Liver,Pancreas(Cancer);Plasma,Blood(Normal)' ] = str(Pvalue_output)
                    ComponentStat_dict[ 'Lung(max_P)' ]  =  max(Pvalue_list)
                elif OneType=='Liver:Cancer':
                    ComponentStat_dict[ 'Liver----Breast,Colon,Eso,Gastric,Lung,Pancreas(Cancer);Plasma,Blood(Normal)' ] = str(Pvalue_output)
                    ComponentStat_dict[ 'Liver(max_P)' ]  =  max(Pvalue_list)
                elif OneType=='Pancreas:Cancer':
                    ComponentStat_dict[ 'Pancreas----Breast,Colon,Eso,Gastric,Lung,Liver(Cancer);Plasma,Blood(Normal)' ] = str(Pvalue_output)
                    ComponentStat_dict[ 'Pancreas(max_P)' ]  =  max(Pvalue_list)
                else:
                    print('Error ==> One type error ' , OneType  )
                ##
                del OneType_sample , OneType_data_list , Pvalue_list , Pvalue_output 
            ##
            Stat_DF.loc[ Component ]  =  pandas.Series(  ComponentStat_dict  )
        ##
        print(  Stat_DF.shape  )
        ##
        ##
        ##
        Component_sortedList = [  ]
        print('Read cluster CDT from: ' , Input_CDT  )
        inputTEXT = open( Input_CDT , 'r' )    #### 只读
        next( inputTEXT )    #### headline 1
        next( inputTEXT )    #### headline 2
        while True:
            try:
                line2list = next( inputTEXT ).split('\n')[0].split('\r')[0].split('\t')
                Component = str(line2list[1])
                Component_sortedList.append(  Component  )
            except StopIteration:
                break
        inputTEXT.close()    #### 读毕
        print('\tChecking ==> Component_sortedList: ' , len(Component_sortedList)  )
        ##
        ##
        ##
        print('Write down: ' , OutputFile  )
        Stat_DF.loc[ Component_sortedList ].to_csv(  OutputFile , sep='\t' )
        ##
        ##
        ##
        print('Write down: ' , OutputRepo  )
        outputTEXT = open( OutputRepo , 'w' )
        outputTEXT.write(  ('\t').join(['CancerType' , 'Min_Pvalue' , 'Component' , 'Fraction' , 'logP*Fraction'])  +'\n'  )
        for CancerType in parameters['CancerType']:
            TOO = CancerType.split(':')[0]
            Min_Pvalue  =  Stat_DF[ str(TOO)+'(max_P)' ].min()
            Component   =  (',').join( Stat_DF[    Stat_DF[ str(TOO)+'(max_P)' ]  ==  Stat_DF[ str(TOO)+'(max_P)' ].min()    ].index.values )
            Fraction    =  Stat_DF.loc[ Component ][ 'Average|'+str(TOO)+':Cancer']
            logP_Frac   =  -math.log2(Min_Pvalue) * float(Fraction)
            ##
            outputTEXT.write(  ('\t').join([ str(CancerType) , str(Min_Pvalue) , str(Component) , str(Fraction) , str(logP_Frac) ])  +'\n'  )
            ##
            del  TOO , Min_Pvalue , Component , Fraction , logP_Frac
        outputTEXT.close()
        ##
        ##
        ##
        print('\n')











if __name__ == '__main__':
    print( '\n' )

    #### 开始执行函数 ####

    main_loop( parameters )
    
    #### 执行函数完毕 ####
    
    print( '\nThe end.\n' )










