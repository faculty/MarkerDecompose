#!/usr/bin/python3
# -*- coding: utf-8 -*-




import sys
from optparse import OptionParser
import os
import numpy
import glob
import pandas








parameters_Breast = {
    'Input_Sample':      '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Breast____.sample.tsv',
    'InputPath':         '/Path/to/Analysis/31_01____NMF_with_InternalReference/Breast/',
    'OutputPath':        '/Path/to/Analysis/31_02____Components/Breast/',
    'Input_Coefficient': None, 
    'OutputFile':        None, 
}

parameters_Colon = {
    'Input_Sample':      '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Colon____.sample.tsv',
    'InputPath':         '/Path/to/Analysis/31_01____NMF_with_InternalReference/Colon/',
    'OutputPath':        '/Path/to/Analysis/31_02____Components/Colon/',
    'Input_Coefficient': None, 
    'OutputFile':        None, 
}

parameters_Eso = {
    'Input_Sample':      '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Eso____.sample.tsv',
    'InputPath':         '/Path/to/Analysis/31_01____NMF_with_InternalReference/Eso/',
    'OutputPath':        '/Path/to/Analysis/31_02____Components/Eso/',
    'Input_Coefficient': None, 
    'OutputFile':        None, 
}

parameters_Gastric = {
    'Input_Sample':      '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Gastric____.sample.tsv',
    'InputPath':         '/Path/to/Analysis/31_01____NMF_with_InternalReference/Gastric/',
    'OutputPath':        '/Path/to/Analysis/31_02____Components/Gastric/',
    'Input_Coefficient': None, 
    'OutputFile':        None, 
}

parameters_Liver = {
    'Input_Sample':      '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Liver____.sample.tsv',
    'InputPath':         '/Path/to/Analysis/31_01____NMF_with_InternalReference/Liver/',
    'OutputPath':        '/Path/to/Analysis/31_02____Components/Liver/',
    'Input_Coefficient': None, 
    'OutputFile':        None, 
}

parameters_Lung = {
    'Input_Sample':      '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Lung____.sample.tsv',
    'InputPath':         '/Path/to/Analysis/31_01____NMF_with_InternalReference/Lung/',
    'OutputPath':        '/Path/to/Analysis/31_02____Components/Lung/',
    'Input_Coefficient': None, 
    'OutputFile':        None, 
}

parameters_Pancreas = {
    'Input_Sample':      '/Path/to/Analysis/21_01____RRBS_Matrix_of_CpG_by_TissueType/____Pancreas____.sample.tsv',
    'InputPath':         '/Path/to/Analysis/31_01____NMF_with_InternalReference/Pancreas/',
    'OutputPath':        '/Path/to/Analysis/31_02____Components/Pancreas/',
    'Input_Coefficient': None, 
    'OutputFile':        None, 
}





class main_loop:
    def __init__( self , parameters ):
        self.Read_Coefficients( parameters )
        self.Read_Sample( parameters )
        self.Write_down_matrix( parameters )
        self.Make_CDT( parameters )
        print('\n')
    

    def Read_Coefficients( self , parameters ):
        print('Read from: ' , parameters['Input_Coefficient'] )
        self.Coefficient_DF  =  pandas.read_csv(  parameters['Input_Coefficient']  ,  header=0  ,  index_col=0  ,  sep='\t'  )
        print(  self.Coefficient_DF.shape  )
    

    def Read_Sample( self , parameters ):
        self.Pathology_Sample____dict  =  {  }
        self.Sample_Pathology____dict  =  {  }
        self.Sample____list  =  [  ]
        ########
        count_line = 0         #### 计数
        print('Read from: ' , parameters['Input_Sample']  )
        inputTEXT = open( parameters['Input_Sample'] , 'r' )    #### 只读
        next( inputTEXT )    #### headline
        while True:
            try:
                line2list = next( inputTEXT ).split('\n')[0].split('\r')[0].split('\t')
                count_line = count_line + 1    #### 计数
                ##
                Sample    = str(line2list[0])
                Pathology = str(line2list[1])
                ##
                if Pathology not in self.Pathology_Sample____dict.keys():
                    self.Pathology_Sample____dict[  Pathology  ]  =  [  ]
                else:
                    pass
                ##
                self.Pathology_Sample____dict[  Pathology  ].append(  Sample  )
                ##
                self.Sample_Pathology____dict[  Sample  ]  =  Pathology
                ##
                self.Sample____list.append(  Sample  )
            except StopIteration:
                break
        inputTEXT.close()    #### 读毕
        print('\tChecking ==> Lines: ' , count_line  )
        print('\tChecking ==> self.Pathology_Sample____dict: ' , sorted( self.Pathology_Sample____dict.keys() )    )
        print('\tChecking ==> self.Sample_Pathology____dict: ' ,    len( self.Sample_Pathology____dict.keys() )    )
        print('\tChecking ==> self.Sample____list: ' , len(self.Sample____list)  )


    def Write_down_matrix( self , parameters ):
        Ranked_Sample_list = [  ]
        Sample_SampleNewName____dict = {  }
        ########
        print('Sorting sample')
        for Sample in self.Sample____list:
            Pathology  =  self.Sample_Pathology____dict[  Sample  ]
            Ranked_Sample_list.append(  Sample  )
            Sample_SampleNewName____dict[  Sample  ]  =  str(Pathology)+' ---- '+str(Sample)
        ##
        print('\tChecking ==> Ranked_Sample_list: ' , len(Ranked_Sample_list)  )
        print('\tChecking ==> Sample_SampleNewName____dict: ' , len(Sample_SampleNewName____dict.keys())  )
        ##
        ##
        ##
        print('Rename column header')
        written_DF  =  self.Coefficient_DF[  Ranked_Sample_list  ]
        written_DF.rename(  columns=Sample_SampleNewName____dict  ,  inplace=True  )
        print(  written_DF.shape  )
        ##
        ##
        ##
        print('Write down: ' , parameters['OutputFile'] )
        written_DF.to_csv(  parameters['OutputFile']  ,  sep='\t'  )
    

    def Make_CDT( self , parameters ):
        print('Read from:  ' , parameters['OutputFile'] )
        print('Write down: ' , parameters['Output_CDT'] )
        inputTEXT  = open( parameters['OutputFile'] , 'r' )    #### 只读
        outputTEXT = open( parameters['Output_CDT'] , 'w' )
        ##
        Sample_list = next(inputTEXT).split('\n')[0].split('\r')[0].split('\t')[ 1 :   ]
        outputTEXT.write(  ('\t').join( ['ID'     ,'NAME','GWEIGHT']  +            Sample_list[:]  )  +'\n'  )
        outputTEXT.write(  ('\t').join( ['EWEIGHT',  ''  ,   ''    ]  +  ['1']*len(Sample_list[:]) )  +'\n'  )
        ##
        while True:
            try:
                line2list = next(inputTEXT).split('\n')[0].split('\r')[0].split('\t')
                Component = str(line2list[0])
                data_list =     line2list[ 1 :   ]
                ##
                if len(Sample_list) == len(data_list):
                    pass
                else:
                    print('Error ==> ' , len(Sample_list) , len(data_list)  )
                ##
                outputTEXT.write(  ('\t').join( [ str(Component) , str(Component) , '1' ]  +  data_list[:] )  +'\n'  )
                ##
                del Component , data_list
            except StopIteration:
                break
        inputTEXT.close()    #### 读毕
        outputTEXT.close()












if __name__ == '__main__':
    print('\n')


    #### 开始执行函数 ####
    '''
    ##for i in [ 2 , 3 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ]:
    for i in [ 6 , 7 , 8 , 9 ]:
        parameters_Breast['Input_Coefficient']  =  parameters_Breast['InputPath'] +str(i)+'____scaled_W____Components_Sample.tsv'
        parameters_Breast['OutputFile']         =  parameters_Breast['OutputPath']+str(i)+'____scaled_W____Components_Sample.txt'
        parameters_Breast['Output_CDT']         =  parameters_Breast['OutputPath']+str(i)+'____scaled_W____Components_Sample.cdt'
        main_loop(  parameters_Breast  )
        del parameters_Breast['Input_Coefficient']
        del parameters_Breast['OutputFile']
        del parameters_Breast['Output_CDT']

    ##for i in [ 2 , 3 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ]:
    for i in [ 6 , 7 , 8 , 9 ]:
        parameters_Colon['Input_Coefficient']  =  parameters_Colon['InputPath'] +str(i)+'____scaled_W____Components_Sample.tsv'
        parameters_Colon['OutputFile']         =  parameters_Colon['OutputPath']+str(i)+'____scaled_W____Components_Sample.txt'
        parameters_Colon['Output_CDT']         =  parameters_Colon['OutputPath']+str(i)+'____scaled_W____Components_Sample.cdt'
        main_loop(  parameters_Colon  )
        del parameters_Colon['Input_Coefficient']
        del parameters_Colon['OutputFile']
        del parameters_Colon['Output_CDT']
    '''
    


    '''
    for i in [ 2 , 3 , 5 , 6 , 7 , 8 , 9 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ]:
        parameters_Eso['Input_Coefficient']  =  parameters_Eso['InputPath'] +str(i)+'____scaled_W____Components_Sample.tsv'
        parameters_Eso['OutputFile']         =  parameters_Eso['OutputPath']+str(i)+'____scaled_W____Components_Sample.txt'
        parameters_Eso['Output_CDT']         =  parameters_Eso['OutputPath']+str(i)+'____scaled_W____Components_Sample.cdt'
        main_loop(  parameters_Eso  )
        del parameters_Eso['Input_Coefficient']
        del parameters_Eso['OutputFile']
        del parameters_Eso['Output_CDT']
    '''



    
    for i in [ 2 , 3 , 5 , 6 , 7 , 8 , 9 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ]:
        parameters_Gastric['Input_Coefficient']  =  parameters_Gastric['InputPath'] +str(i)+'____scaled_W____Components_Sample.tsv'
        parameters_Gastric['OutputFile']         =  parameters_Gastric['OutputPath']+str(i)+'____scaled_W____Components_Sample.txt'
        parameters_Gastric['Output_CDT']         =  parameters_Gastric['OutputPath']+str(i)+'____scaled_W____Components_Sample.cdt'
        main_loop(  parameters_Gastric  )
        del parameters_Gastric['Input_Coefficient']
        del parameters_Gastric['OutputFile']
        del parameters_Gastric['Output_CDT']
    




    '''
    ##for i in [ 2 , 3 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ]:
    for i in [ 6 , 7 , 8 , 9 ]:
        parameters_Liver['Input_Coefficient']  =  parameters_Liver['InputPath'] +str(i)+'____scaled_W____Components_Sample.tsv'
        parameters_Liver['OutputFile']         =  parameters_Liver['OutputPath']+str(i)+'____scaled_W____Components_Sample.txt'
        parameters_Liver['Output_CDT']         =  parameters_Liver['OutputPath']+str(i)+'____scaled_W____Components_Sample.cdt'
        main_loop(  parameters_Liver  )
        del parameters_Liver['Input_Coefficient']
        del parameters_Liver['OutputFile']
        del parameters_Liver['Output_CDT']
    

    ##for i in [ 2 , 3 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ]:
    for i in [ 6 , 7 , 8 , 9 ]:
        parameters_Lung['Input_Coefficient']  =  parameters_Lung['InputPath'] +str(i)+'____scaled_W____Components_Sample.tsv'
        parameters_Lung['OutputFile']         =  parameters_Lung['OutputPath']+str(i)+'____scaled_W____Components_Sample.txt'
        parameters_Lung['Output_CDT']         =  parameters_Lung['OutputPath']+str(i)+'____scaled_W____Components_Sample.cdt'
        main_loop(  parameters_Lung  )
        del parameters_Lung['Input_Coefficient']
        del parameters_Lung['OutputFile']
        del parameters_Lung['Output_CDT']
    

    ##for i in [ 2 , 3 , 5 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ]:
    for i in [ 6 , 7 , 8 , 9 ]:
        parameters_Pancreas['Input_Coefficient']  =  parameters_Pancreas['InputPath'] +str(i)+'____scaled_W____Components_Sample.tsv'
        parameters_Pancreas['OutputFile']         =  parameters_Pancreas['OutputPath']+str(i)+'____scaled_W____Components_Sample.txt'
        parameters_Pancreas['Output_CDT']         =  parameters_Pancreas['OutputPath']+str(i)+'____scaled_W____Components_Sample.cdt'
        main_loop(  parameters_Pancreas  )
        del parameters_Pancreas['Input_Coefficient']
        del parameters_Pancreas['OutputFile']
        del parameters_Pancreas['Output_CDT']
    '''


    #### 执行函数完毕 ####
    

    print('\nThe end.\n')









