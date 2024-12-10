#!/usr/bin/python3
# -*- coding: utf-8 -*-




import sys
from optparse import OptionParser
import os
import numpy
import glob
import pandas








parameters = {
    'InputPath':         '/Path/to/Analysis/51_02____NMF_with_InternalReference/NMF/',
    'OutputPath':        '/Path/to/Analysis/51_03____Components/',
    'InputFile':         None, 
    'OutputFile':        None, 
}











class main_loop:
    def __init__( self , parameters ):
        self.Read_H_Files( parameters )
        self.Write_down_matrix( parameters )
        print('\n')
    

    def Read_H_Files( self , parameters ):
        print('Read from: ' , parameters['InputFile'] )
        self.H_DF  =  pandas.read_csv(  parameters['InputFile']  ,  header=0  ,  index_col=0  ,  sep='\t'  )
        print(  self.H_DF.shape  )
        ##
        print('Remove internal reference ...... ')
        self.H_DF  =  self.H_DF.drop( ['InternalReference:Meth=1','InternalReference:Meth=0']  ,  axis=0 )
        print(  self.H_DF.shape  )
        QuantityBeforeAdj  =  numpy.linalg.norm( self.H_DF )
        ##
        ##
        ##
        print('Adjust sup outlier above:' , 100.0  )
        self.H_DF[  self.H_DF > 100.0  ]    =    100.0
        ##
        print('Adjust inf outlier below:' ,   0.0  )
        self.H_DF[  self.H_DF <   0.0  ]    =      0.0
        ##
        QuantityAfterAdj   =  numpy.linalg.norm( self.H_DF )
        print('Before truncation: ' , QuantityBeforeAdj , '    After truncation: ' , QuantityAfterAdj , '    ==>    ' , QuantityAfterAdj/QuantityBeforeAdj  )


    def Write_down_matrix( self , parameters ):
        print('Write down: ' , parameters['OutputFile'] )
        self.H_DF.to_csv(  parameters['OutputFile']  ,  sep='\t'  )












if __name__ == '__main__':
    print('\n')


    #### 开始执行函数 ####

    for i in [ 5 , 6 , 7 , 8 , 9 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60 , 70 , 80 , 90 , 100 ]:
        parameters['InputFile']  =  parameters['InputPath'] +str(i)+'____scaled_H____Marker_Components.tsv'
        parameters['OutputFile'] =  parameters['OutputPath']+str(i)+'____scaled_H____Marker_Components.Truncated_IntRefRemoved.txt'
        main_loop(  parameters  )
        del parameters['InputFile']
        del parameters['OutputFile']
    

    #### 执行函数完毕 ####
    

    print('\nThe end.\n')









