# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:48:03 2022

@author: lawashburn
"""

import pandas as pd
import numpy as np
import csv
from scipy import stats
import math

working_directory = input('Input path to working directory, where all results will be stored: ')
protein_path = input('Input path to Proteingroup.csv file: ')
peptides_path = input('Input path to Peptidegroup.csv file: ')
#protein_path = r"C:\Users\lawashburn\Documents\LiP-MS_Haiyan\20220527\proteinGroups_abridged.csv"
#peptides_path = r"C:\Users\lawashburn\Documents\LiP-MS_Haiyan\20220527\peptides_abridged.csv"
#working_directory = r"C:\Users\lawashburn\Documents\LiP-MS_Haiyan\20220527\test_out"


protein = pd.read_csv(protein_path)
peptide = pd.read_csv(peptides_path)
protein = protein.rename(columns={'LFQ intensity SRT_L1_T':'SRT1','LFQ intensity SRT_L2_T':'SRT2','LFQ intensity SRT_L3_T':'SRT3',
                                  'LFQ intensity W_L1_T':'Ctrl1',
                                  'LFQ intensity W_L2_T':'Ctrl2','LFQ intensity W_L3_T':'Ctrl3'})

protein_out2 = working_directory + '\\protein_test.csv'
with open(protein_out2,'w',newline='') as filed:
    writerd = csv.writer(filed)
    protein.to_csv(filed,index=False)

SRT = pd.DataFrame()
SRT['SRT1'] = protein['SRT1']
SRT['SRT2'] = protein['SRT2']
SRT['SRT3'] = protein['SRT3']
SRT = SRT.iloc[1:,:]
SRT['mean'] = SRT.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

Ctrl = pd.DataFrame()
Ctrl['Ctrl1'] = protein['Ctrl1']
Ctrl['Ctrl2'] = protein['Ctrl2']
Ctrl['Ctrl3'] = protein['Ctrl3']
Ctrl = Ctrl.iloc[1:,:]
Ctrl['mean'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

protein_calc = pd.DataFrame()
protein_calc['Protein ID'] = protein['Protein IDs']
protein_calc['Protein names'] = protein['Protein names']
protein_calc['Gene name'] = protein['Gene names']
protein_calc['Fasta header'] = protein['Fasta headers']
protein_calc['SRT avg'] = SRT['mean']
protein_calc['Ctrl avg'] = Ctrl['mean']
protein_calc['SRT/Ctrl fold change'] = protein_calc['SRT avg'] / protein_calc['Ctrl avg']

df_common = pd.merge(protein_calc, peptide, on=['Protein names'], how='inner')

SRTvCtrl_out2 = working_directory + '\\common.csv'
with open(SRTvCtrl_out2,'w',newline='') as filed:
    writerd = csv.writer(filed)
    df_common.to_csv(filed,index=False)


df_common.dropna(subset = ["Protein names"], inplace=True)

df_common = df_common.rename(columns={'LFQ intensity W_L1_P':'Pep Ctrl 1','LFQ intensity W_L2_P':'Pep Ctrl 2','LFQ intensity W_L3_P':'Pep Ctrl 3',
                                  'LFQ intensity SRT_L1_P':'Pep SRT 1','LFQ intensity SRT_L2_P':'Pep SRT 2',
                                  'LFQ intensity SRT_L3_P':'Pep SRT 3'})


df_common = df_common.drop(df_common[df_common['Pep Ctrl 1'] == 'CTRL'].index)
df_common = df_common.drop(df_common[df_common['Pep SRT 1'] == 'SRT'].index)


Ctrl = pd.DataFrame()
Ctrl['Pep Ctrl 1'] = df_common['Pep Ctrl 1'].astype(float)
Ctrl['Pep Ctrl 2'] = df_common['Pep Ctrl 2'].astype(float)
Ctrl['Pep Ctrl 3'] = df_common['Pep Ctrl 3'].astype(float)
Ctrl['Ctrl mean'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

SRT = pd.DataFrame()
SRT['Pep SRT 1'] = df_common['Pep SRT 1'].astype(float)
SRT['Pep SRT 2'] = df_common['Pep SRT 2'].astype(float)
SRT['Pep SRT 3'] = df_common['Pep SRT 3'].astype(float)
SRT['SRT mean'] = SRT.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)


SRTvCtrl = pd.DataFrame()
SRTvCtrl['SRT 1 vs Ctrl normalized'] = (df_common['Pep SRT 1'].astype(float)) / (df_common['SRT/Ctrl fold change'].astype(float))
SRTvCtrl['SRT 2 vs Ctrl normalized'] = (df_common['Pep SRT 2'].astype(float)) / (df_common['SRT/Ctrl fold change'].astype(float))
SRTvCtrl['SRT 3 vs Ctrl normalized'] = (df_common['Pep SRT 3'].astype(float)) / (df_common['SRT/Ctrl fold change'].astype(float))
SRTvCtrl['SRT vs Ctrl normalized average'] = SRTvCtrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

Ctrl = pd.DataFrame()
Ctrl['Pep Ctrl 1'] = df_common['Pep Ctrl 1'].astype(float)
Ctrl['Pep Ctrl 2'] = df_common['Pep Ctrl 2'].astype(float)
Ctrl['Pep Ctrl 3'] = df_common['Pep Ctrl 3'].astype(float)
Ctrl['Ctrl average'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

SRTvCtrl['Pep Ctrl 1'] = Ctrl['Pep Ctrl 1']
SRTvCtrl['Pep Ctrl 2'] = Ctrl['Pep Ctrl 2']
SRTvCtrl['Pep Ctrl 3'] = Ctrl['Pep Ctrl 3']
SRTvCtrl['Ctrl average'] = Ctrl['Ctrl average']

SRTvCtrl['Protein IDs'] = df_common['Protein ID']
SRTvCtrl['Protein Names'] = df_common['Protein names']
SRTvCtrl['Gene Names'] = df_common['Gene names']
SRTvCtrl['Fasta headers'] = df_common['Fasta header']
SRTvCtrl['Sequence'] = df_common['Sequence']

SRTvCtrl['SRT vs Ctrl fold change'] = SRTvCtrl['SRT vs Ctrl normalized average'] / Ctrl['Ctrl average']

SRTvCtrl.dropna(subset = ['SRT vs Ctrl fold change'], inplace=True)

MCIvCtrl_out = working_directory + '\\SRTvCtrl_nofilter.csv' 


with open(MCIvCtrl_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    SRTvCtrl.to_csv(filed,index=False)

print('Unfiltered data has been exported to working directory')
