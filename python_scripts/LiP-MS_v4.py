# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 15:58:24 2022

@author: lawashburn
"""

import pandas as pd
import numpy as np
import csv
from scipy import stats
import math
import scipy
from scipy.stats import ttest_ind

working_directory = input('Input path to working directory, where all results will be stored: ')
protein_path = input('Input path to Proteingroup.csv file: ')
peptides_path = input('Input path to Peptidegroup.csv file: ')
#protein_path = r"D:\LiP-MS_Haiyan\20220527\proteinGroups_abridged.csv"
#peptides_path = r"D:\LiP-MS_Haiyan\20220527\peptides_abridged.csv"
#working_directory = r"D:\LiP-MS_Haiyan\20221027"
p_cutoff = input('Enter p-value cutoff (0.05 is recommended): ')
up_thresh = input('Enter upregulated peptide threshold (recommended is 2): ')
down_thresh = input('Enter downregulated peptide threshold (recommended is -0.5): ')


#p_cutoff = 0.05
#up_thresh = 1.5
#down_thresh = 2/3


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
SRTvCtrl['SRT vs Ctrl fold change'] = SRTvCtrl['SRT vs Ctrl fold change'].interpolate(method='polynomial', order=2)
SRTvCtrl.dropna(subset = ['SRT vs Ctrl fold change'], inplace=True)



SRTvCtrl['p value'] = ttest_ind(SRTvCtrl[['SRT 1 vs Ctrl normalized', 'SRT 2 vs Ctrl normalized','SRT 3 vs Ctrl normalized']], SRTvCtrl[['Pep Ctrl 1', 'Pep Ctrl 2', 'Pep Ctrl 3']], axis=1)[1]

p_cutoff = float(p_cutoff)
up_thresh = float(up_thresh)
down_thresh = float(down_thresh)

SRTvCtrl_mask = SRTvCtrl['p value']<= p_cutoff
filtered_SRTvCtrl = SRTvCtrl[SRTvCtrl_mask]

SRTvCtrlUp_mask = filtered_SRTvCtrl['SRT vs Ctrl fold change']>=up_thresh
SRTvCtrlUp = filtered_SRTvCtrl[SRTvCtrlUp_mask]

SRTvCtrlDown_mask = filtered_SRTvCtrl['SRT vs Ctrl fold change']<=down_thresh
SRTvCtrlDown = filtered_SRTvCtrl[SRTvCtrlDown_mask]

SRTvCtrlUp_out = working_directory + '\\SRTvCtrl_Upregulated.csv'

SRTvCtrlDown_out = working_directory + '\\SRTvCtrl_downregulated.csv'

with open(SRTvCtrlUp_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    SRTvCtrlUp.to_csv(filed,index=False)

with open(SRTvCtrlDown_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    SRTvCtrlDown.to_csv(filed,index=False)
    
SRTvCtrlUpSeq = SRTvCtrlUp['Sequence'].values.tolist()
SRTvCtrlDownSeq = SRTvCtrlDown['Sequence'].values.tolist()
SRTvCtrlDownHit = []
SRTvsCtrlUpHit = []

for a in SRTvCtrlDownSeq:
    for b in SRTvCtrlUpSeq:
        if b in a:
            SRTvCtrlDownHit.append(a)
            SRTvsCtrlUpHit.append(b)
            
for a in SRTvCtrlUpSeq:
    for b in SRTvCtrlDownSeq:
        if b in a:
            SRTvCtrlDownHit.append(b)
            SRTvsCtrlUpHit.append(a)
           
SRTvCtrl_Hits = pd.DataFrame()
SRTvCtrl_Hits['Upregulated'] = SRTvsCtrlUpHit
SRTvCtrl_Hits['Downregulated'] = SRTvCtrlDownHit
SRTvCtrl_Hits['Sequence'] = SRTvCtrl_Hits['Upregulated']

SRTvCtrl_Hits_merge_store = pd.DataFrame()
if len(SRTvCtrlUp)>0:
    SRTvCtrl_Hits_merge = pd.merge(SRTvCtrl_Hits, SRTvCtrlUp, on=['Sequence'], how='inner')
    SRTvCtrl_Hits_merge_store = SRTvCtrl_Hits_merge_store.append(SRTvCtrl_Hits_merge)
else:
    pass

if len(SRTvCtrlDown)>0:
    SRTvCtrl_Hits_merge = pd.merge(SRTvCtrl_Hits, SRTvCtrlDown, on=['Sequence'], how='inner')
    SRTvCtrl_Hits_merge_store = SRTvCtrl_Hits_merge_store.append(SRTvCtrl_Hits_merge)
else:
    pass
    
SRTvCtrl_out = pd.DataFrame()
SRTvCtrl_out['Upregulated'] = SRTvCtrl_Hits_merge_store['Upregulated']
SRTvCtrl_out['Downregulated'] = SRTvCtrl_Hits_merge_store['Downregulated']
SRTvCtrl_out['Protein IDs'] = SRTvCtrl_Hits_merge_store['Protein IDs']
SRTvCtrl_out['Protein Names'] = SRTvCtrl_Hits_merge_store['Protein Names']
SRTvCtrl_out['Gene Names'] = SRTvCtrl_Hits_merge_store['Gene Names']
SRTvCtrl_out['Fasta headers'] = SRTvCtrl_Hits_merge_store['Fasta headers']

if len(SRTvCtrl_out)>0:
    hits_out = working_directory + '\\SRTvCtrl_hits.csv'
    with open(hits_out,'w',newline='') as filed:
        writerd = csv.writer(filed)
        SRTvCtrl_out.to_csv(filed,index=False)
else:
    print('No hits')