# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 11:05:28 2022

@author: lawashburn
"""

import pandas as pd
import numpy as np
import csv
import math
import scipy
from scipy.stats import ttest_ind

working_directory = input('Input path to working directory, where all results will be stored: ')
protein_path = input('Input path to Proteingroup.csv file: ')
peptides_path = input('Input path to Peptidegroup.csv file: ')
p_cutoff_protein = input('Enter p-value cutoff for proteins (0.05 is recommended): ')
up_thresh_protein = input('Enter upregulated protein threshold (recommended is 2): ')
down_thresh_protein = input('Enter downregulated protein threshold (recommended is -0.5): ')

p_cutoff_peptide = input('Enter p-value cutoff for peptide (0.05 is recommended): ')
up_thresh_peptide = input('Enter upregulated peptide threshold (recommended is 2): ')
down_thresh_peptide = input('Enter downregulated peptide threshold (recommended is -0.5): ')
p_cutoff_protein = float(p_cutoff_protein)
up_thresh_protein = float(up_thresh_protein)
down_thresh_protein = float(down_thresh_protein)
p_cutoff_peptide = float(p_cutoff_peptide)
up_thresh_peptide = float(up_thresh_peptide)
down_thresh_peptide = float(down_thresh_peptide)
# =============================================================================
# protein_path = r"D:\LiP-MS_Haiyan\20221108_EarlyFilter_Haiyan\Diagnostics\Input\proteinGroups_L - Copy.csv"
# peptides_path = r"D:\LiP-MS_Haiyan\20221108_EarlyFilter_Haiyan\Diagnostics\Input\peptides_L - Copy.csv"
# working_directory = r"D:\LiP-MS_Haiyan\20221108_EarlyFilter_Haiyan\Diagnostics\Input\test_out\test_again"
# p_cutoff_protein = 0.05
# up_thresh_protein = 2
# down_thresh_protein = -0.5
# p_cutoff_peptide = 0.05
# up_thresh_peptide = 2
# down_thresh_peptide = -0.5
# =============================================================================

protein = pd.read_csv(protein_path)
peptide = pd.read_csv(peptides_path)
protein = protein.rename(columns={'LFQ intensity SRT_L1_T':'SRT1_protein',
                                  'LFQ intensity SRT_L2_T':'SRT2_protein',
                                  'LFQ intensity SRT_L3_T':'SRT3_protein',
                                  'LFQ intensity W_L1_T':'Ctrl1_protein',
                                  'LFQ intensity W_L2_T':'Ctrl2_protein',
                                  'LFQ intensity W_L3_T':'Ctrl3_protein'})


peptide = peptide.rename(columns={'LFQ intensity SRT_L1_P':'SRT1_peptide',
                                  'LFQ intensity SRT_L2_P':'SRT2_peptide',
                                  'LFQ intensity SRT_L3_P':'SRT3_peptide',
                                  'LFQ intensity W_L1_P':'Ctrl1_peptide',
                                  'LFQ intensity W_L2_P':'Ctrl2_peptide',
                                  'LFQ intensity W_L3_P':'Ctrl3_peptide'})

SRT = pd.DataFrame()
SRT['SRT1_protein'] = protein['SRT1_protein']
SRT['SRT2_protein'] = protein['SRT2_protein']
SRT['SRT3_protein'] = protein['SRT3_protein']
SRT = SRT.iloc[1:,:]
SRT['mean_protein'] = SRT.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

Ctrl = pd.DataFrame()
Ctrl['Ctrl1_protein'] = protein['Ctrl1_protein']
Ctrl['Ctrl2_protein'] = protein['Ctrl2_protein']
Ctrl['Ctrl3_protein'] = protein['Ctrl3_protein']
Ctrl = Ctrl.iloc[1:,:]
Ctrl['mean_protein'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

protein_calc = pd.DataFrame()
protein_calc['Protein ID'] = protein['Protein IDs']
protein_calc['Protein names'] = protein['Protein names']
protein_calc['Gene name'] = protein['Gene names']
protein_calc['Fasta header'] = protein['Fasta headers']
protein_calc['SRT1_protein'] = SRT['SRT1_protein']
protein_calc['SRT2_protein'] = SRT['SRT2_protein']
protein_calc['SRT3_protein'] = SRT['SRT3_protein']
protein_calc['Ctrl1_protein'] = Ctrl['Ctrl1_protein']
protein_calc['Ctrl2_protein'] = Ctrl['Ctrl2_protein']
protein_calc['Ctrl3_protein'] = Ctrl['Ctrl3_protein']
protein_calc['SRT avg_protein'] = SRT['mean_protein']
protein_calc['Ctrl avg_protein'] = Ctrl['mean_protein']
protein_calc['SRT/Ctrl fold change_protein'] = protein_calc['SRT avg_protein'] / protein_calc['Ctrl avg_protein']

df_common = pd.merge(protein_calc, peptide, on=['Protein names'])
df_common['protein p-value'] = ttest_ind(df_common[['SRT1_protein', 'SRT2_protein','SRT3_protein']], 
                                df_common[['Ctrl1_protein', 'Ctrl2_protein', 'Ctrl3_protein']], axis=1, equal_var=True, alternative='two-sided')[1]
print(df_common)
df_common = df_common.drop(df_common[df_common['Ctrl1_peptide'] == 'CTRL'].index)
df_common = df_common.drop(df_common[df_common['SRT1_peptide'] == 'SRT'].index)
print(df_common)
#%%
Ctrl = pd.DataFrame()
Ctrl['Ctrl1_peptide'] = df_common['Ctrl1_peptide'].astype(float)
Ctrl['Ctrl2_peptide'] = df_common['Ctrl2_peptide'].astype(float)
Ctrl['Ctrl3_peptide'] = df_common['Ctrl3_peptide'].astype(float)
Ctrl['Ctrl mean_peptide'] = Ctrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)

SRTvCtrl_test = df_common["SRT/Ctrl fold change_protein"].astype(float)
SRTvCtrl_p = df_common['protein p-value'].astype(float)
SRTvCtrl_test = np.where(SRTvCtrl_test > p_cutoff_protein, df_common["SRT/Ctrl fold change_protein"],1)
SRTvCtrl_test = np.where(SRTvCtrl_test < (up_thresh_protein), df_common["SRT/Ctrl fold change_protein"],1)
SRTvCtrl_test = np.where(SRTvCtrl_p < down_thresh_protein, df_common["SRT/Ctrl fold change_protein"],1)

SRTvCtrl = pd.DataFrame()
SRTvCtrl['SRT 1 vs Ctrl normalized_peptide'] = (df_common['SRT1_peptide'].astype(float)) / (df_common['SRT/Ctrl fold change_protein'].astype(float))
SRTvCtrl['SRT 2 vs Ctrl normalized_peptide'] = (df_common['SRT2_peptide'].astype(float)) / (df_common['SRT/Ctrl fold change_protein'].astype(float))
SRTvCtrl['SRT 3 vs Ctrl normalized_peptide'] = (df_common['SRT3_peptide'].astype(float)) / (df_common['SRT/Ctrl fold change_protein'].astype(float))
SRTvCtrl['SRT vs Ctrl normalized average_peptide'] = SRTvCtrl.astype(float).replace(0, np.nan).mean(axis=1, skipna=True)
SRTvCtrl['Ctrl1_protein'] = df_common['Ctrl1_protein']
SRTvCtrl['Ctrl2_protein'] = df_common['Ctrl2_protein']
SRTvCtrl['Ctrl3_protein'] = df_common['Ctrl3_protein']
SRTvCtrl['Protein IDs'] = df_common['Protein ID']
SRTvCtrl['Protein Names'] = df_common['Protein names']
SRTvCtrl['Gene Names'] = df_common['Gene names']
SRTvCtrl['Fasta headers'] = df_common['Fasta header']
SRTvCtrl['Sequence'] = df_common['Sequence']
SRTvCtrl['SRT vs Ctrl fold change_peptide'] = SRTvCtrl['SRT vs Ctrl normalized average_peptide'] / Ctrl['Ctrl mean_peptide']
SRTvCtrl['peptide p-value'] = ttest_ind(df_common[['SRT1_peptide', 'SRT2_peptide','SRT3_peptide']], 
                                df_common[['Ctrl1_peptide', 'Ctrl2_peptide', 'Ctrl3_peptide']], axis=1, equal_var=True, alternative='two-sided')[1]

SRTvCtrlunfiltered_out = working_directory + '\\SRTvCtrl_unfiltered_normalized_peptides.csv'

with open(SRTvCtrlunfiltered_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    SRTvCtrl.to_csv(filed,index=False)


SRTvCtrl_mask = SRTvCtrl['peptide p-value']<= p_cutoff_peptide
filtered_SRTvCtrl= SRTvCtrl[SRTvCtrl_mask]


SRTvCtrlUp_mask = filtered_SRTvCtrl['SRT vs Ctrl fold change_peptide']>=up_thresh_peptide
SRTvCtrlUp = filtered_SRTvCtrl[SRTvCtrlUp_mask]

SRTvCtrlDown_mask = filtered_SRTvCtrl['SRT vs Ctrl fold change_peptide']<=down_thresh_peptide
SRTvCtrlDown = filtered_SRTvCtrl[SRTvCtrlDown_mask]

SRTvCtrlUp_out = working_directory + '\\SRTvCtrl_upregulated.csv'
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
SRTvCtrl_Hits_merge = pd.merge(SRTvCtrl_Hits, SRTvCtrlUp, on=['Sequence'], how='inner')

SRTvCtrl_out = pd.DataFrame()
SRTvCtrl_out['Upregulated'] = SRTvCtrl_Hits_merge['Upregulated']
SRTvCtrl_out['Downregulated'] = SRTvCtrl_Hits_merge['Downregulated']
SRTvCtrl_out['Protein IDs'] = SRTvCtrl_Hits_merge['Protein IDs']
SRTvCtrl_out['Protein Names'] = SRTvCtrl_Hits_merge['Protein Names']
SRTvCtrl_out['Gene Names'] = SRTvCtrl_Hits_merge['Gene Names']
SRTvCtrl_out['Fasta headers'] = SRTvCtrl_Hits_merge['Fasta headers']

SRTvCtrl_hits_out = working_directory + '\\SRTvCtrl_hits.csv'
with open(SRTvCtrl_hits_out,'w',newline='') as filed:
    writerd = csv.writer(filed)
    SRTvCtrl_out.to_csv(filed,index=False)
