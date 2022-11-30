# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 09:08:40 2022

@author: lawashburn
"""

import pandas as pd
import csv
import math
import scipy
from scipy.stats import ttest_ind
import numpy as np

protein_path = r"D:\LiP-MS_Haiyan\20221130\proteinGroups_11302022.csv"
peptide_path = r"D:\LiP-MS_Haiyan\20221130\peptides_11302022.csv"
output_path = r"D:\LiP-MS_Haiyan\20221130"
p_cutoff = 0.05
min_FC = 2
min_log_FC = 1

#Determine fold change based on protein level results
protein_report = pd.read_csv(protein_path)
protein_report = protein_report.replace(0, np.NaN) #replace empty values with NaN so mean can be taken without impact from 0s
protein_report['SRT mean protein'] = protein_report[['LFQ intensity SRT_L1_T', 'LFQ intensity SRT_L2_T','LFQ intensity SRT_L3_T']].mean(axis=1,
                                                                                                                                        skipna=True) #the mean is taken without considering 0s
protein_report['Ctrl mean protein'] = protein_report[['LFQ intensity WT_L1_T', 'LFQ intensity WT_L2_T','LFQ intensity WT_L3_T']].mean(axis=1,
                                                                                                                                        skipna=True) #the mean is taken without considering 0s
protein_report['protein fold change'] = protein_report['SRT mean protein']/protein_report['Ctrl mean protein']
protein_report['Log2(protein fold change)'] = abs(np.log2(protein_report['protein fold change']))

#Determine p-value between ctrl and exp peptides
peptide_report = pd.read_csv(peptide_path)
merged_pep_prot_report = peptide_report.merge(protein_report,on=['Protein names','Gene names'],how='left')

print('Number of protein entries: ',len(protein_report))
print('Number of peptide entries: ',len(peptide_report))
print('Number of entries post-protein/peptide merge: ',len(merged_pep_prot_report))

merged_pep_prot_report['protein p-value'] = ttest_ind(merged_pep_prot_report[['LFQ intensity SRT_L1_T', 'LFQ intensity SRT_L2_T','LFQ intensity SRT_L3_T']], 
                                merged_pep_prot_report[['LFQ intensity WT_L1_T', 'LFQ intensity WT_L2_T','LFQ intensity WT_L3_T']], axis=1, equal_var=True, alternative='two-sided',nan_policy='omit')[1]


merged_pep_prot_report['unaltered fold change (SetA)'] = merged_pep_prot_report['protein fold change'] #in this column the FC will not be changed regardless of p-value
merged_pep_prot_report['p-value filtered protein fold change (SetB)'] = merged_pep_prot_report['protein fold change'] #in this column the FC will be equal to 1 if p-value is > 0.05
merged_pep_prot_report['p-value OR log2(FC) filtered protein fold change (SetC)'] = merged_pep_prot_report['protein fold change'] #in this column the FC will be equal to 1 if p-value is > 0.05 OR the FC is insignificant
merged_pep_prot_report['p-value AND log2(FC) filtered protein fold change (SetD)'] = merged_pep_prot_report['protein fold change'] #in this column the FC will be equal to 1 if p-value is > 0.05 AND the FC is insignificant

merged_pep_prot_report.loc[merged_pep_prot_report['protein p-value'] > p_cutoff, ['p-value filtered protein fold change (SetB)']] = 1
merged_pep_prot_report.loc[merged_pep_prot_report['protein p-value'].isnull(), ['p-value filtered protein fold change (SetB)']] = 1

merged_pep_prot_report.loc[merged_pep_prot_report['protein p-value'] > p_cutoff, ['p-value OR log2(FC) filtered protein fold change (SetC)']] = 1
merged_pep_prot_report.loc[merged_pep_prot_report['protein p-value'].isnull(), ['p-value OR log2(FC) filtered protein fold change (SetC)']] = 1
merged_pep_prot_report.loc[abs(merged_pep_prot_report['Log2(protein fold change)']) < min_log_FC, ['p-value OR log2(FC) filtered protein fold change (SetC)']] = 1


merged_pep_prot_report.loc[abs((merged_pep_prot_report['Log2(protein fold change)'])<min_log_FC) & (merged_pep_prot_report['protein p-value']> p_cutoff),
                    ['p-value AND log2(FC) filtered protein fold change (SetD)']] = 1
merged_pep_prot_report.loc[abs((merged_pep_prot_report['Log2(protein fold change)'])<min_log_FC) & (merged_pep_prot_report['protein p-value'].isnull()),
                    ['p-value AND log2(FC) filtered protein fold change (SetD)']] = 1


merged_pep_prot_report = merged_pep_prot_report.replace(0, np.NaN) #replace empty values with NaN so mean can be taken without impact from 0s
merged_pep_prot_report['Ctrl mean peptide'] = merged_pep_prot_report[['LFQ intensity WT_L1_P', 'LFQ intensity WT_L2_P','LFQ intensity WT_L3_P']].mean(axis=1,
                                                                                                                                        skipna=True) #the mean is taken without considering 0s
merged_pep_prot_report = merged_pep_prot_report.replace(0, np.NaN) #replace empty values with NaN so mean can be taken without impact from 0s
merged_pep_prot_report['SetA (LFQ intensity SRT_L1_P normalized)'] = (merged_pep_prot_report['LFQ intensity SRT_L1_P'] / 
                                                                                   merged_pep_prot_report['unaltered fold change (SetA)'])
merged_pep_prot_report['SetA (LFQ intensity SRT_L2_P normalized)'] = (merged_pep_prot_report['LFQ intensity SRT_L2_P'] / 
                                                                                   merged_pep_prot_report['unaltered fold change (SetA)'])
merged_pep_prot_report['SetA (LFQ intensity SRT_L3_P normalized)'] = (merged_pep_prot_report['LFQ intensity SRT_L3_P'] / 
                                                                                   merged_pep_prot_report['unaltered fold change (SetA)'])


merged_pep_prot_report['SetB (LFQ intensity SRT_L1_P normalized)'] = (merged_pep_prot_report['LFQ intensity SRT_L1_P'] / 
                                                                                   merged_pep_prot_report['p-value filtered protein fold change (SetB)'])
merged_pep_prot_report['SetB (LFQ intensity SRT_L2_P normalized)'] = (merged_pep_prot_report['LFQ intensity SRT_L2_P'] / 
                                                                                   merged_pep_prot_report['p-value filtered protein fold change (SetB)'])
merged_pep_prot_report['SetB (LFQ intensity SRT_L3_P normalized)'] = (merged_pep_prot_report['LFQ intensity SRT_L3_P'] / 
                                                                                   merged_pep_prot_report['p-value filtered protein fold change (SetB)'])

merged_pep_prot_report['SetC (LFQ intensity SRT_L1_P normalized)'] = (merged_pep_prot_report['LFQ intensity SRT_L1_P'] / 
                                                                                   merged_pep_prot_report['p-value OR log2(FC) filtered protein fold change (SetC)'])
merged_pep_prot_report['SetC (LFQ intensity SRT_L2_P normalized)'] = (merged_pep_prot_report['LFQ intensity SRT_L2_P'] / 
                                                                                   merged_pep_prot_report['p-value OR log2(FC) filtered protein fold change (SetC)'])
merged_pep_prot_report['SetC (LFQ intensity SRT_L3_P normalized)'] = (merged_pep_prot_report['LFQ intensity SRT_L3_P'] / 
                                                                                   merged_pep_prot_report['p-value OR log2(FC) filtered protein fold change (SetC)'])

merged_pep_prot_report['SetD (LFQ intensity SRT_L1_P normalized)'] = (merged_pep_prot_report['LFQ intensity SRT_L1_P'] / 
                                                                                   merged_pep_prot_report['p-value AND log2(FC) filtered protein fold change (SetD)'])
merged_pep_prot_report['SetD (LFQ intensity SRT_L2_P normalized)'] = (merged_pep_prot_report['LFQ intensity SRT_L2_P'] / 
                                                                                   merged_pep_prot_report['p-value AND log2(FC) filtered protein fold change (SetD)'])
merged_pep_prot_report['SetD (LFQ intensity SRT_L3_P normalized)'] = (merged_pep_prot_report['LFQ intensity SRT_L3_P'] / 
                                                                                   merged_pep_prot_report['p-value AND log2(FC) filtered protein fold change (SetD)'])
merged_pep_prot_report = merged_pep_prot_report.replace(0, np.NaN) #replace empty values with NaN so mean can be taken without impact from 0s
merged_pep_prot_report['SetA SRT mean intensity peptide'] = merged_pep_prot_report[['SetA (LFQ intensity SRT_L1_P normalized)', 
                                                                                                            'SetA (LFQ intensity SRT_L2_P normalized)',
                                                                                                            'SetA (LFQ intensity SRT_L3_P normalized)']].mean(axis=1,
                                                                                                                                        skipna=True) #the mean is taken without considering 0s

merged_pep_prot_report['SetB SRT mean intensity peptide'] = merged_pep_prot_report[['SetB (LFQ intensity SRT_L1_P normalized)', 
                                                                                                            'SetB (LFQ intensity SRT_L2_P normalized)',
                                                                                                            'SetB (LFQ intensity SRT_L3_P normalized)']].mean(axis=1,
                                                                                                                                        skipna=True) #the mean is taken without considering 0s                                                                                                                                                              
    
merged_pep_prot_report['SetC SRT mean intensity peptide'] = merged_pep_prot_report[['SetC (LFQ intensity SRT_L1_P normalized)', 
                                                                                                            'SetC (LFQ intensity SRT_L2_P normalized)',
                                                                                                            'SetC (LFQ intensity SRT_L3_P normalized)']].mean(axis=1,
                                                                                                                                        skipna=True) #the mean is taken without considering 0s    
  
merged_pep_prot_report['SetD SRT mean intensity peptide'] = merged_pep_prot_report[['SetD (LFQ intensity SRT_L1_P normalized)', 
                                                                                                            'SetD (LFQ intensity SRT_L2_P normalized)',
                                                                                                            'SetD (LFQ intensity SRT_L3_P normalized)']].mean(axis=1,
                                                                                                                                        skipna=True) #the mean is taken without considering 0s                                                                                                                                                                 
    
merged_pep_prot_report['SetA SRT mean fold change peptide normalized'] = merged_pep_prot_report['SetA SRT mean intensity peptide']/merged_pep_prot_report['Ctrl mean peptide']                                                                                                                                                             
merged_pep_prot_report['SetB SRT mean fold change peptide normalized'] = merged_pep_prot_report['SetB SRT mean intensity peptide']/merged_pep_prot_report['Ctrl mean peptide'] 
merged_pep_prot_report['SetC SRT mean fold change peptide normalized'] = merged_pep_prot_report['SetC SRT mean intensity peptide']/merged_pep_prot_report['Ctrl mean peptide'] 
merged_pep_prot_report['SetD SRT mean fold change peptide normalized'] = merged_pep_prot_report['SetD SRT mean intensity peptide']/merged_pep_prot_report['Ctrl mean peptide'] 

                                                                                                                                                              
file_out_path = output_path + '\\Updated_report_20221130_3.csv'
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        merged_pep_prot_report.to_csv(filec,index=False)