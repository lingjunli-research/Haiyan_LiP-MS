# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 10:52:07 2022

@author: lawashburn
"""

import pandas as pd
import csv


working_directory = input('Input path to working directory, where all results will be stored: ')
#working_directory = r"C:\Users\lawashburn\Documents\LiP-MS_Haiyan\20220527\test_out"
SRTvCtrl_path = input('Enter path to SRTvCtrl .csv file updated with p values: ')
#SRTvCtrl_path = r"C:\Users\lawashburn\Documents\LiP-MS_Haiyan\20220527\test_out\SRTvCtrl_nofilter.csv"
p_cutoff = input('Enter p-value cutoff (0.05 is recommended): ')
up_thresh = input('Enter upregulated peptide threshold (recommended is 2): ')
down_thresh = input('Enter downregulated peptide threshold (recommended is -0.5): ')


#p_cutoff = 0.05
#up_thresh = 1.5
#down_thresh = 2/3

p_cutoff = float(p_cutoff)
up_thresh = float(up_thresh)
down_thresh = float(down_thresh)

SRTvCtrl = pd.read_csv(SRTvCtrl_path)

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
    MCIvCtrl_hits_out = working_directory + '\\SRTvCtrl_hits.csv'
    with open(SRTvCtrl_out,'w',newline='') as filed:
        writerd = csv.writer(filed)
else:
    print('No results')