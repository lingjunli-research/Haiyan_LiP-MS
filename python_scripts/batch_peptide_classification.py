# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 11:05:37 2022

@author: lawashburn
"""

import pandas as pd
import csv
import numpy as np
import time
start = time.time()


S_file_path = r"D:\LiP-MS_Haiyan\20221206\Sample Data\Lung_significant peptides.csv"
P_file_path = r"D:\LiP-MS_Haiyan\20221206\Sample Data\Lung_peptides_in P.csv"
T_file_path = r"D:\LiP-MS_Haiyan\20221206\Sample Data\Lung_peptides_in T.csv"
Output_path = r"D:\LiP-MS_Haiyan\20221206\Sample Data"

S_file = pd.read_csv(S_file_path)
T_file = pd.read_csv(T_file_path)
P_file = pd.read_csv(P_file_path)

S_file.columns = [str(col) + '_S' for col in S_file.columns]
P_file.columns = [str(col) + '_P' for col in P_file.columns]
T_file.columns = [str(col) + '_T' for col in T_file.columns]

S_file['Sequence_S'] = S_file['Sequence_S'].fillna('NA')
P_file['Sequence_P'] = P_file['Sequence_P'].fillna('NA')
T_file['Sequence_T'] = T_file['Sequence_T'].fillna('NA')

S_peptides = S_file['Sequence_S'].values.tolist()
P_peptides = P_file['Sequence_P'].values.tolist()
T_peptides = T_file['Sequence_T'].values.tolist()

S_matches = []
P_matches = []
T_matches = []

P_matches_minimize = []
T_matches_minimize = []

for f in P_peptides:
    if f not in S_peptides:
        P_matches_minimize.append(f)

for g in T_peptides:
    if g not in S_peptides:
        T_matches_minimize.append(g)


#rule 2
for a in S_peptides:
    for b in T_matches_minimize:
        if a in b:
            S_matches.append(a)
            T_matches.append(b)
            if b in P_matches_minimize:
                P_matches.append(b)
            else:
                P_matches.append(np.nan)



results_df = pd.DataFrame()
results_df['Sequence S'] = S_matches
results_df['Sequence P'] = P_matches
results_df['Sequence T'] = T_matches

results_df = results_df.merge(S_file, left_on='Sequence S',right_on='Sequence_S',how='left')
results_df = results_df.merge(P_file, left_on='Sequence P',right_on='Sequence_P',how='left')
results_df = results_df.merge(T_file, left_on='Sequence T',right_on='Sequence_T',how='left')


file_path = Output_path + '\\' + 'batch_test.csv'
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        results_df.to_csv(filec,index=False) 
end = time.time()
elapsed_time = end-start
print('Elapsed time: ',elapsed_time)