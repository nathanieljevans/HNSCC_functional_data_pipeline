# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 13:38:01 2019

@author: nathaniel evans

This script combines the whole pipeline into one easy-to-run script
"""
from HNSCC_analysis_pipeline_lib import *
import sys
import os
from datetime import datetime as dt
import pandas as pd
import numpy as np
import shutil
from keras.models import load_model
import argparse

##############################################################################
#                              globals                                       #
##############################################################################

REPROCESS = True             # careful with this, all previously processed data will be erased including assay ids.  
DO_RAISE = False          # raise exceptions in sequential processing job. 

PLATEMAP_DIR = '../plate_maps/'
RAW_DATA_DIR = '../data/'

OUTPUT_DIR = '../output'
ASSAYID_PATH = '../.assay_ids'

SENSITIVE_QUANTILE = 0.2
RESISTSNT_QUANTILE = 0.8

PROCESSED_DATA_PATH = '../output/HNSCC_all_functional_data.csv'

ATYP_CLASSIFIER_MODEL_PATH = '../../atypical_doseresponse_classifier/python/classifier/best_model.20-0.10.h5'

##############################################################################
#                         house keeping                                      #
##############################################################################

print('-------------------------------------------------------------')
      
print('Non-critical exceptions during processing will be raised: %s' %str(DO_RAISE))
print(f'removing output directory & reprocessing: {REPROCESS}')

if REPROCESS and os.path.isdir(OUTPUT_DIR): 
    shutil.rmtree(OUTPUT_DIR, ignore_errors=True)
    open(ASSAYID_PATH, 'w').close()

tic = dt.now()
_, map_dir, data_dir = None, PLATEMAP_DIR, RAW_DATA_DIR

dirlist = os.listdir(data_dir)

if not REPROCESS: 
    started = [x for x in os.listdir(OUTPUT_DIR) if os.path.isdir(OUTPUT_DIR + '/' + x)]

    already_processed = []
    for x in started: 
        if "HNSCC_processed_functional_data.csv" in os.listdir(OUTPUT_DIR + '/' + x):
            already_processed.append(x)
        else: 
            print(f'{x} - not completed, removing directory.')
            shutil.rmtree(OUTPUT_DIR + '/' + x, ignore_errors=True)

    print(f'number of files already processed: {len(already_processed)}')
    to_process = set(dirlist) - set(already_processed)
    print(f'number of files to process: {len(to_process)}')
else: 
    already_processed = []
    to_process = dirlist

print('-------------------------------------------------------------')

failures = []
errors = []

for i,p in enumerate(to_process):
    try:
        print('Panels processed: [%d/%d]' %(i + len(already_processed), len(dirlist)))
        process(data_dir + '/' + p, platemap_dir=map_dir, verbose=False, do_raise=DO_RAISE)
    except Exception as e:
        print('Processing failed. See error log for details.')
        failures.append(p)
        errors.append(e)
        if DO_RAISE: raise

print('There were %d failures. \n\t %r' %(len(failures), failures))
print()
print('Error Messages: ')
[print(str(e)) for e in errors]

print(f'Processing complete. Time elapsed: {dt.now() - tic}')
print('###########################################')
print('###########################################')

##############################################################################
#                         combine_processed_panels.py                        #
##############################################################################
print('Combining all processed data in the output directory...')

_, output_dir = None, OUTPUT_DIR
adata = None
for ppath in os.listdir(output_dir):
    try:
        print('joining %s ... \t' %ppath, end='')
        path = output_dir + '/' + ppath + '/HNSCC_processed_functional_data.csv'
        adata = pd.read_csv(path) if adata is None else adata.append(pd.read_csv(path), ignore_index=True, sort=False)
        print('complete.')
    except Exception as e:
        print('Failed to combine %s\n\t%s' %(ppath, str(e)))

adata = adata.drop('Unnamed: 0', axis='columns')
adata.to_csv(output_dir + '/HNSCC_all_functional_data.csv')
print('data combination complete.')
print('###########################################')
print('###########################################')

##############################################################################
#                         patient_sensitivity_assignment.py                  #
##############################################################################
print('Assigning patient sensitivity labels...')

_, data_path = None, PROCESSED_DATA_PATH

data = pd.read_csv(data_path, low_memory=False)
in_shp = data.shape

data.lab_id = data.lab_id.astype(str)
data.auc = data.auc.astype(float)

# in case already has call
if 'call' in data.columns:
    data = data.drop(['call'], axis='columns')

#[print(f' {x}') for x in data['auc']]
data2 = data[~pd.isnull(data['auc'])]

assignments = pd.DataFrame(columns=['inhibitor','auc','lab_id','call'])
for inhib in data2.inhibitor.unique():

    inhib_dat = data2[data2['inhibitor'] == inhib]

    inhib_aucs = inhib_dat.auc #.unique()

    s,r = np.quantile(inhib_aucs, [SENSITIVE_QUANTILE, RESISTSNT_QUANTILE])
    #print(f'quantiles for {inhib}: {s}, {r}')

    #print(inhib_dat[['lab_id','inhibitor','auc']].drop_duplicates().head())

    res = inhib_dat[['lab_id','inhibitor','auc']].drop_duplicates().assign(sens=lambda x: x.auc < s, res=lambda x: x.auc > r)
    res = res.assign(call = ['sens' if s else 'res' if r else 'int' for s,r in zip(res['sens'], res['res'])]).drop(['sens', 'res'], axis='columns')

    #print(res.head())

    assignments = assignments.append(res, ignore_index=True, sort=False)

#print(assignments.head())

data = data.merge(assignments, how='left', on=['lab_id', 'inhibitor', 'auc'])

#print(data[['lab_id', 'inhibitor', 'auc', 'call']].head(20))

data.to_csv(data_path)
out_shp = data.shape

print(f'input shape: {in_shp}\noutput shape: {out_shp}')
print('Sensitivity assignment complete.')
print('###########################################')
print('###########################################')

##############################################################################
#                         python predict_hermetic_transitions.py             #
##############################################################################
      
print('starting atypical predictions...')

# data in
_, data_path, model_path = None, PROCESSED_DATA_PATH, ATYP_CLASSIFIER_MODEL_PATH

data = pd.read_csv(data_path, low_memory=False)
in_shp = data.shape

data.lab_id = data.lab_id.astype(str)
data.auc = data.auc.astype(float)

herm_df = data[['lab_id', 'inhibitor', 'conc_norm', 'cell_viab', 'plate_num', 'panel_id']]

print(herm_df.head())

# filter out all controls
print(f'data length before controls filter: {herm_df.shape[0]}')
herm_df = herm_df[~herm_df.inhibitor.isin(['DMSO', 'F-S-V', 'NONE'])]
print(f'data length after controls filter: {herm_df.shape[0]}')

model = load_model(model_path)

n = 0
i = 0
res = {x:[] for x in ['lab_id', 'inhibitor', 'plate_num', 'panel_id', 'conc_norm', 'atyp_prob']}
for inhib in herm_df.inhibitor.unique():
    #print(f'processing inhibitor: {inhib}')
    inhib_df = herm_df[herm_df.inhibitor == inhib]
    for pat in inhib_df.lab_id.unique():
        #print(f'processing patient: {pat}')
        pat_df = inhib_df[inhib_df.lab_id == pat]
        for pid in pat_df.panel_id.unique():
            #print(f'processing panel id: {pid}')
            df_pid = pat_df[pat_df.panel_id == pid]
            for pnum in df_pid.plate_num.unique():
                #print(f'processing panel number: {pnum}')
                df = df_pid[df_pid.plate_num == pnum]
                #print(df)
                assert df.shape[0] == 7, f'wrong number of doses, should only be 7, got {df.shape[0]} [inhibitor: {inhib} | lab_id: {pat} | panel_id: {pid}] \n {df}'
                df = df.sort_values(by=['conc_norm'], axis=0)
                conc = df['conc_norm'].values
                viab = df['cell_viab'].values

                X = np.array([conc, viab]).reshape(-1, 2, 7, 1)
                #print(X.shape)
                #print(X)
                Y = model.predict(X).reshape(2,7)

                #print(Y)
                n += 1
                #print(Y[0,:])
                #print(Y[1,:])
                #print(conc)
                #x = 3/0

                for typ, atyp, dose in zip(Y[0,:], Y[1,:], conc):
                    res['lab_id'].append(pat)
                    res['inhibitor'].append(inhib)
                    res['plate_num'].append(pnum)
                    res['panel_id'].append(pid)
                    res['conc_norm'].append(dose)
                    res['atyp_prob'].append(atyp / (atyp + typ))
                    i += 1
print()
print('processing complete.')

res = pd.DataFrame(res)

print(f'processing counter (i): {i}')
print(f'processing counter (n): {n}')

print(f'Results length: {res.shape[0]} [input data length: {data.shape[0]}]')

data = data.merge(res, how='left', on=['lab_id', 'inhibitor', 'plate_num', 'panel_id', 'conc_norm'])
out_shp = data.shape

print(f'input shape: {in_shp}\noutput shape: {out_shp}')

data.to_csv(data_path)

print('atypical predictions complete.')
print('###########################################')
print('###########################################')

##############################################################################
#                         add annotations                                    #
##############################################################################



##############################################################################
#                         final data cleaning                                #
##############################################################################