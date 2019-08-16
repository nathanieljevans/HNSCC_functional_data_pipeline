'''


'''
from HNSCC_analysis_pipeline_lib import *
import sys
import os
from datetime import datetime
import pandas as pd

if __name__ == '__main__':

    _, output_dir = sys.argv
    adata = None
    for ppath in [x for x in os.listdir(output_dir) if '.' not in x]:
        try:
            print('joining %s ... \t' %ppath, end='')
            path = output_dir + '/' + ppath + '/HNSCC_processed_functional_data.csv'
            adata = pd.read_csv(path) if adata is None else adata.append(pd.read_csv(path), ignore_index=True, sort=False)
            print('complete.')
        except Exception as e:
            print('Failed to combine %s\n\t%s' %(ppath, str(e)))

    adata.to_csv(output_dir + '/HNSCC_all_functional_data.csv')
