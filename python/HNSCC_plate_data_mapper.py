'''
@author: Nate Evans
@date: 7/16/2019
@title: HNSCC Plate Mapper

This script is intended to be run from command line by:
    $ python HNSCC_plate_data_mapper.py ./path/to/plate_maps/dir ./path/to/photospec/output/distribution

./path/to/plate_maps/dir should contain all [and do not include other file types] relevant plate maps saved with naming convention:
    HNSCC_plate_map-version_id=XXX.xlsx

./path/to/photospec/output/distribution should contain [only] files with naming convention:
    lab_id=XXXXX-norm=XXXX-plate_version_id=XXX.xlsx
        where XXX signifys varied input

This script will save a single csv file in the folder (it will create one if it does not exist)
    ./output/HNSCC_FUNCTIONAL_ASSAY_MAPPED_DOSE_RESPONSE_DATA
with the naming convention
    HNSCC_funcData_MM_DD_YYYY.csv

'''

from HNSCC_analysis_pipeline_lib import *
import sys
import os
from datetime import datetime

_, map_dir, data_dir = sys.argv

# create dictionary, version_id -> map_data
print('loading plate maps from: %s ... ' %map_dir)
map_data = {}
for v in os.listdir(map_dir):
    version_id = v[:-5].split('=')[-1]
    print(v)
    map_data[version_id] = get_plate_map(map_dir + '/' + v)
print('\n')

print('mapping panel data...')
# map the data
panels = []
failures = []
for p in os.listdir(data_dir):

    try:
        pdat = get_plate_data(data_dir + '/' + p)
        plate_map = map_data[pdat['assay_version_id'].values[0]]
        panels.append( pdat.merge(plate_map, how='left', left_on=['plate_row','plate_col','plate_num'], right_on=['row','col','plate_number']).drop(['row','col','plate_number'], axis='columns') )
    except:
        failures.append(p)
        raise
        #rint('Failure, file name: %s' %p)

print('%d failure(s) [%r]' %(len(failures), failures))
# combine all data
data = panels[0]
for p in panels[1:]:
    data = data.append(p, ignore_index=True)

if not os.path.exists('./output/HNSCC_FUNCTIONAL_ASSAY_MAPPED_DOSE_RESPONSE_DATA'):
    os.makedirs('./output/HNSCC_FUNCTIONAL_ASSAY_MAPPED_DOSE_RESPONSE_DATA')

data.to_csv('./output/HNSCC_FUNCTIONAL_ASSAY_MAPPED_DOSE_RESPONSE_DATA/HNSCC_funcData_%s.csv' %datetime.today().strftime('%m_%d_%Y'))

#print( data.groupby(['inhibitor']).size() )
