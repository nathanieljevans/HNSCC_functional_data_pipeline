'''
@author: Nate Evans
@date: 7/16/2019
@title: HNSCC Functional Data Analysis Pipeline

The purpose of this script is to provide the resources to enable an automated
analysis pipeline of functional drug response data in the HNSCC project at
OHSU led by Shannon McWeeney and Molly Kulesz-Martin.

'''

import pandas as pd
import numpy as np
pd.options.display.width = 0

def parse_pathname(path, verbose=False):
    '''
    parse relevant information from data path, must return
        - lab_id [unique 5 digit patient identifier]
        - norm [Blank490, 490]
        - version_id <str>

    file name must be in form:
        lab_id=XXXXX-norm=XXXXXX-plate_version_id=XXXXX.xlsx

    [TODO] exceptions should be thrown if file name is atypical
    '''

    name = path.split('/')[-1][:-5].split('-')

    labels = [n.split('=')[0] for n in name]
    values = [n.split('=')[-1] for n in name]

    if verbose: print('lables: %r' %labels)
    if verbose: print('values: %r' %values)

    return values



def get_plate_data(data_path, verbose=True):
    '''
    each plate has 16 rows of drug data, then a empty row, then the next plate. Varying number of plates

    inputs
        data_path <str> path to single patient plate data

    outputs
        dataframe
    '''

    lab_id, norm, version_id = parse_pathname(data_path)

    if verbose: print('---------------------------------------------------------')
    if verbose: print( 'please double check file name parsing is accurate')
    if verbose: print( 'file name: %s' %data_path.split('/')[-1])
    if verbose: print('lab_id: %s' %lab_id)
    if verbose: print('norm: %s' %norm)
    if verbose: print('version_id: %s' %version_id)
    if verbose: print('---------------------------------------------------------')

    allplates = pd.read_excel(data_path, header=None)

    nplates = (allplates.shape[0] + 1) / 18

    assert nplates%1.0 == 0, 'atypical number of rows, check data format'

    if verbose: print('assay has %d plates' %int(nplates))

    plates = []
    i = 1 # skip original header
    for p in range(int(nplates)):
        dat = pd.DataFrame( allplates.values[i:(i+16),:])
        dat = dat.assign(plate_row = dat[0]).assign(norm_type = dat[25], plate_num = p+1, lab_id = lab_id, assay_version_id=version_id).drop(labels = [0,25], axis='columns')
        dat = pd.melt(dat, id_vars=['lab_id', 'norm_type', 'plate_num', 'plate_row','assay_version_id'], value_vars=None, var_name='plate_col', value_name='optical_density', col_level=None)

        plates.append( dat )
        i += 16 + 2 # skip empty row + header

    # Frame = Frame.append(pandas.DataFrame(data = SomeNewLineOfData), ignore_index=True)
    df = plates[0]
    for p in plates[1:]:
        df = df.append(p, ignore_index=True)

    return df

def get_plate_map(map_path, verbose=True):
    '''
    get the plate mapping data from the excel plate map

    input
        map_path <str> local or absolute path to the plate map

    output
        <dataframe> plate mapping data. header = [plate_number, row, col, drug, conc, version_id]
    '''
    meta = pd.read_excel(map_path, header=None, sheet_name='meta')

    plate_version = meta[meta[0]=='version_id'].values[0][1]
    num_sheets = int( meta[meta[0]=='num_plates'][1] )

    if verbose: print('plate version id is: %s' %plate_version)
    if verbose: print('There are %d plates in this plate map' %num_sheets)

    map_data = pd.read_excel(map_path, header=0, sheet_name=list(range(1,((num_sheets*2)+1))))

    plates = []
    for p in range(num_sheets):
        plate_num = p+1
        concs = pd.melt( map_data[p*2+1], id_vars=['row'], var_name='col', value_name='conc')
        drugs = pd.melt( map_data[p*2+2], id_vars=['row'], var_name='col', value_name='inhibitor')
        plate_map = concs.merge(drugs, on=['row','col'])
        plate_map = plate_map.assign(map_version_id = plate_version, plate_number=plate_num)

        plates.append(plate_map)

    dat = plates[0]
    for p in plates[1:]:
        dat = dat.append(p, ignore_index=True)

    return dat

if __name__ == '__main__':

    # diagnostics
    plate_path = './data/HNSCC_panels/lab_id=10348-norm=Blank490-plate_version_id=OHSU_HNSCC_derm002.xlsx'

    plate_dat = get_plate_data(plate_path)

    #print( plate_dat )

    map_path = './HNSCC_vs1-COMBINATION_panel_map-nje.xlsx'

    plate_map = get_plate_map(map_path)

    print( plate_dat.head() )

    print( plate_map.head() )

    plate_map.to_csv('./for_nate_test.csv')

    dr_dat = plate_dat.merge(plate_map, how='left', left_on=['plate_row','plate_col','plate_num'], right_on=['row','col','plate_number'])

    dr_dat = dr_dat.drop(['row','col','plate_number'], axis='columns')

    print( dr_dat.head() )

    print(plate_map.groupby(['inhibitor']).size() )

    print(dr_dat.groupby(['inhibitor']).size())
