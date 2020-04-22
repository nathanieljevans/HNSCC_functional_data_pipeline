'''

''' 

import sys
sys.path.append('../../python/')
import HNSCC_analysis_pipeline_lib as lib 


lab_id          =  'tst01'
norm            =  'Blank 490'
plate_num_range =  (1,1)
plate_row_range =  ('A', 'P')
plate_col_range =  (1, 24)
platemap_id     =  'test001_platemap'


if __name__ == '__main__': 

    platemap_path = 'test001_platemap.xlsx'
    plate_data_path = 'lab_id=tst01-norm=blank490-plate_version_id=test001_platemap-note=test001.xlsx'

    p = lib.panel(plate_path=plate_data_path, platemap_dir = './', verbose=True)

    p.map_data()

    print()
    print('-------------------------------------------------------------------------')
    print('processing:', p.plate_path)
    print('total number of rows in data after mapping:', p.data.shape[0])
    print('number of unique inhibitors:', p.data.inhibitor.unique().shape[0])
    print('platemap being used:', p.platemap_path)
    print('-------------------------------------------------------------------------')
    print()

    print(p.data[['inhibitor', 'conc', 'optical_density']])

    print('testing mapping...')

    for i, row in p.data[['inhibitor', 'conc', 'optical_density']].iterrows(): 
        inh = int(row.inhibitor.split('-')[-1])
        conc = int(row.conc.split('-')[-1])
        OD = int(row.optical_density.split('-')[-1])
        assert inh == conc, 'mapping is incorrect'
        assert inh == OD, 'mapping is incorrect'

    print('mapping is verified correct by test001')