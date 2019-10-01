'''


'''

from HNSCC_analysis_pipeline_lib import *
import sys
import os
from datetime import datetime as dt

DO_RAISE = True

if __name__ == '__main__':
    tic = dt.now()
    _, map_dir, data_dir = sys.argv

    dirlist = os.listdir(data_dir)
    failures = []
    errors = []

    for i,p in enumerate(dirlist):
        try:
            print('Panels processed: [%d/%d]' %(i, len(dirlist)))
            process(data_dir + '/' + p, platemap_dir=map_dir, verbose=False, do_raise=DO_RAISE)
        except Exception as e:
            print('Processing failed. See error log for details.')
            failures.append(p)
            errors.append(e)
            if DO_RAISE: raise

    print('There were %d failures. \n\t %r' %(len(failures), failures))
    [print(str(e)) for e in errors]

    print(f'Total time elapsed: {dt.now() - tic}')
