'''
This script adds annotations from `../annotations/error_annotation.xlsx` to the summary output file `HNSCC_all_functional_data.csv`. 

For more information, please refer to `../annotations/readme.md` 
'''


import pandas as pd 
import argparse
import sys


def get_assay_id(assay_name, assayids): 
    '''
    '''
    for i, ass in assayids.iterrows():
        if assay_name in ass[0]: 
            return ass[1]
    
    print()
    print('Assay ID not found:', assay_name) 
    print('Please double check that the assay name is correct and has been processed.')
    return None


if __name__ == '__main__': 

    parser = argparse.ArgumentParser(description='This script adds annotations from `../annotations/error_annotation.xlsx` to the summary output file `HNSCC_all_functional_data.csv`. \nFor more information, please refer to `../annotations/readme.md`')
    parser.add_argument('--no-verbose', dest='verbose', action='store_const',
                       const=False, default=True,
                       help='whether to print output messages to console')
    parser.add_argument('--annotation-path', dest='annotation_path', metavar='<pth>', type=str,
                        help='path to annotations', default='../annotation/error_annotation.xlsx')
    parser.add_argument('--data-path', dest='data_path', metavar='<pth>', type=str,
                        help='path to all data output file', default='../output/HNSCC_all_functional_data.csv')
    parser.add_argument('--assayID-path', dest='assayid_path', metavar='<pth>', type=str,
                        help='path to .assayid file', default='../.assay_ids')
    parser.add_argument('--output-path', dest='output_path', metavar='<pth>', type=str,
                        help='path to output file', default='../output/HNSCC_all_functional_data-annotated.csv')

    args = parser.parse_args() 

    if args.verbose: print('annotation file path:', args.annotation_path)
    if args.verbose: print('output data file path:', args.data_path)

    annotations = pd.read_excel(args.annotation_path)
    data = pd.read_csv(args.data_path, low_memory=False)
    assayids = pd.read_csv(args.assayid_path, sep=' : ',header=None, engine='python')

    data = data.assign(manual_flag=False)
    data = data.assign(manual_remove=False)

    before_note_nonNA = data[~data.note.isna()].shape[0]

    print('adding annotations... ', end='')
    for idx, ann in annotations.iterrows():
        try: 
            assay_id = get_assay_id(ann.assay_name, assayids)
            assert assay_id is not None, 'No panel ID found.'
            data.loc[(data.panel_id == assay_id) & 
                        (data.plate_num == ann.plate_number) &
                        (data.plate_row >= ann.row_from) &
                        (data.plate_row <= ann.row_to) & 
                        (data.plate_col >= ann.col_from) &
                        (data.plate_col <= ann.col_to) & 
                        ~(ann.annotations in data.note), 'note'] = ann.annotations
            data.loc[(data.panel_id == assay_id) & 
                        (data.plate_num == ann.plate_number) &
                        (data.plate_row >= ann.row_from) &
                        (data.plate_row <= ann.row_to) & 
                        (data.plate_col >= ann.col_from) &
                        (data.plate_col <= ann.col_to), 'manual_remove'] = ann.remove
            data.loc[(data.panel_id == assay_id) & 
                        (data.plate_num == ann.plate_number) &
                        (data.plate_row >= ann.row_from) &
                        (data.plate_row <= ann.row_to) & 
                        (data.plate_col >= ann.col_from) &
                        (data.plate_col <= ann.col_to), 'manual_flag'] = True
        except Exception as err: 
            print()
            print('-----------------------------------------------------------')
            print('failed to add annotation:')
            print(ann)
            print()
            print('Exception traceback: ')
            print(err)
            print('-----------------------------------------------------------')

    after_note_nonNA = data[~data.note.isna()].shape[0]

    print(f'{after_note_nonNA - before_note_nonNA} annotations added.')

    print('saving data...', end='')
    data.to_csv(args.output_path)
    print('data saved.')

    print('script completed.')