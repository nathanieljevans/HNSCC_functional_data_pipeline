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
    for panel_name, panel_id in assayids.values:
        if assay_name in panel_name: 
            return panel_id
    
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

    print('before data shape:', data.shape)
    before_note_nonNA = data[~data.note.isna()].shape[0]

    print('adding annotations... ', end='')
    nerr = 0
    success = 0
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
                        ~(ann.annotations in data.note), 'note'] = data.loc[(data.panel_id == assay_id) & 
                        (data.plate_num == ann.plate_number) &
                        (data.plate_row >= ann.row_from) &
                        (data.plate_row <= ann.row_to) & 
                        (data.plate_col >= ann.col_from) &
                        (data.plate_col <= ann.col_to) & 
                        ~(ann.annotations in data.note), 'note'].apply(str) + " | " + str(ann.annotations)

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
            success+=1
        except Exception as err: 
            nerr += 1
            print()
            print('-----------------------------------------------------------')
            print('failed to add annotation:')
            print(ann)
            print()
            print('Exception traceback: ')
            print(err)
            print('-----------------------------------------------------------')

    after_note_nonNA = data[~data.note.isna()].shape[0]
    new_annotations = data[[('|' in str(x)) for x in data.note]].shape[0]
    print(f'{new_annotations} annotations added [{success} inputs].')
    print('Number of non-NA notes before annotating:', before_note_nonNA)
    print('Number of non-NA notes after annotating:', after_note_nonNA)

    print('removing "unnamed" columns...')
    data = data.loc[:, ~data.columns.str.contains('^Unnamed')]

    print('adding "summary" column for QC criteria (to remove or not - does not include manual additions')

    '''
    for i,x in data.iterrows(): 
        print()
        print('PAC', not x.low_PAC_flag)
        print('within-repl', not x.is_within_plate_repl==True) 
        print('across-repl', not x.is_across_plate_repl==True) 
        print('across-flag val:', x.across_plate_repl_flag)
        print('across-flag', not x.across_plate_repl_flag==True)
        print('AIC-flag', not x.AIC_flag)
        print('DEV-flag', not x.DEV_flag)
        print('overfit-flag', not x.overfit_flag)
    '''

    data = data.assign(passed_QC = [(not x.low_PAC_flag==True and not x.is_within_plate_repl==True and not x.is_across_plate_repl==True and not x.across_plate_repl_flag==True and not x.AIC_flag==True and not x.DEV_flag==True and not x.overfit_flag==True) for i,x in data.iterrows()])

    print('number of annotations that failed:', nerr)
    print('after data shape:', data.shape)

    print('saving data...', end='')
    data.to_csv(args.output_path, index=False)
    print('data saved.')

    print('script completed.')