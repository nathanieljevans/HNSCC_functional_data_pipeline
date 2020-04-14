'''
This script takes the ouput data and: 
1. Aggregates 7-point dose-response observations into a single observation
2. Filters observations based on quality control metrics

'''

# TODO : doc testing... 

import pandas as pd 
import argparse
import sys
import numpy as np

def QC(assay): 
    '''

    '''
    assay = assay[assay.low_PAC_flag != True]               #? Remove assays with plate average control less than ____
    assay = assay[assay.is_within_plate_repl != True]       #? Remove assays that are aggregated within plate - how to test that this is done correctly?
    assay = assay[assay.is_across_plate_repl != True]       #? Remove assays that have been aggregated across plates - how to test? 
    assay = assay[assay.across_plate_repl_flag != True]      #! Remove assays that have been aggregated across plates and have a difference greater than 1. 
    assay = assay[assay.AIC_flag != True]                   #? Remove assays that have AIC value of ____ 
    assay = assay[assay.DEV_flag != True]                   #? Remove assays that have Deviance value of ____ 
    assay = assay[assay.overfit_flag != True]               #? Remove assays that have probit_BIC > poly_BIC
    assay = assay[assay.manual_remove != True]              #? Remove assays that have been manually annotated for removal

    #TODO: add in functionality to remove full well range if only *some* of the wells were manually flagged
    #if assay.shape[0] % 7 != 0: 
    #    print('----- adfd-----')
    #    print(f'there are wrong number of doses [{assay.shape[0]}] for assay: {(assay.lab_id.values[0], assay.inhibitor.values[0])}')
    #    print(assay)

    assert assay.shape[0] % 7 == 0, f'there are wrong number of doses [{assay.shape[0]}] for assay: {(assay.lab_id.values[0], assay.inhibitor.values[0])}'

    return assay

def format_plate_location( assay ): 
    '''
    '''
    panel = assay['panel_id'].unique()
    plate = assay['plate_num'].unique()
    row = assay['plate_row'].unique()
    col_min = assay['plate_col'].min()
    col_max = assay['plate_col'].max()
    
    return f"{' '.join(panel.astype(str))}:{' '.join(plate.astype(str))}:{' '.join(row.astype(str))}:{col_min}-{col_max}"

def get_max_conc( assay ): 
    '''
    '''
    conc = assay['conc'].values
    if ';' in conc[0]: 
        temp = [float(c.split(';')[0]) for c in conc]
        conc_max = conc[np.argmax(temp)]
    else:
        conc_max = conc.astype(np.float).max()
    return conc_max


if __name__ == '__main__': 

    parser = argparse.ArgumentParser(description='This script This script takes the ouput data and:\n 1. Aggregates 7-point dose-response observations into a single observation \n2. Filters observations based on quality control metrics')
    parser.add_argument('--no-verbose', dest='verbose', action='store_const',
                       const=False, default=True,
                       help='whether to print output messages to console')
    parser.add_argument('--data-path', dest='data_path', metavar='<pth>', type=str,
                        help='path to all data output file', default='../output/HNSCC_all_functional_data-annotated.csv')
    parser.add_argument('--output-path', dest='output_path', metavar='<pth>', type=str,
                        help='path to output file', default='../output/HNSCC_cleaned_data.csv')
                        
    args = parser.parse_args() 

    data = pd.read_csv(args.data_path, low_memory=False)
    
    lab_ids = data.lab_id.unique() 
    inhibitors = data[~data.inhibitor.isin(['F-S-V','DMSO','NONE'])].inhibitor.unique()

    print('\n-------------------------------------------------------------------------------------')
    print('Before filtering and aggregation...')
    print(f'\tnumber of inhibitors:\t\t\t {len(inhibitors)}')
    print(f'\tnumber of patients:\t\t\t {len(lab_ids)}')
    print(f'\tnumber of dose-response assays:\t\t {data[["inhibitor", "lab_id"]].drop_duplicates().shape[0] }')
    print('-------------------------------------------------------------------------------------')

    print()
    nerr = 0
    cleaned = {'lab_id':[], 'inhibitor':[], 'AUC':[], 'plate_loc':[], 'flagged':[], 'max_conc':[], 'call':[], 'replicates':[]}
    for i, inhib in enumerate(inhibitors): 
        inhib_data = data[data.inhibitor == inhib]
        print(f'progress: {100*i/len(inhibitors):.2f}% - processing:   {inhib} \t\t', end='\r')
        for lab_id in lab_ids: 
            try: 
                assay = inhib_data[inhib_data.lab_id == lab_id]
                repl = assay.shape[0] / 7
                assay = QC(assay)

                # TODO: fix all the nans - this is happening in the aggregation step across/within plates I think. 
                if assay.plate_num.isna().any(): 
                    continue

                if assay.shape[0] == 0: continue 

                aucs = assay.auc.unique()
                if len( aucs ) > 1: 

                    # TODO: This should have been aggregated earlier in the process - why is it not? 
                    #print()
                    #print(f'User Warning - aggregating AUC ({lab_id}, {inhib}); may be an issue in processing pipeline')
                    #print('AUCS:', aucs)
                    #print('repls', repl)
                    #print(assay)
                    assert len(aucs) == repl, f'There are {repl} replicates but only {len(aucs)} unique auc values'
                    auc = np.mean(aucs)
                else: 
                    auc = aucs[0]

                plate_loc = format_plate_location(assay)
                
                cleaned['lab_id'].append( lab_id )
                cleaned['inhibitor'].append( inhib )
                cleaned['AUC'].append( auc )
                cleaned['plate_loc'].append( plate_loc )
                cleaned['flagged'].append( assay['manual_flag'].any() )
                cleaned['max_conc'].append( get_max_conc(assay) )
                cleaned['call'].append( ' '.join(assay['call'].unique()) )
                cleaned['replicates'].append( repl )

            except Exception as err: 
                nerr += 1
                print('\t\t-----------------------------------------------------------------------------')
                print('\t\tError |     lab_id           inhibitor                         plate_loc')
                print(f'\t\t           {lab_id}         {inhib}       \t\t\t  {plate_loc}')
                print(f'\t\t{err}')
                print('\t\t-----------------------------------------------------------------------------')
                #print(assay)
                #raise

    cleaned = pd.DataFrame(cleaned)

    print('\n\n-------------------------------------------------------------------------------------')
    print('After filtering and aggregation...')
    print(f'\tnumber of inhibitors:\t\t\t {len(cleaned.inhibitor.unique())}')
    print(f'\tnumber of patients:\t\t\t {len(cleaned.lab_id.unique())}')
    print(f'\tnumber of dose-response assays:\t\t {cleaned[["inhibitor", "lab_id"]].drop_duplicates().shape[0] }')
    print('\tNumber of errors:', nerr)
    print('-------------------------------------------------------------------------------------')

    print('saving cleaned data to file...', end='')
    cleaned.to_csv(args.output_path, index=False )
    print('done')