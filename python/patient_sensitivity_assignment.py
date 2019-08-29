'''
This file's intended purpose is to assign patient sensitivity designations based
on the traditional (beatAML) method of 20-80 quantile tail assignment.

'''
import pandas as pd
import sys
import numpy as np

SENSITIVE_QUANTILE = 0.2
RESISTSNT_QUANTILE = 0.8

if __name__ == '__main__':

    _, data_path = sys.argv

    data = pd.read_csv(data_path, low_memory=False)
    data.lab_id = data.lab_id.astype(str)
    data.auc = data.auc.astype(float)

    #[print(f' {x}') for x in data['auc']]
    data2 = data[~pd.isnull(data['auc'])]

    assignments = pd.DataFrame(columns=['inhibitor','auc','lab_id','call'])
    for inhib in data2.inhibitor.unique():

        inhib_dat = data[data['inhibitor'] == inhib]

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
