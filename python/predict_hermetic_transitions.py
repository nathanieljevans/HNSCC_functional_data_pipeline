'''
The purpose of this script is to use `atypical_doseresponse_classifier` to
predict hermetic transition points in this data. More information can be found
at: https://github.com/nathanieljevans/atypical_doseresponse_classifier
'''

import pandas as pd
import sys
import numpy as np
from keras.models import load_model

if __name__ == '__main__':

    # data in
    _, data_path, model_path = sys.argv

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


'''
## DEPRECATED ##
REGRESSION MODEL VVV
if __name__ == '__main__':

    _, data_path, model_path = sys.argv

    data = pd.read_csv(data_path, low_memory=False)
    in_shp = data.shape

    in_shp = data.shape

    data.lab_id = data.lab_id.astype(str)
    data.auc = data.auc.astype(float)

    herm_df = data[['lab_id', 'inhibitor', 'conc_norm', 'cell_viab', 'panel_id']]
    herm_df = herm_df.assign(conc_norm = lambda x: np.round(x.conc_norm, 4))

    herm_df.conc_norm = herm_df.conc_norm.astype(str)
    #herm_wide = herm_df.piv(['lab_id', 'inhibitor', 'panel_id']).unstack()

    herm_wide = herm_df.pivot_table(index = ['lab_id', 'inhibitor', 'panel_id'], columns='conc_norm', values='cell_viab', aggfunc=np.mean)
    herm_wide = herm_wide.reset_index()

    good_doses = [str(x) for x in reversed([np.round(10/(3**i), 4) for i in range(0,7)])]

    herm_wide = herm_wide.dropna(subset=good_doses)

    herm_wide = herm_wide[good_doses + ['lab_id', 'inhibitor', 'panel_id']]

    #print (good_doses)

    #print(herm_wide.head())

    #print( herm_wide.columns )

    #print( herm_wide[good_doses + ['lab_id', 'inhibitor', 'panel_id']].head(50) )

    X = herm_wide[good_doses].values

    model = load_model(model_path)

    #print('X ---------')
    #print(X[0,:])
    #print('-----------')

    yhat = model.predict(X)

    herm_wide = herm_wide.assign(hermetic_transition = yhat)

    #print(herm_wide.head(50))

    dat2 = herm_wide[ ['lab_id', 'inhibitor', 'panel_id', 'hermetic_transition'] ]

    data = data.merge(dat2, how='left', on=['lab_id', 'inhibitor', 'panel_id'])
    out_shp = data.shape

    print(f'input shape: {in_shp}\noutput shape: {out_shp}')

    data.to_csv(data_path)

'''
