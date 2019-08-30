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
