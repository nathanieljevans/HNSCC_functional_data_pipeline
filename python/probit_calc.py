'''
@author: Nate Evans
@date: 7/16/2019
@title: HNSCC Plate Mapper

Takes a csv file with [lab_id, inhibitor, conc, avg_opt_density]

Assumes normalization steps have already been applied to dataset

XXX

Run the file by
    $ python probit_calc.py ./path/to/data.csv

results are saved to ./data/single_drug_probit_fit_results.csv


TODO
- Perfect separation fails, need work around


'''

import sys
import os
import pandas as pd
import statsmodels.api as sm
import numpy as np
import seaborn as sbn
from matplotlib import pyplot as plt

# comment this to False - will plot failures
DIAGNOSTICS = True

if __name__ == '__main__':

    path = sys.argv[1]

    print('Data path: %s' %path)

    data = pd.read_csv(path, sep=',')

    failures = []
    res = {x:[] for x in ['lab_id', 'inhibitor','beta0', 'beta1', 'auc']}
    i = 0
    for patient in set(data['lab_id'].values):
        pat_dat = data[data['lab_id'] == patient]
        print('Fitting probit | n = %d | lab_id = %s | n_failures = %d' %(i, patient, len(failures) ), end='\r')
        for inhib in set(pat_dat['inhibitor'].values):
            i+=1

            df = pat_dat[pat_dat['inhibitor'] == inhib]

            assert df.shape[0] == 7, 'wrong number of doses'

            try:
                x = sm.add_constant(np.log10(df['conc'].values))
                y = df['avg.opt.density'].values

                pr = sm.GLM(y, x, family=sm.families.Binomial(link=sm.families.links.probit()))
                glm_res = pr.fit(disp=False)

                # AUC calculation -----------------------------------------------------
                # left rectangle auc estimate
                delta = 0.001
                x2 = np.arange(np.log10(min(df['conc'].values)), np.log10(max(df['conc'].values)), delta)
                yhat = glm_res.predict(sm.add_constant(x2))
                auc = np.sum(yhat*delta)

                # beta0 = intercept
                # beta1 = slope
                (beta0,beta1) = glm_res.params

                # update results
                [res[var].append(val) for var,val in zip(['lab_id', 'inhibitor','beta0', 'beta1', 'auc'], [patient, inhib, beta0 ,beta1 , auc])]

            except:
                plt.figure()
                plt.plot(df['conc'].values, df['avg.opt.density'].values, 'b-')
                plt.show()

                failures.append( (patient, inhib) )
                [res[var].append(val) for var,val in zip(['lab_id', 'inhibitor','beta0', 'beta1', 'auc'], [patient, inhib, 'NA' ,'NA' , 'NA'])]
                if DIAGNOSTICS:
                    print('FAILURE: %s, %s' %(patient, inhib))
                    print(df.head(7))

                    f, ax = plt.subplots(1,1, figsize = (10,10))
                    ax.set_title('FAILURE: %s, %s' %(patient, inhib))
                    #sbn.scatterplot(x=x2, y=yhat , ax=ax)
                    plt.xscale('log')
                    sbn.scatterplot(x=np.log10(df['conc'].values), y=df['avg.opt.density'].values, ax=ax)
                    plt.show()

    print('Failures [%d]: %r' %(len(failures),failures))

    res = pd.DataFrame(res)
    print(res.head())

    res.to_csv('./data/single_drug_probit_fit_results.csv')
