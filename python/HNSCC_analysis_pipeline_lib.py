'''
@author: Nate Evans
@date: 7/16/2019
@title: HNSCC Functional Data Analysis Pipeline Library

The purpose of this script is to provide the resources to enable an automated
analysis pipeline of functional drug response data in the HNSCC project at
OHSU led by Shannon McWeeney and Molly Kulesz-Martin.

'''

import pandas as pd
import numpy as np
from datetime import datetime as dt
from matplotlib import pyplot as plt
import seaborn as sbn
import statsmodels.api as sm
import statsmodels
import os
import shutil
from sklearn.preprocessing import PolynomialFeatures
import warnings

def parse_pathname(path, verbose=False):
    '''
    parse relevant information from data path, must return
        - lab_id [unique 5 digit patient identifier]
        - norm [Blank490, 490]
        - version_id <str>

    file name must be in form:
        lab_id=XXXXX-norm=XXXXXX-plate_version_id=XXXXX.xlsx

    '''

    name = path.split('/')[-1][:-5].split('-')

    labels = [n.split('=')[0] for n in name]
    values = [n.split('=')[-1] for n in name]

    assert len(labels) == 4, 'parsing error: atypical pathname - labels require: lab_id, norm, version_id, note'
    assert len(values) == 4, 'parsing error: atypical pathname - values require 4 inputs for each label.'

    if verbose: lab_id, norm, version_id, notes = values
    if verbose: print('---------------------------------------------------------')
    if verbose: print( 'please double check file name parsing is accurate')
    if verbose: print( 'file name: %s' %data_path.split('/')[-1])
    if verbose: print('lab_id: %s' %lab_id)
    if verbose: print('norm: %s' %norm)
    if verbose: print('version_id: %s' %version_id)
    if verbose: print('notes: %s' %notes)
    if verbose: print('---------------------------------------------------------')

    return values

# ------------------------------------------------------------------------------
# deprecated
def get_plate_data(data_path, verbose=False):
    '''
    DEPRECATED
    each plate has 16 rows of drug data, then a empty row, then the next plate. Varying number of plates

    inputs
        data_path <str> path to single patient plate data

    outputs
        dataframe
    '''
    print('this method is deprecated, use panel._get_plate_data() in the future')

    lab_id, norm, version_id, notes = parse_pathname(data_path)

    if verbose: print('---------------------------------------------------------')
    if verbose: print( 'please double check file name parsing is accurate')
    if verbose: print( 'file name: %s' %data_path.split('/')[-1])
    if verbose: print('lab_id: %s' %lab_id)
    if verbose: print('norm: %s' %norm)
    if verbose: print('version_id: %s' %version_id)
    if verbose: print('notes: %s' %notes)
    if verbose: print('---------------------------------------------------------')

    allplates = pd.read_excel(data_path, header=None)

    nplates = (allplates.shape[0] + 1) / 18

    assert nplates%1.0 == 0, 'atypical number of rows, check data format'

    if verbose: print('assay has %d plates' %int(nplates))

    plates = []
    i = 1 # skip original header
    warned = False
    for p in range(int(nplates)):
        dat = pd.DataFrame( allplates.values[i:(i+16),:])

        # check for blank490 column ... confidence in proper documentation
        if dat.shape[1] < 26:
            if not warned: print('WARNING: This assay [lab_id=%s,notes=%s] does not have a "blank490" column (last col), please double check that the data has been normalized by the positive controls.' %(lab_id,notes))
            warned = True
            dat = dat.assign(plate_row = dat[0]).assign(norm_type = 'none', plate_num = p+1, lab_id = lab_id, assay_version_id=version_id, note=notes).drop(labels = [0], axis='columns')
        else:
            dat = dat.assign(plate_row = dat[0]).assign(norm_type = dat[25], plate_num = p+1, lab_id = lab_id, assay_version_id=version_id, note=notes).drop(labels = [0,25], axis='columns')

        dat = pd.melt(dat, id_vars=['lab_id', 'norm_type', 'plate_num', 'plate_row','assay_version_id', 'note'], value_vars=None, var_name='plate_col', value_name='optical_density', col_level=None)

        plates.append( dat )
        i += 16 + 2 # skip empty row + header

    # Frame = Frame.append(pandas.DataFrame(data = SomeNewLineOfData), ignore_index=True)
    df = plates[0]
    for p in plates[1:]:
        df = df.append(p, ignore_index=True)

    return df

# ------------------------------------------------------------------------------

def get_plate_map(map_path, verbose=False):
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

class panel:
    '''
    This class represents a single HNSCC drug panel.
    '''
    def __init__(self, plate_path, platemap_dir, verbose=True):
        '''
        parse path name, load plate data, load plate map, check positive control
        '''

        self.plate_path = plate_path
        self.platemap_dir = platemap_dir
        self.verbose = verbose

        lab_id, norm, version_id, notes = parse_pathname(plate_path)

        self.msg_log = ''
        self._log('\n--------------------------------------------------------------------------\nThis is the message log for: \n\tlab_id: %s\n\tversion_id: %s\n\tnotes: %s\n--------------------------------------------------------------------------' %(lab_id, version_id, notes))

        assert len(lab_id) <= 5, 'lab_id is expected to be 5 or less characters long'
        if len(lab_id) < 5: 
            print('Warning: lab_id is less than 5 characters long; appending 0 prefix')
            lab_id = '0'*(5-len(lab_id)) + lab_id
        self.lab_id = lab_id

        #assert version_id in ['OHSU_HNSCC_derm001', 'OHSU_HNSCC_derm002'], 'unknown platemap ID, please double check pathname or update script with acceptable plate IDs'
        self.version_id = version_id

        assert norm in ['blank490', '490', 'Blank490'], 'improper pathname value for "norm" [blank490, 490]'
        self.blank490 = (norm in ['blank490', 'Blank490'])
        self._log('Data has already been normalized by positive controls: %s' %str(self.blank490))

        self.notes = notes

        self._raw = self._get_plate_data()

        platemap_path = self.platemap_dir + 'HNSCC_plate_map-version_id=' + self.version_id + '.xlsx'
        self.platemap = get_plate_map(platemap_path)

    def _get_plate_data(self):
        '''
        each plate has 16 rows of drug data, then a empty row, then the next plate. Varying number of plates

        file name specifies if the spectrometer software has already normalized the data by positive controls
        (eg. set optical_density = 0 when all cells are dead). This is noted by norm='Blank490'
        To double check this however, we can make sure that there is a 26th column of the raw plate data
        which is filled with 'blank490' if this is the case and empty otherwise.

        inputs
            None

        outputs
            dataframe
        '''

        allplates = pd.read_excel(self.plate_path, header=None)

        nplates = (allplates.shape[0] + 1) / 18

        assert nplates%1.0 == 0, 'atypical number of rows, check data format'

        self._log('This assay has %d plates' %int(nplates))

        plates = []
        i = 1 # skip original header
        warned = False
        for p in range(int(nplates)):
            dat = pd.DataFrame( allplates.values[i:(i+16),:])

            # check for blank490 column ... confidence in proper documentation
            if dat.shape[1] < 26:
                if not warned: self._log('WARNING: This assay [lab_id=%s,notes=%s] does not have a "blank490" column (last col), please double check that the data has been normalized by the positive controls.' %(self.lab_id,self.notes))
                warned = True
                if self.blank490: self._log('WARNING: file name specifies blank490, but data does not appeaer to be postive control normalized.')
                self.blank490=False
                dat = dat.assign(plate_row = dat[0]).assign(norm_type = 'none', plate_num = p+1, lab_id = self.lab_id, assay_version_id=self.version_id, note=self.notes).drop(labels = [0], axis='columns')
            else:
                dat = dat.assign(plate_row = dat[0]).assign(norm_type = dat[25], plate_num = p+1, lab_id = self.lab_id, assay_version_id=self.version_id, note=self.notes).drop(labels = [0,25], axis='columns')

            dat = pd.melt(dat, id_vars=['lab_id', 'norm_type', 'plate_num', 'plate_row','assay_version_id', 'note'], value_vars=None, var_name='plate_col', value_name='optical_density', col_level=None)

            plates.append( dat )
            i += 16 + 2 # skip empty row + header

        df = plates[0]
        for p in plates[1:]:
            df = df.append(p, ignore_index=True)

        return df

    def _log(self, msg):
        '''

        '''
        if self.verbose: print(msg)
        self.msg_log += '>> ' + str(dt.now()) + ': ' + msg + '\n'

    def map_data(self):
        '''

        '''
        self._log('mapping data... [%s]' %self.version_id)
        self.data = self._raw.merge(self.platemap, how='left', left_on=['plate_row','plate_col','plate_num'], right_on=['row','col','plate_number']).drop(['row','col','plate_number'], axis='columns')

        # remove any leading or trailing spaces
        self.data = self.data.assign(inhibitor = [x.strip() for x in self.data.inhibitor])

        self._log('mapping complete. \n\t mapped data shape: %s\n\t column names: [%r]' %(str(self.data.shape), self.data.columns.tolist()))

    def normalize_combinationagent_concentrations(self):
        '''
        normalize combination drug assays concentrations by mapping into euclidean distance.

        conc_norm = ( conc1**2 + conc2**2 )**(0.5) if combination_agent else conc

        '''
        self._log('normalizing combination agent concentrations...')

        self.data = self.data.assign(iscomb= [';' in str(conc) for conc in self.data['conc']])

        self.data = self.data.assign(conc_norm = [ [np.sqrt(np.sum(float(y)**2 for y in str(x).split(';')))].pop() for x in self.data['conc']])


    def normalize_cell_viability_by_negative_controls(self):
        '''
        The zero value (positive control) of optical density is set by the p.spec software [blank490] but to interpret assay value as cell viability, we need to scale optical density by plate controls:
        Cell_viab = opt_dens / plate_average_control (PAC)
        [This sets PAC control value as cell viability of 1]
        '''
        self._log('normalizing cell viability by negative controls...')

        # double check the data type is numeric
        self.data = self.data.assign(optical_density = [float(x) for x in self.data.optical_density])

        # get plate average controls
        PAC = self.data[self.data['inhibitor'] == 'NONE'].groupby(by='plate_num')['optical_density'].mean().to_frame().assign(PAC = lambda x: x.optical_density).drop(['optical_density'], axis='columns')

        self._log('plate average controls: \n%s' %str(PAC))

        self.data = self.data.merge(PAC, how='left', on='plate_num')

        self.data = self.data.assign(cell_viab = [od / pac for od, pac in zip(self.data['optical_density'], self.data['PAC'])])

        # set a flag for low_pac (0.03)
        self.data = self.data.assign(low_PAC_flag = self.data.PAC < 0.03)

    def set_floor(self):
        '''
        There is still some contention here, Dan's protocol shifts the PAC up by the smallest value on plate.
        [and therefore this step would need to come before negative control normalization]

        I don't think that's a good idea though, rather, I'm going to round any negatives up to 0. These will be marked with 'is_adj'
        FOR NOW - talk to the group
        '''
        self.data = self.data.assign(is_adj = self.data.cell_viab < 0, cell_viab = [0 if cv < 0 else cv for cv in self.data.cell_viab] )


    def replicate_QC(self, method=['within','across'], flag_threshold = 1):
        '''
        linear regression model fit to each within-plate replicate separately
        and AUC values are calculated. Replicates with AUC differences greater
        than 100 are flagged. ?flag_name? [TODO] replicates are not normalized,
        future probit calculations should use all data to fit the probit curve.

        inputs
            flag_threshold <float> : threshold to set flag for QC removal
            method <str> : whether to average within plate or across plates.

        outputs
            none
        '''



    ### DEPRECATED ###
    def avg_plate_replicates(self, method=['within','across'], flag_threshold = 1):
        '''
        A ‘curve-free’ AUC (integration based on fine linear interpolation between
        the 7 data points themselves) was calculated for those runs with within-panel
        replicates after applying a ceiling of 100 for the normalized_viability.
        The maximum change in AUC amongst the replicates was noted and those runs
        with differences > 1 were removed.

        Across-plate replicates are fit with linear regression and replicates with AUC differences > 0.75
         are removed. Other replicates are averaged.


        The HNSCC plates don't appear to have any within plate replicates, so this is an essentially unused method.
        As such, it has been only lightly tested. Be wary of performance and robustness.

        inputs
            flag_threshold <float> : threshold to set flag for QC removal
            method <str> : whether to average within plate or across plates.

        outputs
            none
        '''
        warnings.warn('This method is deprecated and will be removed in future release. Use `replicate_QC()`` in the future.')

        n = 0

        replicates = []
        toflag = []

        fitdat = self.data[~self.data['inhibitor'].isin(['NONE', 'F-S-V', 'DMSO'])]

        if method == 'within':
            avgd_obs = []
            self._log('averaging within plate replicates...')

            for plate in set(fitdat['plate_num'].values):
                platedat = fitdat[fitdat['plate_num'] == plate]

                for inhib in set(platedat['inhibitor'].values):
                    assay = platedat[platedat['inhibitor'] == inhib]

                    if assay.shape[0] > 7:
                        n+=1

                        aucs = []
                        # groupby row - TODO: should double check that no replicates are on same row
                        for row in set(assay.plate_row):
                            repl = assay[assay['plate_row']==row]

                            aucs.append(self._get_lin_auc(repl))

                        if max(aucs) - min(aucs) > flag_threshold: toflag.append(inhib)

                        replicates.append(inhib)
                        assay.groupby(['lab_id','inhibitor','conc_norm', 'plate_num'])['cell_viab'].mean().to_frame().reset_index(level=['lab_id', 'inhibitor', 'conc_norm', 'plate_num']).assign(within_plate_repl_flag = (max(aucs) - min(aucs)) > flag_threshold)
                        avgd_obs.append(new_obs)

            self._log('There were %d within plate replicates [%r]' %(n, replicates))
            self._log('within plate replicates flagged: %r' %toflag)

            if len(replicates) > 0:
                # set is_repl flag, signifies that it has been averaged
                repls = pd.DataFrame({'inhibitor': replicates, 'is_within_plate_repl': [True]*len(replicates)})
                self.data = self.data.merge(repls, how='left', on='inhibitor').assign(is_within_plate_repl = lambda x: [True if f else False for f in x.is_within_plate_repl])
                # then add the averaged replicate observation to the bottom of dataframe
                for nobs in avgd_obs:
                    self.data = self.data.append(nobs, ignore_index=True)

            else:
                self.data = self.data.assign(is_within_plate_repl = False)

            # add auc diff flag for removal
            if (len(toflag) > 0):
                flags = pd.DataFrame({'inhibitor': toflag, 'within_plate_repl_flag': [True]*len(toflag)})
                self.data = self.data.assign(within_plate_repl_flag = lambda x: [True if f else False for f in x.withinplate_repl_flag]).merge(flags, how='left', by='inhibitor')
            else:
                self.data.assign(within_plate_repl_flag = False)

        # ------------- Across plate replicate averaging VVV ------------------

        elif method == 'across':
            self._log('averaging across plate replicates...')
            avgd_obs = []
            replicates = []

            for inhib in set(fitdat['inhibitor'].values):
                assay = fitdat[fitdat['inhibitor'] == inhib]

                #print('assay shape: %s' %str(assay.shape))
                if assay.shape[0] > 7:
                    n+=1

                    aucs = []
                    # groupby row - TODO: should double check that no replicates are on same row
                    for plate in set(assay.plate_num):
                        repl = assay[assay['plate_num']==plate]

                        aucs.append(self._get_lin_auc(repl))

                    self._log('[%s] linear AUCs: %r -- range: (%.2f, %.2f)' %(inhib, aucs, min(aucs), max(aucs)))
                    if (max(aucs) - min(aucs)) > flag_threshold: toflag.append(inhib)

                    replicates.append(inhib)
                    new_obs = assay.groupby(['lab_id','inhibitor','conc_norm'])['cell_viab'].mean().to_frame().reset_index(level=['lab_id', 'inhibitor', 'conc_norm']).assign(across_plate_repl_flag = (max(aucs) - min(aucs)) > flag_threshold)
                    avgd_obs.append(new_obs)

            self._log('There were %d across plate replicates [%r]' %(n, replicates))
            self._log('across plate replicates flagged: %r' %toflag)

            if len(replicates) > 0:
                # set is_repl flag, signifies that it has been averaged
                repls = pd.DataFrame({'inhibitor': replicates, 'is_across_plate_repl': [True]*len(replicates)})
                self.data = self.data.merge(repls, how='left', on='inhibitor')
                self.data = self.data.assign(is_across_plate_repl = lambda x: [False if f != True else True for f in x.is_across_plate_repl])
                # then add the averaged replicate observation to the bottom of dataframe
                for nobs in avgd_obs:
                    self.data = self.data.append(nobs, ignore_index=True, sort=False)
            else:
                self.data = self.data.assign(is_across_plate_repl = False)

            if (len(toflag) > 0):
                flags = pd.DataFrame({'inhibitor': toflag, 'across_plate_repl_flag': [True]*len(toflag)})
                self.data = self.data.assign(across_plate_repl_flag = lambda x: [True if f else False for f in x.across_plate_repl_flag]).merge(flags, how='left', on='inhibitor')
            else:
                self.data.assign(across_plate_repl_flag = False)

        else:
            raise 'choose a proper averaging method [within, across] plates'
    
    def _get_rectangle_auc(self, x, y): 
        '''
        This is used when regression fit fails due to perfect separation. Using left-shifted rectangles.  
        
        x should be log10 transformed concentration 
        y should be cell viab [0,1]
        '''
        delta = (np.max(x) - np.min(x))/(len(y)-1)
        auc = np.sum([delta * yy for yy in  y[1:]])
        #print(f'rectangular auc calc: %s' %auc)
        return auc

    def _get_lin_auc(self, df, plot=False, return_fig=False):
        '''
        fit a linear regression to the data and calculate auc

        inputs
            df <dataframe> 7 points to plot
            plot <bool> option to plot the fit; for diagnostics

        outputs
            auc <float> area under the curve
        '''
        assert df.shape[0] == 7, 'There should only be 7 observations per replicate'
        skip_auc_calc=False
        
        x = [np.log10(x) for x in df['conc_norm'].values]
        y = df['cell_viab'].values

        pr = sm.GLM(y, sm.add_constant(x))
        
        try: 
            glm_res = pr.fit()
        except ValueError: 
            try: 
                glm_res = pr.fit(method='newton')
            except sm.tools.sm_exceptions.PerfectSeparationError: 
                auc = self._get_rectangle_auc(x,y)
                self._log('Perfect separtion in linear regression: auc calculated by rectangular approximation. AUC=%.3f' %auc)
                skip_auc_calc=True

        # AUC calculation -----------------------------------------------------
        # left rectangle auc estimate
        if not skip_auc_calc: 
            delta = 0.001
            x2 = np.arange(np.log10(min(df['conc_norm'].values)), np.log10(max(df['conc_norm'].values)), delta)
            yhat = glm_res.predict(sm.add_constant(x2))
            auc = np.sum(yhat*delta)

            f = plt.figure(figsize = (10,10))
            plt.title('Linear Regression [AUC=%.2f]' %auc)
            plt.plot(x, y, 'ro', label='replicate')
            plt.plot(x2, yhat, 'g-', label='fit')
            plt.legend()
            if plot: plt.show()

            if return_fig: return auc, f
            else:
                plt.close(f)
        
        return auc
        
        

    def set_ceiling(self):
        '''
        Apply ceiling of 1
        (Dan’s protocol uses 100 – note for AUC threshold adjustments)

        inputs
            none

        outputs
            none
        '''
        self._log('Applying a ceiling of 1 to cell viability...')
        self.data = self.data.assign(cell_viab = [1 if cv > 1 else cv for cv in self.data.cell_viab])

    def fit_probit(self, inhib, x, y, df, res, failures, plot=True):
        '''

        '''
        try:
            pr = sm.GLM(y, sm.add_constant(x), family=sm.families.Binomial(link=sm.families.links.probit()))
            glm_res = pr.fit(disp=False)

            auc, x2, yhat = self.calculate_auc(x, glm_res)

            f, ax = plt.subplots(1,1, figsize=(10,10))
            plt.title('inhib: %s [AUC= %.2f]' %(inhib, auc))
            plt.plot(x2, yhat, 'r-', label='probit_fit')
            plt.plot(x, y, 'bo', label='replicates')
            plt.legend()
            if plot: plt.show()
            self._save_plot(f, inhib)
            plt.close(f)

            # beta0 = intercept
            # beta1 = slope
            intercept,beta1 = glm_res.params
            probit_AIC, probit_BIC = glm_res.aic, glm_res.bic
            probit_Deviance = glm_res.deviance
            probit_pval = glm_res.pvalues

            self._write_summary_to_file(glm_res.summary(), 'probit', inhib)

            # update results
            [res[var].append(val) for var,val in zip(['inhibitor','intercept', 'beta1', 'auc', 'prob_conv', 'prob_AIC', 'prob_BIC', 'prob_deviance', 'prob_pval'], [inhib, intercept ,beta1 , auc, True, probit_AIC, probit_BIC, probit_Deviance, probit_pval])]

        except statsmodels.tools.sm_exceptions.PerfectSeparationError as e:
            '''
            Perfect separation occurs when all cell_viab values are identical.
            This most commonly occurs with cell_viab = 1
            '''
            self._log('Perfect separation has occured during probit fitting [%s]\n\t cell_viab values: %r' %(inhib, df['cell_viab'].values))
            cv = df['cell_viab'].unique()
            #assert len(cv) == 1, 'perfect separation has occured, however, there are multiple cell_viab values. Please investigate issue.'

            if (len(cv) == 1):
                self._log('cell_viab values identical; AUC will be calculated by conc-range * cell_viab')
                cv = cv[0]
                auc = cv * (np.log10(max(df['conc_norm'])) - np.log10(min(df['conc_norm'])))
            else:
                self._log('cell_viab values non-identical; AUC valule will be calculated by linear regression.')
                auc, f = self._get_lin_auc(df, plot=plot, return_fig=True)
                self._save_plot(f, inhib)
                plt.close('all')

            [res[var].append(val) for var,val in zip(['inhibitor','intercept', 'beta1', 'auc', 'prob_conv', 'prob_AIC', 'prob_BIC', 'prob_deviance', 'prob_pval'], [inhib, None, None, auc, False, None, None, None, None])]

        except:
            #print('some other exception...')
            failures.append( inhib )
            #[res[var].append(val) for var,val in zip(['inhibitor','intercept', 'beta1', 'auc'], [inhib, 'NA' ,'NA' , 'NA'])]
            if plot:
                f, ax = plt.subplots(1,1, figsize = (10,10))
                ax.set_title('FAILURE: %s' %(inhib))
                #sbn.scatterplot(x=x2, y=yhat , ax=ax)
                plt.xscale('log')
                sbn.scatterplot(x=x, y=y, ax=ax)
                plt.show()
            raise

    def fit_poly(self, inhib, x, y, res, failures, degree=5, plot=True):
        '''

        '''
        try:
            xp = PolynomialFeatures(degree=degree).fit_transform(x.reshape((-1,1)))

            pr = sm.GLM(y, xp)
            poly_res = pr.fit(disp=False)

            x2 = np.arange(min(x), max(x), 0.01)
            x2p = PolynomialFeatures(degree=degree).fit_transform(x2.reshape((-1,1)))
            yhat = poly_res.predict(sm.add_constant(x2p))

            f, ax = plt.subplots(1,1, figsize=(10,10))
            plt.title('POLY FIT (degree=%d) - inhib: %s' %(degree, inhib))
            plt.plot(x2, yhat, 'r-', label='probit_fit')
            plt.plot(x, y, 'bo', label='replicates')
            plt.legend()
            if plot: plt.show()

            self._save_plot(f, inhib, suffix='-poly5')
            plt.close(f)

            poly_AIC, poly_BIC = poly_res.aic, poly_res.bic
            poly_Deviance = poly_res.deviance
            poly_pval = poly_res.pvalues

            self._write_summary_to_file(poly_res.summary(), 'poly5', inhib)

            # update results
            [res[var].append(val) for var,val in zip(['inhibitor', 'poly_degree', 'poly_AIC', 'poly_BIC', 'poly_deviance', 'poly_pval'], [inhib, degree, poly_AIC, poly_BIC, poly_Deviance, poly_pval])]

        except ValueError as e:
            self._log('WARNING! [%s] - %s' %(inhib, str(e)))

        except:
            print('I see the error of my ways :|')
            failures.append( inhib )

            raise

    def _write_summary_to_file(self, summary, model_type, inhib, output_dir='../output/'):
        '''

        '''
        dirout = './%s/%s' %(output_dir, self.plate_path[:-5].split('/')[-1])
        with open('%s/dose-response-plots/%s/GLM_%s_summary.txt' %(dirout, inhib, model_type), 'w') as f:
            f.write(str(summary))

    def calculate_auc(self, x, model, delta=0.001):
        '''
        left rectangle auc estimate. Note that x is already in log10 space.

        inputs
            x <numpy array> endogenous variable
            model <sm.GLMmodel.results> model
            delta <float> the rectangular width to use for each rectangle; smaller = more accurate
        outputs
            auc
        '''
        x2 = np.arange(min(x), max(x), delta)
        yhat = model.predict(sm.add_constant(x2))
        return np.sum(yhat*delta), x2, yhat

    def fit_regressions(self, plot=True):
        '''
        A probit regression was fit to all possible run groups using the model:

        $$\frac{normalized\_viability}{100} \ = \ 1 + log10(concentration) $$

        For all groups there were N=7 dose-response measurements and the following statistics were generated:

        - Regression coefficients labeled as ‘intercept’ and ‘beta’ (slope)
        - Tests of significance for the slope coefficient: *z-statistic* [beta_z] and *p-value* [beta_p]
        - model AIC
        - Summary measures of fit from residuals (Hosmer & Lemeshow 2000)
            - Pearsons chi square
            - Deviance
        - Indication of whether model converged - *fitting method and convergence threshold?*
            - Otherwise values are from last iteration of fitting alg.
        - AUC
            - Area under the curve from numerical integration
        - no_curve_auc
            - Auc from a ‘no-curve’ approach as described in the ‘Preprocessing’ section.
        - linsp_ic50 and linsp_auc
            - IC50 and AUC from an ‘overfit’ curve using a linear 5th degree spline

        inputs
            none

        outputs
            none
        '''
        self._log('fitting probit regressions...')

        failures = []
        probit_res = {x:[] for x in ['inhibitor','intercept', 'beta1', 'auc', 'prob_conv', 'prob_AIC', 'prob_BIC', 'prob_deviance', 'prob_pval']}
        poly_res = {x:[] for x in ['inhibitor', 'poly_degree', 'poly_AIC', 'poly_BIC', 'poly_deviance', 'poly_pval']}

        i = 0

        # Filter controls
        pat_dat = self.data[~self.data['inhibitor'].isin(['NONE', 'F-S-V', 'DMSO'])]

        # Filter replicates that have been averaged into new observation
        pat_dat = pat_dat[pat_dat['is_across_plate_repl'] != True]
        pat_dat = pat_dat[pat_dat['is_within_plate_repl'] != True]

        ntofit = int(pat_dat.shape[0] / 7)
        self._log('number of assays to fit: %d' %ntofit)

        for inhib in set(pat_dat['inhibitor'].values):
            if i%1==0: print('\t\tFitting regressions...Progress: %.1f%% \t[%d/%d]' %(i/ntofit*100, i, ntofit), end='\t\t\t\t\t\r')
            i+=1

            df = pat_dat[pat_dat['inhibitor'] == inhib]

            assert df.shape[0] == 7, 'should have exactly 7 doses! [has %d]' %df.shape[0]

            x = np.log10( df['conc_norm'].values )
            y = df['cell_viab'].values

            # probit_res is an object and so should be passed by reference, hence modified inplace
            self.fit_probit(inhib, x, y, df, probit_res, failures, plot=plot)

            # fit polynomial for comparison to overfitting
            self.fit_poly(inhib, x, y, poly_res, failures, plot=plot)

        self._log('Failures [%d]: \n\t%r' %(len(failures),failures))

        # add probit features
        self.data = self.data.merge(pd.DataFrame(probit_res), how='left', on='inhibitor')

        # add poly features
        self.data = self.data.merge(pd.DataFrame(poly_res), how='left', on='inhibitor')

    def _save_plot(self, fig, inhibitor, output_dir='../output', suffix='-probit'):
        '''

        '''
        try:
            dirout = './%s/%s' %(output_dir, self.plate_path[:-5].split('/')[-1])
            if not os.path.exists(dirout):
                self._log('creating dose response curve plot directory: %s' %dirout)
                os.makedirs(dirout)

            if not os.path.exists(dirout + '/dose-response-plots/' + inhibitor): os.makedirs(dirout + '/dose-response-plots/' + inhibitor)
            fig.savefig('%s/dose-response-plots/%s/dose-response-curve%s.PNG' %(dirout, inhibitor, suffix))
            plt.close(fig)
        except Exception as e:
            self._log('failed to save plot [%s] \n\tError:  %s' %(inhibitor, str(e)))

    def write_log(self, output_dir='../output'):
        '''

        '''
        try:         
            dirout = './%s/%s' %(output_dir, self.plate_path[:-5].split('/')[-1])
            with open(dirout + '/output.logs', 'w') as f:
                f.write(self.msg_log)
        except FileNotFoundError: 
            if not os.path.isdir('./' + output_dir):
                os.mkdir(output_dir)
                os.mkdir(dirout)
                with open(dirout + '/output.logs', 'w') as f:
                    f.write(self.msg_log)
                

    def write_data_to_file(self, output_dir='../output'):
        '''

        '''
        dirout = './%s/%s' %(output_dir, self.plate_path[:-5].split('/')[-1])
        self._log('writing data to file: %s/HNSCC_processed_functional_data.csv' %dirout)
        self.data.to_csv('%s/HNSCC_processed_functional_data.csv' %dirout)

    def post_processing_set_flags(self, aic_lim = 12, deviance_lim = 2):
        '''

        '''
        self.data = self.data.assign(AIC_flag = lambda x: x.prob_AIC > aic_lim,
                                     DEV_flag = lambda x: x.prob_deviance > deviance_lim,
                                     overfit_flag = lambda x: x.poly_AIC > x. prob_AIC)

    def predict_hermetric_transition(self, model_dir='../../atypical_doseresponse_classifier/best_model/best_model.159-55.24.h5'):
        '''
        This requires dependencies, the github repo can be found here: https://github.com/nathanieljevans/atypical_doseresponse_classifier.git

        This isn't functional yet... maybe it should come in later honestly in the pipeline.
        '''
        if os.path.exists(model_path):
            try:
                model = load_model(model_path)
                pred_data = self.data[['inhibitor', 'conc_norm', 'cell_viab']].pivot(index='inhibitor', columns='conc_norm', values='cell_viab')
                print(pred_data.head())

            except:
                self._log('hermetic transition prediction failed')
                self.data = self.data.assign(predicted_hermetic_transition = 'NA')
        else:
            #print('model path does not exist')
            self._log('Invalid model_path, github repo can be found here: https://github.com/nathanieljevans/atypical_doseresponse_classifier.git')
            self.data = self.data.assign(predicted_hermetic_transition = 'NA')

    def assign_panel_id(self):
        '''
        assay id's are recorded in the file: .assay_ids

        '''
        dlim = ' : '

        with open('../.assay_ids', 'r') as f:
            assay_ids = {name:id for name, id in [x for x in [x.rstrip().split(dlim) for x in f.readlines()] if len(x) > 1]}
            #print(assay_ids)

        if self.plate_path in assay_ids:
            assay_id = assay_ids[self.plate_path]
            self._log(f'This panel has been processed previously, using previously assigned panel_id [{assay_id}]')
            #print(f'This panel has been processed previously, using previously assigned panel_id [{assay_id}]')

        elif (len(assay_ids.keys()) > 0):
            nxt = max([int(x) for x in assay_ids.values()]) + 1
            with open('../.assay_ids', 'a') as f:
                f.write(f'{self.plate_path}{dlim}{nxt}\n')
            assay_id = nxt
            self._log(f'This panel has not been processed previously, assigning new panel_id [{assay_id}]')
            #print(f'This panel has not been processed previously, assigning new panel_id [{assay_id}]')
        else:
            #print('first assay id entry [1]')
            assay_id = 1
            with open('../.assay_ids', 'a') as f:
                f.write(f'{self.plate_path}{dlim}{assay_id}\n')

        self.data = self.data.assign(panel_id =  assay_id)

def process(plate_path, platemap_dir = '../plate_maps/', verbose=False, do_raise=False):
    '''
    This method loads a file representing a single plate assay into memory, maps plate locations to
    concentration, inhibitor and then preforms the processing defined by the data pipeline.

    inputs
        plate_path <str> file path to plate, must be a .xlsx file in proper naming format.
        platemap_dir <str>
    '''
    try:
        print('\t\tinitializing panel...', end='\t\t\t\t\t\r')
        p = panel(plate_path=plate_path, platemap_dir = platemap_dir, verbose=verbose)
        print('\t\tmapping data...', end='\t\t\t\t\t\r')
        p.map_data()
        print('\t\tassigning assay identifier...', end='\t\t\t\t\t\r')
        p.assign_panel_id()
        print('\t\tnormalizing combination agent concentrations...', end='\t\t\t\t\t\r')
        p.normalize_combinationagent_concentrations()
        print('\t\tnormalizing cell viability by negative controls...', end='\t\t\t\t\t\r')
        p.normalize_cell_viability_by_negative_controls()
        print('\t\tsetting floor of zero...', end='\t\t\t\t\t\r')
        p.set_floor()
        print('\t\tsetting ceiling of one...', end='\t\t\t\t\t\r')
        p.set_ceiling()
        print('\t\taveraging within plate replicates...', end='\t\t\t\t\t\r')
        p.avg_plate_replicates(method='within', flag_threshold=1)
        print('\t\taveraging across plate replicates...', end='\t\t\t\t\t\r')
        p.avg_plate_replicates(method='across', flag_threshold=0.75)
        print('\t\tfitting dose response curve...', end='\t\t\t\t\t\r')
        p.fit_regressions(plot=False)
        print('\t\tsetting post processing flags...', end='\t\t\t\t\t\r')
        p.post_processing_set_flags(aic_lim = 12, deviance_lim = 2)
        print('\t\twriting data to file...', end='\t\t\t\t\t\r')
        p.write_data_to_file()
        print('\t\twriting logs to file...', end='\t\t\t\t\t\r')
        p.write_log()
        print('\t\tcomplete.', end='\t\t\t\t\t\n')

    except Exception as e:
        p._log('Processing Failed: \n\t%s' %str(e))
        p.write_log()
        if do_raise: raise

if __name__ == '__main__':
    pd.options.display.width = 0

    # testing
    plate_path = '../data/lab_id=10139-norm=Blank490-plate_version_id=OHSU_HNSCC_derm002-note=NA.xlsx' #
    #'../data/lab_id=10004-norm=Blank490-plate_version_id=OHSU_HNSCC_derm002-note=NA.xlsx'
    platemap_dir = '../plate_maps/'

    process(plate_path, platemap_dir, verbose=False, do_raise=True)
