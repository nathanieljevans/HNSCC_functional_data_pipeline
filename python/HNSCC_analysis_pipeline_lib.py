'''
@author: Nate Evans
@date: 7/16/2019
@title: HNSCC Functional Data Analysis Pipeline

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

pd.options.display.width = 0

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

def get_plate_map(map_path, verbose=True):
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

        assert len(lab_id) == 5, 'lab_id is expected to be 5 characters long'
        self.lab_id = lab_id

        assert version_id in ['OHSU_HNSCC_derm001', 'OHSU_HNSCC_derm002'], 'unknown platemap ID, please double check pathname or update script with acceptable plate IDs'
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
                if not warned: self._log('WARNING: This assay [lab_id=%s,notes=%s] does not have a "blank490" column (last col), please double check that the data has been normalized by the positive controls.' %(lab_id,notes))
                warned = True
                if blank490: self._log('WARNING: file name specifies blank490, but data does not appeaer to be postive control normalized.')
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
        self.msg_log += str(dt.now()) + ': ' + msg + '\n'

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
        print(self.data.head(25))

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

        #print(self.data[['lab_id', 'optical_density', 'cell_viab', 'PAC']].head())

    def set_floor(self):
        '''
        There is still some contention here, Dan's protocol shifts the PAC up by the smallest value on plate.
        [and therefore this step would need to come before negative control normalization]

        I don't think that's a good idea though, rather, I'm going to round any negatives up to 0. These will be marked with 'is_adj'
        FOR NOW - talk to the group
        '''

        self.data = self.data.assign(is_adj = self.data.cell_viab < 0, cell_viab = [0 if cv < 0 else cv for cv in self.data.cell_viab] )

        print(self.data[['lab_id', 'optical_density', 'cell_viab', 'PAC']].head(25))

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
                        replicates.append(inhib)

                        new_obs = assay.groupby(['lab_id','inhibitor','conc_norm', 'plate_num'])['cell_viab'].mean().to_frame().reset_index(level=['lab_id', 'inhibitor', 'conc_norm', 'plate_num'])

                        avgd_obs.append(new_obs)

                        aucs = []
                        # groupby row - TODO: should double check that no replicates are on same row
                        for row in set(assay.plate_row):
                            repl = assay[assay['plate_row']==row]

                            aucs.append(self._get_lin_auc(repl))

                        if max(aucs) - min(aucs) > flag_threshold: toflag.append(inhib)

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
                self.data = self.data.merge(flags, how='left', by='inhibitor').assign(within_plate_repl_flag = lambda x: [True if f else False for f in x.withinplate_repl_flag])
            else:
                self.data.assign(within_plate_repl_flag = False)

        # ------------- Across plate replicate averaging VVV ------------------

        elif method == 'across':
            self._log('averaging across plate replicates...')
            avgd_obs = []

            for inhib in set(fitdat['inhibitor'].values):
                assay = fitdat[fitdat['inhibitor'] == inhib]

                #print('assay shape: %s' %str(assay.shape))
                if assay.shape[0] > 7:
                    n+=1
                    replicates.append(inhib)

                    new_obs = assay.groupby(['lab_id','inhibitor','conc_norm'])['cell_viab'].mean().to_frame().reset_index(level=['lab_id', 'inhibitor', 'conc_norm'])

                    avgd_obs.append(new_obs)

                    aucs = []
                    # groupby row - TODO: should double check that no replicates are on same row
                    for plate in set(assay.plate_num):
                        repl = assay[assay['plate_num']==plate]

                        aucs.append(self._get_lin_auc(repl))

                    self._log('[%s] linear AUCs: %r -- range (%.2f, %.2f)' %(inhib, aucs, min(aucs), max(aucs)))
                    if max(aucs) - min(aucs) > flag_threshold: toflag.append(inhib)

            self._log('There were %d across plate replicates [%r]' %(n, replicates))
            self._log('across plate replicates flagged: %r' %toflag)

            if len(replicates) > 0:
                # set is_repl flag, signifies that it has been averaged
                repls = pd.DataFrame({'inhibitor': replicates, 'is_across_plate_repl': [True]*len(replicates)})
                self.data = self.data.merge(repls, how='left', on='inhibitor').assign(is_across_plate_repl = lambda x: [True if f else False for f in x.is_across_plate_repl])
                # then add the averaged replicate observation to the bottom of dataframe
                for nobs in avgd_obs:
                    self.data = self.data.append(nobs, ignore_index=True)
            else:
                self.data = self.data.assign(is_across_plate_repl = False)

            if (len(toflag) > 0):
                flags = pd.DataFrame({'inhibitor': toflag, 'across_plate_repl_flag': [True]*len(toflag)})
                self.data = self.data.merge(flags, how='left', on='inhibitor').assign(across_plate_repl_flag = lambda x: [True if f else False for f in x.withinplate_repl_flag])
            else:
                self.data.assign(across_plate_repl_flag = False)

        else:
            raise 'choose a proper averaging method [within, across] plates'


    def _get_lin_auc(self, df, plot=False):
        '''
        fit a linear regression to the data and calculate auc

        inputs
            df <dataframe> 7 points to plot
            plot <bool> option to plot the fit; for diagnostics

        outputs
            auc <float> area under the curve
        '''
        assert df.shape[0] == 7, 'There should only be 7 observations per replicate'

        #print()
        #print(df[['inhibitor','conc_norm', 'cell_viab', 'plate_num', 'plate_row', 'plate_col']])

        x = [np.log10(x) for x in df['conc_norm'].values]
        y = df['cell_viab'].values

        pr = sm.GLM(y, sm.add_constant(x))
        glm_res = pr.fit()

        # AUC calculation -----------------------------------------------------
        # left rectangle auc estimate
        delta = 0.001
        x2 = np.arange(np.log10(min(df['conc_norm'].values)), np.log10(max(df['conc_norm'].values)), delta)
        yhat = glm_res.predict(sm.add_constant(x2))
        auc = np.sum(yhat*delta)

        if plot:
            plt.figure()
            plt.plot(x, y, 'ro', label='replicate')
            plt.plot(x2, yhat, 'g-', label='fit')
            plt.legend()
            plt.show()

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
        res = {x:[] for x in ['lab_id', 'inhibitor','intercept', 'beta1', 'auc', 'prob_conv']}
        i = 0

        # Filter controls
        pat_dat = self.data[~self.data['inhibitor'].isin(['NONE', 'F-S-V', 'DMSO'])]

        # check
        #print('----------------------------------------------------')
        #who = pat_dat[~pat_dat['is_across_plate_repl'].isin([True, False])]
        #print(who)
        #print('----------------------------------------------------')

        # Filter replicates that have been averaged into new observation
        pat_dat = pat_dat[pat_dat['is_across_plate_repl'] != True]
        pat_dat = pat_dat[pat_dat['is_within_plate_repl'] != True]

        self._log('number of assays to fit: %d' %pat_dat.shape[0])

        for inhib in set(pat_dat['inhibitor'].values):
            i+=1

            df = pat_dat[pat_dat['inhibitor'] == inhib]

            print()
            print('shape: %s' %str(df.shape))
            print('inhib: %s' %inhib)
            assert df.shape[0] == 7, 'should have exactly 7 doses! [has %d]' %df.shape[0]

            try:
                x = np.log10( df['conc_norm'].values )
                y = df['cell_viab'].values

                pr = sm.GLM(y, sm.add_constant(x), family=sm.families.Binomial(link=sm.families.links.probit()))
                glm_res = pr.fit(disp=False)

                # AUC calculation -----------------------------------------------------
                # left rectangle auc estimate
                delta = 0.001
                x2 = np.arange(np.log10(min(df['conc_norm'].values)), np.log10(max(df['conc_norm'].values)), delta)
                yhat = glm_res.predict(sm.add_constant(x2))
                auc = np.sum(yhat*delta)

                if (plot):
                    plt.subplots(1,1, figsize=(10,10))
                    plt.title('inhib: %s [AUC= %.2f]' %(inhib, auc))
                    plt.plot(x2, yhat, 'r-', label='probit_fit')
                    plt.plot(x, y, 'bo', label='replicates')
                    plt.legend()
                    plt.show()

                # beta0 = intercept
                # beta1 = slope
                (intercept,beta1) = glm_res.params

                # update results
                [res[var].append(val) for var,val in zip(['inhibitor','intercept', 'beta1', 'auc', 'prob_conv'], [inhib, intercept ,beta1 , auc, True])]

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
                    auc = self._get_lin_auc(df, plot=plot)

                [res[var].append(val) for var,val in zip(['inhibitor','intercept', 'beta1', 'auc', 'prob_conv'], [inhib, None, None, auc, False])]

            except:
                print('some other exception...')
                failures.append( inhib )
                [res[var].append(val) for var,val in zip(['inhibitor','intercept', 'beta1', 'auc'], [inhib, 'NA' ,'NA' , 'NA'])]
                if plot:
                    print('----------------------------------------------------')
                    print('FAILURE: %s' %inhib)
                    print(df.head(7))
                    print('----------------------------------------------------')


                    f, ax = plt.subplots(1,1, figsize = (10,10))
                    ax.set_title('FAILURE: %s' %(inhib))
                    #sbn.scatterplot(x=x2, y=yhat , ax=ax)
                    plt.xscale('log')
                    sbn.scatterplot(x=x, y=y, ax=ax)
                    plt.show()
                raise

        print('Failures [%d]: %r' %(len(failures),failures))


if __name__ == '__main__':

    # testing
    plate_path = '../data/lab_id=10004-norm=Blank490-plate_version_id=OHSU_HNSCC_derm002-note=NA.xlsx'
    platemap_dir = '../plate_maps/'

    p = panel(plate_path=plate_path, platemap_dir = platemap_dir, verbose=True)
    p.map_data()
    p.normalize_combinationagent_concentrations()
    p.normalize_cell_viability_by_negative_controls()
    p.set_floor()
    p.avg_plate_replicates(method='within', flag_threshold=1)
    p.avg_plate_replicates(method='across', flag_threshold=0.75)
    p.set_ceiling()
    p.fit_regressions()


    #print(p.msg_log)
