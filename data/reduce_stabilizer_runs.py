from utils import *


def make_internal_standard_correction(data, var='CPS'):
    def get_int_std_data(data):
        sample_list = data.complete_sample_name.unique()
        element_description_list = data.element_description.unique()
        
        ref_run_number = data.loc[data['Sample Name'] == 'blank-1'].iloc[0].run_number
        ref_data = data[(data.run_number==ref_run_number)][['element_description', var]].values
        ref_cps_dict = {k:v for k,v in ref_data}

        data_int_std = []
        for sample in sample_list:
            data_sample = data[data.complete_sample_name == sample]
            for ed in element_description_list:
                if 'I' in ed:
                    int_std = '128  Te  [ No Gas solution ]'
                else:
                    int_std = '103  Rh  [ No Gas solution ]'

                ref_cps_int_std = ref_cps_dict[int_std] 
                cps_int_std = data_sample[(data_sample.element_description == int_std)][var].values[0]
                data_int_std.append([sample, ed, ref_cps_int_std, cps_int_std])
        columns = ['complete_sample_name', 'element_description', f'int_std_ref_{var}', 
            f'int_std_{var}']
        data_int_std = pd.DataFrame(data_int_std, columns=columns)
        return data_int_std

    # remove samples without internal std
    no_intstd_samples = ['JLs-1-test', 'A1-test']
    data = data[~data['Sample Name'].isin(no_intstd_samples)]

    data_int_std = get_int_std_data(data)
    data = pd.merge(data, data_int_std, on=['complete_sample_name', 'element_description'])


    data[f'{var}_intstdc'] = data[var] * (data[f'int_std_ref_{var}']/data[f'int_std_{var}'])
    return data

def make_calibration_curves(data, run_num, var, only_Ca=False):
    # read standards
    standards = read_standards202405(run=run_num, only_Ca=only_Ca)

    data = data[data.element.isin(standards.element.unique())]

    # make calibration curves
    data_std = data[data['Sample Name'].isin(standards.Name.unique())]
    standards = standards.rename(columns={'Name': 'Sample Name'})
    data_std = pd.merge(data_std, standards, on=['Sample Name', 'element'], 
        suffixes=('','_known'))
      
    def linear_fit_zero_interc_slope(element_description):
        data_std_el = data_std[data_std.element_description == element_description]
        if len(data_std_el['know_mol/mol']) >= 2:
            a, _, _, _ = np.linalg.lstsq(data_std_el['know_mol/mol'].to_numpy()[:,np.newaxis],
                data_std_el[var])
            return a[0]
        else:
            return np.nan

    data['zero_interc_slope'] = data.element_description.apply(
        linear_fit_zero_interc_slope)

    data['inferred_mol/mol'] =\
        data[var]/data['zero_interc_slope']
    return data

def normalize_to_Ca_Mg_iodine_runs(data, var, only_Ca=False):
    mg_description = '25  Mg  [ No Gas solution ]'
    ca_description = '44  Ca  [ No Gas solution ]'
    data_Ca = data[data.element_description == ca_description]
    data_Mg = data[data.element_description == mg_description]
    
    for el in ['Ca', 'Mg']:
        data_el = data_Ca if el == 'Ca' else data_Mg
        data_el = data_el[[var, 'complete_sample_name']]

        suffix = '_44_Ca' if el == 'Ca' else '_25_Mg'
        data = pd.merge(data, data_el, on='complete_sample_name',
            how='left', suffixes=('', suffix))
    
    # normalize to Ca + Mg
    if only_Ca:
        norm = data[f'{var}_44_Ca'] + data[f'{var}_25_Mg']
    else:
        norm = data[f'{var}_44_Ca'] + data[f'{var}_25_Mg']
    data[f'{var}_CaMgnorm'] = data[var] / norm
    return data

def make_calibration_curves_ppb(data, run_num, var='CPS_Te_normalized'):
    # read standards
    standards = read_standards202405_ppb(run=run_num)

    data = data[data.element.isin(standards.element.unique())]

    # make calibration curves
    data_std = data[data['Sample Name'].isin(standards.Name.unique())]
    standards = standards.rename(columns={'Name': 'Sample Name'})
    data_std = pd.merge(data_std, standards, on=['Sample Name', 'element'], 
        suffixes=('','_known'))
      
    def linear_fit_zero_interc_slope(element_description):
        data_std_el = data_std[data_std.element_description == element_description]
        if len(data_std_el['ppb']) >= 2:
            a, _, _, _ = np.linalg.lstsq(data_std_el['ppb'].to_numpy()[:,np.newaxis],
                data_std_el[var])
            return a[0]
        else:
            return np.nan

    data['zero_interc_slope'] = data.element_description.apply(
        linear_fit_zero_interc_slope)

    data['inferred_ppb_zero_intercept'] =\
        data[var]/data['zero_interc_slope']
    
     # apply dilution factor to cpnvert from ppb in solutions to ppb in samples
    df = lambda x : x * 4000
    data['inferred_ppb_zero_intercept'] = \
        data['inferred_ppb_zero_intercept'].apply(df)
    return data

def make_drift_correction(data, run_num=8, var='CPS'):
     # find drift solution samples for run
    drift_solution = 'JLs-1 May 30, 2024'
    if run_num == 10:
        drift_solution = 'JLs-1 June 5, 2024'
    if run_num == 11:
        drift_solution = 'JLs-1 June 6, 2024'

    blank1_run_num = data.loc[data['Sample Name'] == 'blank-1'].iloc[0].run_number

    drift_evolution = data[data['Sample Name'] == drift_solution]
    drift_evolution = drift_evolution[drift_evolution['run_number'] > blank1_run_num]
    drift_evolution = data[data['Sample Name'] == drift_solution]

    first_drift_solution = drift_evolution.run_number.min()
    ref_drift_solution = drift_evolution[drift_evolution.run_number==first_drift_solution]
    ref_drift_solution = ref_drift_solution[['element_description', var]]

    drift_evolution = pd.merge(drift_evolution, ref_drift_solution,
        on=['element_description'], suffixes=('', '_drift_ref'))

    drift_evolution['drift_correction_factor'] = drift_evolution[var]/ drift_evolution[var+'_drift_ref']

    # make dict of correction curves
    drift_polynomials = {}
    for ed in drift_evolution.element_description.unique():
        drift_evolution_el = drift_evolution[drift_evolution.element_description == ed]
        p = np.polyfit(drift_evolution_el.run_number,
            drift_evolution_el.drift_correction_factor, 3)
        drift_polynomials[ed] = p

    def get_drift_factors(row):
        return np.poly1d(drift_polynomials[row['element_description']])(row['run_number'])
    

    data['drift_correction_poly_coefs'] = data.element_description.map(drift_polynomials)
    data['drift_correction_factor'] = data.apply(get_drift_factors,axis=1)
    data[f'{var}_driftc'] = data[var]/ data.drift_correction_factor
    return data


def merge_data(ppb=False, only_Ca=False):
    run8 = read_and_format_data(run_number=8)
    run9 = read_and_format_data(run_number=9)
    run10 = read_and_format_data(run_number=10)
    run11 = read_and_format_data(run_number=11)

    def add_iodine_cutoff_values(run, run_num):
        # find drift solution samples for run
        drift_solution = 'JLs-1 May 30, 2024'
        if run_num == 10:
            drift_solution = 'JLs-1 June 5, 2024'
        if run_num == 11:
            drift_solution = 'JLs-1 June 6, 2024'
        
        # if drift_correction:
            # run = make_drift_correction(run)
        blank1_run_num = run.loc[run['Sample Name'] == 'blank-1'].iloc[0].run_number
        run = run[run['run_number'] > blank1_run_num]
        
        drift = run[run['Sample Name'] == drift_solution]
        # drift = drift[drift['run_number'] > blank1_run_num]
        # find sigma of iodine concentration
        drift = drift[drift['element'] == 'I']
        conc = 'CPS_driftc'
        std = drift[conc].std()
        # cutoff = max_sample_blank + 2 sigma
        max_sample_blank = run[(run.is_sample_blank)&(run.element=='I')][conc].max()
        cutoff = max_sample_blank + 3*std

        # run = run [run.element == 'I']
        run['iodine_cutoff'] = cutoff
        iodine_data = run[run.element=='I'][['complete_sample_name', conc]]
        iodine_data = iodine_data.rename(columns={conc: 'iodine_cps'})
        run = pd.merge(run, iodine_data, on='complete_sample_name', how='left')
        run['has_iodine'] = run['iodine_cps'] > run['iodine_cutoff']
        return run

    def add_precision(run, run_num):
        drift_solution = 'JLs-1 May 30, 2024'
        if run_num == 10:
            drift_solution = 'JLs-1 June 5, 2024'
        if run_num == 11:
            drift_solution = 'JLs-1 June 6, 2024'
        drift = run[run['Sample Name']==drift_solution]
        drift = drift[drift.element=='I']
        std_raw = drift['CPS'].std()
        mean_raw = drift['CPS'].mean()
        std_corr = drift['CPS_driftc_intstdc_sample_blank_corrected'].std()
        mean_corr = drift['CPS_driftc_intstdc_sample_blank_corrected'].mean()
        run['drift_precision_corr'] = std_corr/mean_corr
        run['drift_precision_raw'] = std_raw/mean_raw
        return run

    def find_concentrations(run, run_number):
        run = read_and_format_data(run_number=run_number)
        run = make_drift_correction(run, run_number, var='CPS')
        run = make_internal_standard_correction(run, var='CPS_driftc')
        run = blank_correction(run, var='CPS_driftc_intstdc')
        run = normalize_to_Ca_Mg_iodine_runs(run, var='CPS_driftc_intstdc_sample_blank_corrected', only_Ca=only_Ca)
        run = make_calibration_curves(run, run_number, var='CPS_driftc_intstdc_sample_blank_corrected_CaMgnorm', only_Ca=only_Ca)
        run = add_iodine_cutoff_values(run, run_number)
        run = add_precision(run, run_number)
        run['run_id'] = str(run_number)
        return run
    
    def find_concentrations_ppb(run, run_number):
        run = read_and_format_data(run_number=run_number)
        run = make_drift_correction(run, run_number, var='CPS')
        run = make_internal_standard_correction(run, var='CPS_driftc')
        run = blank_correction(run, var='CPS_driftc_intstdc')
        run = make_calibration_curves_ppb(run, run_number, 
            var='CPS_driftc_intstdc_sample_blank_corrected')
        run = add_iodine_cutoff_values(run, run_number)
        run['run_id'] = str(run_number)
        return run
            

    runs_list = [run8, run9, run10, run11]
    run_numbers = [8,9,10,11]

    if ppb:
        runs = [find_concentrations_ppb(run, run_numbers[i]) for i, run in 
            enumerate(runs_list)]
    else:
        runs = [find_concentrations(run, run_numbers[i]) for i, run in 
            enumerate(runs_list)]
    runs = pd.concat(runs)
    
    zero_conc = lambda x : 0 if x <0 else x
    if ppb:
        runs['inferred_ppb_zero_intercept'] = \
            runs['inferred_ppb_zero_intercept'].apply(zero_conc)
        # drop duplicated and keep minimum concentration
        runs = runs.sort_values(['inferred_ppb_zero_intercept']).drop_duplicates(
            subset=['Sample Name', 'element_description'], keep='first')
    else:
        runs['inferred_mol/mol'] = \
            runs['inferred_mol/mol'].apply(zero_conc)
        # drop duplicated and keep minimum concentration
        runs = runs.sort_values(['inferred_mol/mol']).drop_duplicates(
            subset=['Sample Name', 'element_description'], keep='first')

    runs = runs.reset_index(drop=True)
    return runs


data = merge_data()
data_Ca_norm = merge_data(only_Ca=True)
data_ppb = merge_data(ppb=True)
data.to_csv('stabilizer_runs_reduced.csv', index=0)
data_Ca_norm.to_csv('stabilizer_runs_reduced_Ca_norm.csv', index=0)
data_ppb.to_csv('stabilizer_runs_reduced_ppb.csv', index=0)