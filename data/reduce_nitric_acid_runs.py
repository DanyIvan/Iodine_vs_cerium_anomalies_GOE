from utils import *

pd.options.mode.chained_assignment = None 


def merge_icpms_data():
    run4 = read_and_format_data(run_number=4)
    run5 = read_and_format_data(run_number=5)
    run6 = read_and_format_data(run_number=6)

    def add_iodine_cutoff_values(run, outlier_quantile=0.85):
        # find drift solution samples for run
        drift_solution = 'DG JLs-1 Drift' 
        drift = run[run['Sample Name'] == drift_solution]
        # find sigma of iodine concentration
        drift = drift[drift['element'] == 'I']
        ql = drift.CPS.quantile(outlier_quantile)
        drift = drift[drift.CPS <= ql]
        # conc = 'inferred_mol/mol_zero_intercept'
        conc = 'CPS'
        std = drift[conc].std()
        # cutoff = max_sample_blank + 2 sigma
        max_sample_blank = run[(run.is_sample_blank)&(run.element=='I')][conc].max()
        cutoff = max_sample_blank + 2*std
        run['iodine_cutoff'] = cutoff
        iodine_data = run[run.element=='I'][['Sample Name', 'run_number', conc]]
        iodine_data = iodine_data.rename(columns={conc: 'iodine_cps'})
        run = pd.merge(run, iodine_data, on=['Sample Name', 'run_number'],
            how='left')
        run['has_iodine'] = run['iodine_cps'] > run['iodine_cutoff']
        return run

    def find_concentrations(run, run_number, run7_Mg_fill_values=None):
        run = blank_correction(run, var='CPS')
        if run_number != 7:
            run = make_oxide_correction(run, var='CPS_sample_blank_corrected')
        run = normalize_to_Ca(run, run_num=run_number)
        run = make_calibration_curves_Ca(run,
            var='CPS_Ca_normalized_sample')
        run = add_iodine_cutoff_values(run)
        run['run_id'] = str(run_number)
        
        # normalize to Ca + Mg
        run = run[run.is_sample]
        run = run[~run.is_sample_blank]
        
        if run_number == 7:
            run = normalize_to_Ca_Mg(run, run_num=run_number,
                run7_Mg_fill_values=run7_Mg_fill_values)
        else:
            run = normalize_to_Ca_Mg(run, run_num=run_number)
        return run

    runs_list = [run4, run5, run6]
    run_numbers = [4,5,6]

    runs = [find_concentrations(run, run_numbers[i]) for i, run in enumerate(runs_list)]
    runs = pd.concat(runs)


    zero_conc = lambda x : 0 if x <0 else x
    runs['inferred_mol/mol'] = \
        runs['inferred_mol/mol'].apply(zero_conc)

    runs = runs.reset_index(drop=True)

    # save csv
    # filter columns of interest
    runs = runs[['is_sample', 'Sample Name', 'element',
        'element_description', 'complete_sample_name', 'is_blank', 'is_sample_blank',
        'is_standard', 'has_iodine', 'inferred_mol/mol', 'run_id']]
    runs.to_csv('nitric_acid_runs_reduced.csv')

    
def merge_icpms_data_ppb(jls=False):
    run4 = read_and_format_data(run_number=4)
    run5 = read_and_format_data(run_number=5)
    run6 = read_and_format_data(run_number=6)
    run7 = read_and_format_data(run_number=7)

    def find_concentrations(run, run_number):
        run = blank_correction(run, var='CPS')
        if run_number != 7:
            run = make_oxide_correction(run, var='CPS_sample_blank_corrected')
        run = normalize_to_In(run, var='CPS_sample_blank_corrected')
        make_calibration_curves(run, var='CPS_In_normalized_sample')
        run['run_id'] = str(run_number)
        return run

    runs = [run4, run5, run6, run7]
    run_numbers = [4,5,6,7]

    runs = [find_concentrations(run, run_numbers[i]) for i, run in enumerate(runs)]
    runs = pd.concat(runs)
    runs = runs[runs.is_sample]
    runs = runs[~runs.is_sample_blank]

    zero_conc = lambda x : 0 if x <0 else x
    runs['inferred_ppb_zero_intercept'] = \
        runs['inferred_ppb_zero_intercept'].apply(zero_conc)
    
    # apply dilution factor to cpnvert from ppb in solutions to ppb in samples
    df = lambda x : x * 4000
    runs['inferred_ppb_zero_intercept'] = \
        runs['inferred_ppb_zero_intercept'].apply(df)

    runs = runs.reset_index(drop=True)

    # drop duplicates
    # TODO: these duplictaes might be the result of mislabeling. Check vials
    # to make sure

    duplicates = runs[['Sample Name', 'element_description', 'run_id']]
    duplicates = duplicates[duplicates[
        ['Sample Name', 'element_description']].duplicated()]
    duplicate_names = duplicates['Sample Name'].unique()
    duplicates = runs[runs['Sample Name'].isin(duplicate_names)]
    duplicates_dict = {}
    for name in duplicate_names:
        duplicates_dict[name] =\
            duplicates[duplicates['Sample Name']==name].run_id.unique()

    for duplicate in duplicates_dict.keys():
        # print(duplicate)
        delete_idxs = []
        if '7' in duplicates_dict[duplicate]:
            iodine_delete = runs[(
                runs['Sample Name'] == duplicate) &
                (runs.element == 'I') &
                (runs.run_id != '7')].index
            mg_delete = runs[(
                runs['Sample Name'] == duplicate) &
                (runs.element == 'Mg') & 
                (runs.run_id == '7')].index
            delete_idxs += iodine_delete.to_list()
            delete_idxs += mg_delete.to_list()
        if len(duplicates_dict[duplicate]) > 1:
            first_run = duplicates_dict[duplicate][0]
            idxs = runs[(
                runs['Sample Name'] == duplicate) &
                (runs.run_id != first_run) &
                (runs.run_id != '7')].index
            delete_idxs += idxs.to_list()
        delete_idxs = set(delete_idxs)
        # print(delete_idxs)
        runs.drop(delete_idxs, inplace=True)

    # change data to long format
    runs = runs.pivot_table(values='inferred_ppb_zero_intercept',
        columns='element_description', index=['Sample Name'])
    
    runs = runs.reset_index()
    change_names = lambda x: x.replace('18-2-', '')
    runs['Sample Name'] = runs['Sample Name'].apply(change_names)

    def assign_formation(sample_name):
        if sample_name.startswith('SR'):
            return 'Steep Rock'
        elif sample_name.startswith('B'):
            return 'Gordon Lake'
        else:
            return 'Espanola'
        
    # assign formation names
    runs['formation'] = runs['Sample Name'].apply(assign_formation)
    runs.to_csv('nitric_acid_runs_reduced_ppb.csv')
    # return runs


def normalize_to_In(data, run=1, var='CPS_HNO3_blank_corrected'):
    # normalize to In
    
    def get_In_data(data):
        sample_list = data.complete_sample_name.unique()
        element_description_list = data.element_description.unique()

        ref_run_number = 16 if run==3 else 1
        In_list = ['115  In  [ No Gas ]', '115 -> 115  In  [ H2 ]', '115 -> 115  In  [ O2 ]']
        In_ref_cps = data[(data.run_number==ref_run_number)&(
            data.element_description.isin(In_list))][['element_description', 'CPS']].values
        In_ref_cps_dict = {k:v for k,v in In_ref_cps}

        data_In = []
        for sample in sample_list:
            data_sample = data[data.complete_sample_name == sample]
            for ed in element_description_list:
                if 'No Gas' in ed:
                    In_ed = '115  In  [ No Gas ]'
                elif 'H2' in ed:
                    In_ed = '115 -> 115  In  [ H2 ]'
                else:
                    In_ed = '115 -> 115  In  [ O2 ]'

                ref_cps_In = In_ref_cps_dict[In_ed] 
                cps_In = data_sample[(data_sample.element_description == In_ed)]['CPS'].values[0]
                data_In.append([sample, ed, ref_cps_In, cps_In])
        columns = ['complete_sample_name', 'element_description', 'In_ref_CPS', 
            'In_CPS']
        data_In = pd.DataFrame(data_In, columns=columns)
        return data_In
    
    if not 'In_ref_CPS' in data:
        data_In = get_In_data(data)
        data = pd.merge(data, data_In, on=['complete_sample_name', 'element_description'])

    if 'HNO3' in var:
        id_str = '_HNO3'
    elif 'sample' in var:
        id_str = '_sample'
    else:
        id_str = ''

    data[f'CPS_In_normalized{id_str}'] = data[var] * (data['In_ref_CPS']/data['In_CPS'])
    return data

def make_drift_correction(data, run=1, var='CPS_In_normalized'):
    if run==1:
        drift_solution = 'DG JLs-1 DRIFT'
    elif run == 4:
        drift_solution = 'DG JLs-1 Drift '
    else:
        drift_solution = 'JLS-1 1/1000'

    drift_evolution = data[data['Sample Name'] == drift_solution]
    if run ==1:
        drift_evolution = drift_evolution[drift_evolution.run_number_int > 8]
    first_drift_solution = drift_evolution.run_number_int.min()
    ref_drift_solution = drift_evolution[drift_evolution.run_number_int==first_drift_solution]
    ref_drift_solution = ref_drift_solution[['element_description', var]]

    drift_evolution = pd.merge(drift_evolution, ref_drift_solution,
        on=['element_description'], suffixes=('', '_drift_ref'))

    drift_evolution['drift_correction_factor'] = drift_evolution[var]/ drift_evolution[var+'_drift_ref']

    # make dict of correction curves
    drift_polynomials = {}
    for ed in drift_evolution.element_description.unique():
        drift_evolution_el = drift_evolution[drift_evolution.element_description == ed]
        p = np.polyfit(drift_evolution_el.run_number_int,
            drift_evolution_el.drift_correction_factor, 3)
        drift_polynomials[ed] = p

    def get_drift_factors(row):
        return np.poly1d(drift_polynomials[row['element_description']])(row['run_number_int'])
    

    data['drift_correction_poly_coefs'] = data.element_description.map(drift_polynomials)
    data['drift_correction_factor'] = data.apply(get_drift_factors,axis=1)
    data[f'{var}_drift_corrected'] = data[var]/ data.drift_correction_factor
    return data

    

def normalize_to_Ca(data, run_num=7):
    mg_description = '25  Mg  [ No Gas ]' if run_num > 7 else '25 -> 25  Mg  [ H2 ]'
    ca_description = '44 -> 44  Ca  [ H2 ]' if run_num ==1 else '44  Ca  [ No Gas ]'
    data_Ca = data[data.element_description == ca_description]
    data_Mg = data[data.element_description == mg_description]

    for el in ['Ca']:
        data_el = data_Ca if el == 'Ca' else data_Mg
        data_el_hno3 = data_el[['CPS_HNO3_blank_corrected', 'complete_sample_name']]
        data_el_spl = data_el[['CPS_sample_blank_corrected', 'complete_sample_name']]
        data_el_cps = data_el[['CPS', 'complete_sample_name']]

        suffix = '_44_Ca' if el == 'Ca' else '_25_Mg'
        data = pd.merge(data, data_el_hno3, on='complete_sample_name',
            how='left', suffixes=('', suffix))
        data = pd.merge(data, data_el_spl, on='complete_sample_name',
            how='left', suffixes=('', suffix))
        data = pd.merge(data, data_el_cps, on='complete_sample_name',
            how='left', suffixes=('', suffix))
    
    # normalize to Ca
    data['CPS_Ca_normalized_HNO3'] =\
        data['CPS_HNO3_blank_corrected'] / data['CPS_HNO3_blank_corrected_44_Ca']
    data['CPS_Ca_normalized_sample'] =\
        data['CPS_sample_blank_corrected'] / data['CPS_sample_blank_corrected_44_Ca']
    data['CPS_Ca_normalized'] =\
        data['CPS'] / data['CPS_44_Ca']
    
    # # normalize to Ca + Mg
    # Ca_Mg_hno3 = data['CPS_HNO3_blank_corrected_44_Ca'] +\
    #     data['CPS_HNO3_blank_corrected_25_Mg']
    # data[f'CPS_Ca_Mg_normalized_HNO3'] =\
    #     data['CPS_HNO3_blank_corrected'] / Ca_Mg_hno3
    # Ca_Mg_sample = data['CPS_sample_blank_corrected_44_Ca'] +\
    #     data['CPS_sample_blank_corrected_25_Mg']
    # data[f'CPS_Ca_Mg_normalized_sample'] =\
    #     data['CPS_sample_blank_corrected'] / Ca_Mg_sample
    # Ca_Mg = data['CPS_44_Ca'] + data['CPS_25_Mg']
    # data[f'CPS_Ca_Mg_normalized'] = data['CPS'] / Ca_Mg
    return data

def make_calibration_curves(data, var='CPS_In_normalized_HNO3', standards='new'):
    # read standards
    if standards == 'new':
        standards = read_standards202311_ppb()
    elif standards == 'old':
        standards = read_standards202309(with_Ca=False)

    # make ppm calibration curves
    data_std = data[data['Sample Name'].isin(standards.Name.unique())]
    standards = standards.rename(columns={'Name': 'Sample Name'})
    standards = standards.astype({'ppb':'float'})
    data_std = pd.merge(data_std, standards, on=['Sample Name', 'element'], 
        suffixes=('','_known'))
    data_std.ppb = data_std.ppb.fillna(0)

    def linear_fit_slope(element_description):
        # print(element_description)
        data_std_el = data_std[data_std.element_description == element_description]
        if len(data_std_el.ppb) >= 2:
            a, b = np.polyfit(data_std_el.ppb, data_std_el[var], deg=1)
            return a
        else:
            return np.nan

    def linear_fit_intercept(element_description):
        data_std_el = data_std[data_std.element_description == element_description]
        # print(data_std_el)
        if len(data_std_el.ppb) >= 2:
            a, b = np.polyfit(data_std_el.ppb, data_std_el[var], deg=1)
            return b
        else:
            return np.nan
        
    def linear_fit_zero_interc_slope(element_description):
        data_std_el = data_std[data_std.element_description == element_description]
        if len(data_std_el.ppb) >= 2:
            a, _, _, _ = np.linalg.lstsq(data_std_el.ppb.to_numpy()[:,np.newaxis],
                data_std_el[var])
            return a[0]
        else:
            return np.nan

    data['slope'] = data.element_description.apply(linear_fit_slope)
    data['intercept'] = data.element_description.apply(
        linear_fit_intercept)
    data['zero_interc_slope'] = data.element_description.apply(
        linear_fit_zero_interc_slope)

    data['inferred_ppb'] =\
        (data[var] - data['intercept'])/data['slope'] 

    data['inferred_ppb_zero_intercept'] =\
        data[var]/data['zero_interc_slope']

def make_calibration_curves_Ca(data, var='CPS_Ca_normalized',
        standards_type='new'):
    # read standards
    if standards_type == 'new':
        standards = read_standards202311()
    else:
        standards = read_standards202309(with_Ca=True)

    data = data[data.element.isin(standards.element.unique())]

    # make calibration curves
    data_std = data[data['Sample Name'].isin(standards.Name.unique())]
    standards = standards.rename(columns={'Name': 'Sample Name'})
    data_std = pd.merge(data_std, standards, on=['Sample Name', 'element'], 
        suffixes=('','_known'))

    # is difficult to fit to a line Ca normalized data for I and Ca in run1
    # so I just exclude these elements. Otherwise linear solver does not converge
    if standards_type == 'old':
        data_std = data_std[~data_std.element.isin(['I', 'Ca'])]
      
    def linear_fit_zero_interc_slope(element_description):
        data_std_el = data_std[data_std.element_description == element_description]
        if len(data_std_el['mol/mol_Ca']) >= 2:
            a, _, _, _ = np.linalg.lstsq(data_std_el['mol/mol_Ca'].to_numpy()[:,np.newaxis],
                data_std_el[var])
            return a[0]
        else:
            return np.nan

    data['zero_interc_slope'] = data.element_description.apply(
        linear_fit_zero_interc_slope)

    data['inferred_mol/mol_zero_intercept'] =\
        data[var]/data['zero_interc_slope']
    return data

def normalize_to_Ca_Mg(data, run_num=7, run7_Mg_fill_values=None):
    mg_description = '25  Mg  [ No Gas ]' if run_num == 7 else '25 -> 25  Mg  [ H2 ]'

    if run_num == 7:
        Mg_per_Ca = run7_Mg_fill_values 
    else:
        Mg_per_Ca = data[data.element_description == mg_description]\
            [['Sample Name', 'inferred_mol/mol_zero_intercept']]
        Mg_per_Ca.columns = ['Sample Name', 'Mg_Ca_mol/mol']
    
    data = pd.merge(data, Mg_per_Ca, on='Sample Name',
            how='left')
    data['inferred_mol/mol_Mg'] =\
        data['inferred_mol/mol_zero_intercept'] / data['Mg_Ca_mol/mol']
    
    a = data['inferred_mol/mol_zero_intercept']
    b = data['inferred_mol/mol_Mg']
    conc_normalized_Ca_Mg = a*b/(a+b)
    data['inferred_mol/mol'] = conc_normalized_Ca_Mg.fillna(0)
    return data
    

def make_oxide_correction(data, var='CPS_sample_blank_corrected'):
    standards = read_standards202311()
    data_std =  data[data['Sample Name'].isin(standards.Name)]

    # solution without heavy REE, so no oxide interference
    no_hree_solutions = ['DGcalA 1/10']
    data_no_hree = data_std[data_std['Sample Name'].isin(no_hree_solutions)]

    def apply_correction(data, oxide_producer, element):
        # cps coming from oxide producer (ex. Ba)
        cps_oxide_prod = data_no_hree[data_no_hree.element_description
            ==oxide_producer][var].to_numpy()
        # cps coming from element to be corrected (ex. Eu)
        cps_element = data_no_hree[data_no_hree.element_description
            ==element][var].to_numpy()
        correction_factor = cps_element/cps_oxide_prod
        print(correction_factor)
        correction_factor = np.min(correction_factor)

        # measured cps of oxide producer in sample
        sample_oxide_prod = data[data.element_description
            ==oxide_producer][var].to_numpy()
        # measured cps coming from oxide interference
        oxide_interference = correction_factor * sample_oxide_prod
        
        # make correction
        data_oxide = data[data.element_description == element].copy()
        data_oxide[var] = data_oxide[var]  - oxide_interference
        data_oxide['element_description'] = element + ' oxide_corrected'
        return pd.concat([data, data_oxide])

    data = apply_correction(data, '137  Ba  [ No Gas ]', '151  Eu  [ No Gas ]')
    data = apply_correction(data, '141  Pr  [ No Gas ]', '157  Gd  [ No Gas ]')
    return data

# merge_icpms_data()
# merge_icpms_data_ppb()