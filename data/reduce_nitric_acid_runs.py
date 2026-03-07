from utils import *

pd.options.mode.chained_assignment = None 


def merge_icpms_data():
    run4 = read_and_format_data(run_number=4)
    run5 = read_and_format_data(run_number=5)
    run6 = read_and_format_data(run_number=6)
    run7 = read_and_format_data(run_number=7)

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

    # I have to fill the run7 Mg values with values from other runs because I 
    # do not have a calibration curve for Mg in run 7
    run7_Mg_fill_values = runs[runs.element_description == '25 -> 25  Mg  [ H2 ]']\
        [['Sample Name', 'inferred_mol/mol_zero_intercept']]
    run7_Mg_fill_values.columns = ['Sample Name', 'Mg_Ca_mol/mol']
    run7 = find_concentrations(run7, 7, run7_Mg_fill_values=run7_Mg_fill_values)
    runs = pd.concat([runs, run7])

    # filter aragonite and JLS samples
    is_arag_or_jls = lambda x: x.startswith('A') or 'JLs' in x
    arag_or_jls = runs[runs['Sample Name'].apply(is_arag_or_jls)]


    zero_conc = lambda x : 0 if x <0 else x
    runs['inferred_mol/mol'] = \
        runs['inferred_mol/mol'].apply(zero_conc)

    runs = runs.reset_index(drop=True)

    # drop duplicates
    # TODO: these duplictaes might be the result of mislabeling. Check vials
    # to make sure

    def drop_duplicates(runs):
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
                # if iodine is duplicated take the measurement from run 7
                iodine_delete = runs[(
                    runs['Sample Name'] == duplicate) &
                    (runs.element == 'I') &
                    (runs.run_id != '7')].index
                # if Mg is duplicated delete the measurement from run 7
                mg_delete = runs[(
                    runs['Sample Name'] == duplicate) &
                    (runs.element == 'Mg') & 
                    (runs.run_id == '7')].index
                # if Ca is duplicated delete the measurement from run 7
                ca_delete = runs[(
                    runs['Sample Name'] == duplicate) &
                    (runs.element == 'Ca') & 
                    (runs.run_id == '7')].index
                delete_idxs += iodine_delete.to_list()
                delete_idxs += mg_delete.to_list()
                delete_idxs += ca_delete.to_list()
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
    
    drop_duplicates(runs)

    # save samples with iodine
    iodine_samples = runs[runs.has_iodine]
    iodine_samples=iodine_samples[iodine_samples.is_sample]['Sample Name'].unique()
    np.savetxt('iodine_samples', iodine_samples, fmt="%s", delimiter=',')

    # print(runs.columns)
    # change data to long format
    # runs = runs.pivot_table(values='inferred_mol/mol',
    #     columns='element_description', index=['Sample Name'])
    # print(runs.columns)
    # arag_or_jls = arag_or_jls.pivot_table(values='inferred_mol/mol',
    #     columns='element_description', index=['Sample Name'])

    # save csv
    # filter columns of interest
    runs = runs[['is_sample', 'Sample Name', 'element',
        'element_description', 'complete_sample_name', 'is_blank', 'is_sample_blank',
        'is_standard', 'has_iodine', 'inferred_mol/mol', 'run_id']]
    runs.to_csv('nitric_acid_runs_reduced.csv')
    # arag_or_jls.to_csv('arag_jls_mol_per_mol_Ca_Mg.csv')
    
def get_data_mol_mol():
    run4 = read_and_format_data(run_number=4)
    run5 = read_and_format_data(run_number=5)
    run6 = read_and_format_data(run_number=6)
    run7 = read_and_format_data(run_number=7)

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
    
    def add_precision(run):
        drift_solution = 'DG JLs-1 Drift'
        blank1_run_num = run.loc[run['Sample Name'] == 'blank-1'].iloc[0].run_number
        run = run[run['run_number'] > blank1_run_num]
        drift = run[run['Sample Name']==drift_solution]
        drift = drift[drift.element=='I']
        std_raw = drift['CPS'].std()
        mean_raw = drift['CPS'].mean()
        std_corr = drift['inferred_mol/mol_zero_intercept'].std()
        mean_corr = drift['inferred_mol/mol_zero_intercept'].mean()
        run['drift_precision_corr'] = std_corr/mean_corr
        run['drift_precision_raw'] = std_raw/mean_raw
        return run

    def find_concentrations(run, run_number, run7_Mg_fill_values=None):
        run = blank_correction(run, var='CPS')
        if run_number != 7:
            run = make_oxide_correction(run, var='CPS_sample_blank_corrected')
        run = normalize_to_Ca(run, run_num=run_number)
        run = make_calibration_curves_Ca(run,
            var='CPS_Ca_normalized_sample')
        run = add_iodine_cutoff_values(run)
        run = add_precision(run)
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

    # I have to fill the run7 Mg values with values from other runs because I 
    # do not have a calibration curve for Mg in run 7
    run7_Mg_fill_values = runs[runs.element_description == '25 -> 25  Mg  [ H2 ]']\
        [['Sample Name', 'inferred_mol/mol_zero_intercept']]
    run7_Mg_fill_values.columns = ['Sample Name', 'Mg_Ca_mol/mol']
    run7 = find_concentrations(run7, 7, run7_Mg_fill_values=run7_Mg_fill_values)
    runs = pd.concat([runs, run7])

    # filter aragonite and JLS samples
    is_arag_or_jls = lambda x: x.startswith('A') or 'JLs' in x
    arag_or_jls = runs[runs['Sample Name'].apply(is_arag_or_jls)]


    zero_conc = lambda x : 0 if x <0 else x
    runs['inferred_mol/mol'] = \
        runs['inferred_mol/mol'].apply(zero_conc)

    runs = runs.reset_index(drop=True)
    return runs

def get_data_ppm(include_run1=False):
    if include_run1:
        run1 = read_and_format_data(run_number=1)
    run4 = read_and_format_data(run_number=4)
    run5 = read_and_format_data(run_number=5)
    run6 = read_and_format_data(run_number=6)
    run7 = read_and_format_data(run_number=7)


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
    
    def find_concentrations(run, run_number):
        # print(run_number)
        run = blank_correction(run, var='CPS')
        run = normalize_to_In(run, var='CPS_sample_blank_corrected')
        standards= 'old' if run_number == 1 else 'new'
        make_calibration_curves(run, var='CPS_In_normalized_sample',
            standards=standards)
        run = add_iodine_cutoff_values(run)
        run['run_id'] = str(run_number)
        return run

    runs = [run4, run5, run6, run7]
    run_numbers = [4,5,6,7]

    if include_run1:
        runs = [run1] + runs
        run_numbers = [1] + run_numbers

    runs = [find_concentrations(run, run_numbers[i]) for i, run in enumerate(runs)]
    runs = pd.concat(runs)
    # runs = runs[runs.is_sample]
    runs = runs[~runs.is_sample_blank]

    zero_conc = lambda x : 0 if x <0 else x
    runs['inferred_ppb_zero_intercept'] = \
        runs['inferred_ppb_zero_intercept'].apply(zero_conc)
    
    def assign_formation(sample_name):
        if sample_name.startswith('SR'):
            return 'Steep Rock'
        elif sample_name.startswith('B'):
            return 'Gordon Lake'
        elif sample_name.startswith('18'):
            return 'Espanola'
        else:
            return 'Other'
        
    # assign formation names
    runs['formation'] = runs['Sample Name'].apply(assign_formation)

    runs = runs.reset_index(drop=True)
    return runs


def merge_ICPMS_data_long():
    run4 = read_and_format_data(run_number=4)
    run5 = read_and_format_data(run_number=5)
    run6 = read_and_format_data(run_number=6)
    run7 = read_and_format_data(run_number=7)

    def find_concentrations(run, run_number):
        run = blank_correction(run, var='CPS')
        run = normalize_to_Ca(run, run_num=run_number)
        run = make_calibration_curves_Ca(run, var='CPS_Ca_normalized_sample')
        run['run_id'] = str(run_number)
        return run

    runs = [run4, run5, run6, run7]
    run_numbers = [4,5,6,7]
    runs = [find_concentrations(run, run_numbers[i]) for i, run in enumerate(runs)]
    runs = pd.concat(runs)
    runs = runs[runs.is_sample]
    runs = runs[~runs.is_sample_blank]

    zero_conc = lambda x : 0 if x <0 else x
    runs['inferred_mol/mol_zero_intercept'] = \
        runs['inferred_mol/mol_zero_intercept'].apply(zero_conc)
    
    def assign_formation(sample_name):
        if sample_name.startswith('SR'):
            return 'Steep Rock'
        elif sample_name.startswith('B'):
            return 'Gordon Lake'
        else:
            return 'Espanola'
        
    # assign formation names
    runs['formation'] = runs['Sample Name'].apply(assign_formation)

    runs.to_csv('concentrations_mol_per_mol_long.csv')


def merge_icpms_data_ppb():
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

def get_dataset(formation='Gordon Lake'):
    # read concentration data
    conc_data = pd.read_csv('concentrations_mol_per_mol_Ca_Mg.csv')
    # cahnge espanola sample names
    change_names = lambda x: x.replace('18-2-', '')
    conc_data['Sample Name'] = conc_data['Sample Name'].apply(change_names)

    def assign_formation(sample_name):
        if sample_name.startswith('SR'):
            return 'Steep Rock'
        elif sample_name.startswith('B'):
            return 'Gordon Lake'
        else:
            return 'Espanola'
        
    # assign formation names
    conc_data['formation'] = conc_data['Sample Name'].apply(assign_formation)

    # calculate total REE
    REE = ['La','Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
        'Tm', 'Yb', 'Lu']
    
    REE_cols = []
    for col in conc_data.columns:
        for ree in REE:
            if ree in col:
                REE_cols.append(col)

    total_REE = 0
    for col in REE_cols:
        total_REE = total_REE + conc_data[col]

    conc_data['total_REE'] = total_REE

    # calculate Ce anomalies
    Ce, Pr, Nd = '140  Ce  [ No Gas ]', '141  Pr  [ No Gas ]', '146  Nd  [ No Gas ]'

    # normalize to PAAS
    paas = get_PAAS_data()
    paas_Ce = paas[paas.element == 'Ce']['paas/Ca+Mg_mol_per_mol'].values[0]
    paas_Nd = paas[paas.element == 'Nd']['paas/Ca+Mg_mol_per_mol'].values[0]
    paas_Pr = paas[paas.element == 'Pr']['paas/Ca+Mg_mol_per_mol'].values[0]
    Ce_norm = conc_data[Ce]/paas_Ce
    Nd_morm = conc_data[Nd]/paas_Nd
    Pr_morm = conc_data[Pr]/paas_Pr

    # calculate cerium anomalies. Eqs 7, 9 from Barrat 2023
    # geometric extrapolation
    conc_data['Ce_anomaly_geom'] = Ce_norm / (Pr_morm**2 / Nd_morm) 
    # linear extrapolation
    conc_data['Ce_anomaly_lin'] = Ce_norm / (2 * Pr_morm - Nd_morm) 
    conc_data['neg_Ce_anomaly'] = conc_data['Ce_anomaly_geom'] < 0.9
    conc_data['pos_Ce_anomaly'] = conc_data['Ce_anomaly_geom'] > 1.1

    # iodine flag
    iodine_samples = np.loadtxt('iodine_samples', dtype='str', delimiter=',')
    has_iodine = lambda x: True if x in iodine_samples else False
    conc_data['has_iodine'] = conc_data['Sample Name'].apply(has_iodine)

    # carbonate richness flag
    Al_str = '27  Al  [ No Gas ]'
    def is_carbonate_rich(Al_conc):
        if formation == 'Gordon Lake':
            if Al_conc < 0.01:
                return True
        elif formation == 'Espanola':
            if Al_conc < 0.04:
                return True
        elif formation == 'Steep Rock':
            return True
        return False
    conc_data['is_carbonate_rich'] = conc_data[Al_str].apply(is_carbonate_rich)

    # Carbonate classification based on Mg/Ca ratio

    conc_data['Mg/Ca'] = conc_data['25 -> 25  Mg  [ H2 ]'] /\
        conc_data['44  Ca  [ No Gas ]']
    
    def classify_carbonates(MgCa_ratio):
        if MgCa_ratio >= 0.8:
            return 'Dolomite'
        elif MgCa_ratio >= 0.3 and MgCa_ratio < 0.8:
            return 'Ankerite'
        else:
            return 'Calcite'
    
    conc_data['carbonate_type'] = conc_data['Mg/Ca'].apply(classify_carbonates)

    # get isotope data
    isotope_data =get_isotope_data()
    isotope_data['formation'] = isotope_data['Sample Name'].apply(assign_formation)

    # filter formation data
    conc_data = conc_data[conc_data.formation == formation]
    isotope_data = isotope_data[isotope_data.formation == formation]

    if formation == 'Steep Rock':
        data = conc_data[conc_data.formation == 'Steep Rock']
    else:
        # merge datafiles
        data = pd.merge(isotope_data, conc_data,
            on=['Sample Name', 'formation'], how='inner')
    return data



def normalize_to_In_old(data):
    # normalize to In
    In_ref_cps_hno3 = data[data.element_description ==\
        '115  In  [ No Gas ]']['CPS_HNO3_blank_mean'].iloc[0]
    data['CPS_115_In_NoGas_HNO3_blank_mean'] = In_ref_cps_hno3

    In_ref_cps_spl = data[data.element_description ==\
        '115  In  [ No Gas ]']['CPS_sample_blank_mean'].iloc[0]
    data['CPS_115_In_NoGas_sample_blank_mean'] = In_ref_cps_spl
    
    data_In = data[data.element_description == '115  In  [ No Gas ]']
    data_In = data_In[['CPS', 'complete_sample_name']]
    data = pd.merge(data, data_In, on='complete_sample_name',
        how='left', suffixes=('','_115_In_No_Gas'))

    data['CPS_In_normalized_HNO3'] =\
        data['CPS_HNO3_blank_corrected'] * (In_ref_cps_hno3/data['CPS_115_In_No_Gas'])
    data['CPS_In_normalized_sample'] =\
        data['CPS_sample_blank_corrected'] * (In_ref_cps_spl/data['CPS_115_In_No_Gas'])
    return data


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

    

# def normalize_to_Ca(data):
#     data_Ca = data[data.element_description == '44  Ca  [ No Gas ]']
#     data_Ca_hno3 = data_Ca[['CPS_HNO3_blank_corrected', 'complete_sample_name']]
#     data_Ca_spl = data_Ca[['CPS_sample_blank_corrected', 'complete_sample_name']]
#     data_Ca_cps = data_Ca[['CPS', 'complete_sample_name']]

#     data = pd.merge(data, data_Ca_hno3, on='complete_sample_name',
#         how='left', suffixes=('','_44_Ca'))
#     data = pd.merge(data, data_Ca_spl, on='complete_sample_name',
#         how='left', suffixes=('','_44_Ca'))
#     data = pd.merge(data, data_Ca_cps, on='complete_sample_name',
#         how='left', suffixes=('','_44_Ca'))
    
#     data['CPS_Ca_normalized_HNO3'] =\
#         data['CPS_HNO3_blank_corrected'] / data['CPS_HNO3_blank_corrected_44_Ca']
#     data['CPS_Ca_normalized_sample'] =\
#         data['CPS_sample_blank_corrected'] / data['CPS_sample_blank_corrected_44_Ca']
#     data['CPS_Ca_normalized'] =\
#         data['CPS'] / data['CPS_44_Ca']
#     return data

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


    







merge_icpms_data()
merge_icpms_data_ppb()

