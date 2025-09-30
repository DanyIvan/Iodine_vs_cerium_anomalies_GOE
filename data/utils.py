import pandas as pd
import numpy as np
import os
from molar_mass import molar_masses

def read_and_format_data(run_number=3, start=0):
    if run_number == 1:
        data = pd.read_excel('raw_data/ICPMS-RUN-DG-20230919.xlsx',
            sheet_name='cps averages')
    if run_number == 2:
        data = pd.read_excel('raw_data/ICPMS-RUN-DG-20230926.xlsx',
            sheet_name='cps averages')
    if run_number == 3:
        data = pd.read_excel('raw_data/ICPMS-RUN-DG-20230928.xlsx',
            sheet_name='cps averages', skiprows=list(range(2,17)))
    if run_number == 4:
        data = pd.read_excel('raw_data/ICPMS-RUN-DG-20231118.xlsx',
            sheet_name='cps averages')
    if run_number == 5:
        data = pd.read_excel('raw_data/ICPMS-RUN-DG-20231123.xlsx',
            sheet_name='cps averages')
    if run_number == 6:
        data = pd.read_excel('raw_data/ICPMS-RUN-DG-20231124.xlsx',
            sheet_name='cps averages')
    if run_number == 7:
        data = pd.read_excel('raw_data/ICPMS-RUN-DG-20231201.xlsx',
            sheet_name='cps averages')
    if run_number == 8:
        data = pd.read_excel('raw_data/ICPMS-RUN-DG-20240530.xlsx',
            sheet_name='cps averages')
    if run_number == 9:
        data = pd.read_excel('raw_data/ICPMS-RUN-DG-20240605.xlsx',
            sheet_name='cps averages')
    if run_number == 10:
        data = pd.read_excel('raw_data/ICPMS-RUN-DG-20240606.xlsx',
            sheet_name='cps averages')
    if run_number == 11:
        data = pd.read_excel('raw_data/ICPMS-RUN-DG-20240607.xlsx',
            sheet_name='cps averages')
        

    # format data
    elements = [x.strip() for x in data.columns if 'Unnamed' not in x][1:]
    info_cols = [2,5,6,7] if run_number < 8 else [3,4]
    info = data.iloc[1:, info_cols]
    info.columns = data.iloc[0,info_cols].values
    info['run_number'] = info.index

    # stack data to have a long format
    columns = ['CPS', 'CPS_RSD', 'Det.'] if run_number != 8 else ['CPS', 'CPS_RSD']
    ncols = len(columns)
    long_data = []
    i = 9 if run_number < 8 else 6
    j=0
    while i < data.shape[1]:
        # print(data.columns[i])
        data_element = data.iloc[1:,i:i+ncols]
        data_element.columns = columns
        data_element = pd.concat([info, data_element], axis=1)
        atomic_number, element, gas =  elements[j].split('  ')
        data_element['atomic_number'] = atomic_number
        data_element['element'] = element
        data_element['gas'] = gas
        data_element['element_description'] = elements[j]
        long_data.append(data_element)
        i += ncols
        j += 1

    long_data = pd.concat(long_data)
    # long_data['run_number'] = long_data['Data File'].apply(lambda x:x[:3])
    # long_data['run_number_int'] = long_data['run_number'].astype(int)
    long_data['complete_sample_name'] = long_data['Sample Name'] +\
        '_' + long_data['Comment'].fillna('') + '_' + long_data['run_number'].astype(str)

    def is_blank(name):
        if 'blank' in name or '2% HNO3' in name or 'Stabilizer' in name:
            return True
        else:
            return False
        
    def is_sample_blank(name):
        if ('DG' in name and 'blank' in name) or 'blank' in name:
            return True
        else:
            return False
    
    long_data['is_blank'] = long_data.complete_sample_name.apply(is_blank)
    long_data['is_sample_blank'] = long_data.complete_sample_name.apply(is_sample_blank)
    long_data = long_data.astype({'CPS':float, 'CPS_RSD':float})

    s1 = read_standards202309()
    s2 = read_standards202311()
    s3 = read_standards202405()
    standards = pd.concat([s1, s2, s3]).Name.unique()

    is_standard = lambda x: True if x in standards else False
    is_sample = lambda x: not (is_standard(x) or is_blank(x))
    long_data['is_standard'] = long_data['Sample Name'].apply(is_standard)
    long_data['is_sample'] = long_data['Sample Name'].apply(is_sample)
    return long_data

# read standards
def read_standards202309(with_Ca=False):
    path = 'calibation_standards/standards_202309/'
    folder = f'{path}/with_Ca' if with_Ca else f'{path}/without_Ca/'
    files = os.listdir(folder)
    standards = []
    for f in files:
        data = pd.read_csv(folder + f, sep='\t')

        if with_Ca:
             # convert to mol/L
            ppb_to_mol_l = lambda el: 1e6 * data[el]/molar_masses[el]
            for col in data.columns[1:]:
                data[col] = ppb_to_mol_l(col)
            # normalize to Ca 
            # data['Ca+Mg']  = data['Ca'] + data['Mg']
            data[data.columns[1:]] = data[data.columns[1:]].div(data['Ca'], axis=0)
        data = pd.melt(data, id_vars=['Name'],
            value_vars=data.columns[1:],
            var_name='element', value_name='ppb')
        standards.append(data)
    standards = pd.concat(standards)

    if not with_Ca:
        CaCal = pd.DataFrame(
            {'Name': ['Ca Cal'], 'element':['Ca'], 'ppb':[411720.0304]})
        MgCal = pd.DataFrame(
            {'Name': ['Mg Cal'], 'element':['Mg'], 'ppb':[253148.3321]})
        standards = pd.concat([standards, CaCal], ignore_index=True)
        standards = pd.concat([standards, MgCal], ignore_index=True)
    else:
        standards = standards.rename(columns={'ppb':'mol/mol_Ca'})
    return standards


def read_standards202311():
    folder = 'calibation_standards/standards_202311/' 
    files = os.listdir(folder)
    standards = []
    for f in files:
        data = pd.read_csv(folder + f, sep='\t')
        # convert to mol/L. ppb = 1 micro gram/ l
        ppb_to_mol_l = lambda el: 1e6 * data[el]/molar_masses[el]
        for col in data.columns[1:]:
            data[col] = ppb_to_mol_l(col)
        # normalize to Ca
        # data['Ca+Mg']  = data['Ca'] + data['Mg']
        data[data.columns[1:]] = data[data.columns[1:]].div(data['Ca'], axis=0)
        data = pd.melt(data, id_vars=['Name'],
            value_vars=data.columns[1:],
            var_name='element', value_name='ppb')
        standards.append(data)
    standards = pd.concat(standards)

    standards = standards.rename(columns={'ppb':'mol/mol_Ca'})
    return standards

def read_standards202311_ppb():
    folder = 'calibation_standards/standards_202311/'
    files = os.listdir(folder)
    standards = []
    for f in files:
        data = pd.read_csv(folder + f, sep='\t')
        data = pd.melt(data, id_vars=['Name'],
            value_vars=data.columns[1:],
            var_name='element', value_name='ppb')
        standards.append(data)
    standards = pd.concat(standards)
    return standards


def read_standards202405(run=8, only_Ca=False):
    if run == 8:
        filename = 'iodine_std_20240530.tsv'
    elif run == 9:
        filename = 'iodine_std_20240604.tsv'
    elif run == 10:
        filename = 'iodine_std_20240605.tsv'
    elif run == 11:
        filename = 'iodine_std_20240606.tsv'
    path = f'calibation_standards/standards_202405/{filename}' 

    data = pd.read_csv(path, sep='\t')
    # convert to mol/L. ppb = 1 micro gram/ l
    ppb_to_mol_l = lambda el: 1e6 * data[el]/molar_masses[el]
    for col in data.columns[1:]:
        data[col] = ppb_to_mol_l(col)
    
    if only_Ca:
        # normalize to Ca 
        data['norm']  = data['Ca']
    else:
        # normalize to Ca + Mg
        data['norm']  = data['Ca'] + data['Mg']
    data[data.columns[1:]] = data[data.columns[1:]].div(data['norm'], axis=0)
    data = pd.melt(data, id_vars=['Name'],
        value_vars=data.columns[1:],
        var_name='element', value_name='know_mol/mol')
    return data

def read_standards202405_ppb(run=8):
    if run == 8:
        filename = 'iodine_std_20240530.tsv'
    elif run == 9:
        filename = 'iodine_std_20240604.tsv'
    elif run == 10:
        filename = 'iodine_std_20240605.tsv'
    elif run == 11:
        filename = 'iodine_std_20240606.tsv'
    path = f'calibation_standards/standards_202405/{filename}' 

    data = pd.read_csv(path, sep='\t')
    data = pd.melt(data, id_vars=['Name'],
        value_vars=data.columns[1:],
        var_name='element', value_name='ppb')
    return data

def blank_correction(data, var='CPS'):
    blanks_HNO3 = data[data['Sample Name'] == '2% HNO3']
    blanks_samples = data[data.is_sample_blank]

    blanks_HNO3_mean_cps = blanks_HNO3[['element_description', var]].groupby('element_description').mean().reset_index()
    blanks_HNO3_mean_cps = blanks_HNO3_mean_cps[[ 'element_description', var]]

    blanks_spl_mean_cps = blanks_samples[['element_description', var]].groupby('element_description').mean().reset_index()
    blanks_spl_mean_cps = blanks_spl_mean_cps[[ 'element_description', var]]

    # merge blank average
    data = pd.merge(data, blanks_HNO3_mean_cps, on='element_description',
        how='left', suffixes=('','_HNO3_blank_mean'))
    data = pd.merge(data, blanks_spl_mean_cps, on='element_description',
        how='left', suffixes=('','_sample_blank_mean'))

    # make blank corrections
    data[f'{var}_HNO3_blank_corrected'] = data[var] - data[f'{var}_HNO3_blank_mean']
    data[f'{var}_sample_blank_corrected'] = data[var] - data[f'{var}_sample_blank_mean']
    return data

def get_isotope_data():
    data = pd.read_csv('isotope_data.csv')
    data_steep = pd.read_csv('steep_rock_isotopes.csv')
    # drop duplicates keeping first:
    duplicates = data.drill_core[data.drill_core.duplicated()]
    duplicates_idxs = duplicates.index
    data = data.drop(duplicates_idxs)

    data = data[['drill_core', 'd13c_tp', 'd18o_tp', 'texture']]
    data = data.rename(columns={'drill_core': 'Sample Name'})

    data_steep = data_steep[['Sample', 'd13C', 'd18O', 'Textures']]
    data_steep.columns = ['Sample Name', 'd13c_tp', 'd18o_tp', 'texture']

    data = pd.concat([data, data_steep])
    return data

def get_PAAS_data():
    paas = pd.read_csv('PAAS_TaylorMcLennan1985.csv')

    # get CaO concentration
    CaO_pct = float(paas[paas['var']=='CaO']['value'].iloc[0]) # units: %
    # convert to ppm (1pct = 1e4 ppm)
    CaO_ppm = CaO_pct * 1e4
    # get Ca ppm
    Ca_ppm = CaO_ppm * (molar_masses['Ca']/(molar_masses['Ca'] + molar_masses['O']))

    # get MgO concentration
    MgO_pct = float(paas[paas['var']=='MgO']['value'].iloc[0]) # units: %
    # convert to ppm
    MgO_ppm = MgO_pct * 1e4
    # get Mg ppm
    Mg_ppm = MgO_ppm * (molar_masses['Mg']/(molar_masses['Mg'] + molar_masses['O']))

    # get Al2O3 concentration
    Al2O3_pct = float(paas[paas['var']=='Al2O3']['value'].iloc[0]) # units: %
    # convert to ppm (1pct = 1e4 ppm)
    Al2O3_ppm = Al2O3_pct * 1e4
    # get Al ppm
    Al_ppm = Al2O3_ppm * ((2*molar_masses['Al'])/(2*molar_masses['Al'] + 3*molar_masses['O']))

    # filter REE data
    elements = ['La','Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb',
                'Dy', 'Y', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    element_order = {elements[i]:i for i in range(len(elements))}

    paas = paas[paas['var'].isin(elements)]
    paas = paas[['var', 'value']]
    paas = paas.rename(columns={'var':'element', 'value':'ppm'})
    paas = paas.astype({'ppm':float})

    # convert to mol/L. 1ppm = 1 mg/l; 1mg/l = 1e-3g/l
    Ca_mol_per_l =  1e-3 * Ca_ppm/molar_masses['Ca']
    Mg_mol_per_l = 1e-3 * Mg_ppm/molar_masses['Mg']
    Al_mol_per_l = 1e-3 * Al_ppm/molar_masses['Al']
 
    paas['molar_mass'] = paas.element.map(molar_masses)
    paas['mol_per_l'] = 1e-3 * paas.ppm / paas.molar_mass

    # normalize to Ca + Mg
    paas['paas/Ca+Mg_mol_per_mol'] = paas.mol_per_l / (Ca_mol_per_l + Mg_mol_per_l)
    Al_over_Ca_Mg= Al_mol_per_l / (Ca_mol_per_l + Mg_mol_per_l)

    # normalize to Al
    paas['paas/Ca+Mg/Al_mol_per_mol'] = paas['paas/Ca+Mg_mol_per_mol']/ Al_over_Ca_Mg

    paas['Al_ppm'] = Al_ppm
    paas['Ca_ppm'] = Ca_ppm
    paas['Mg_ppm'] = Mg_ppm

    # sort elements
    paas['element_oder'] = paas['element'].map(element_order)
    paas = paas.sort_values('element_oder')
    return paas