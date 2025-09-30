import pandas as pd

# read data in mol per mol units
data_molpermol = pd.read_csv('../sample_data_mol_per_mol.csv')

excel_file = '../sample_datasheet.xlsx'
writer = pd.ExcelWriter(excel_file, engine='xlsxwriter')

def get_element(element_description):
    p2 = element_description.find('[')
    p1 = p2-5
    return element_description[p1:p2].strip()

# elements to include in datasheet
included_elements = ['127  I  [ No Gas solution ]', '137  Ba  [ No Gas ]', 
        '139  La  [ No Gas ]', '140  Ce  [ No Gas ]',
       '141  Pr  [ No Gas ]', '146  Nd  [ No Gas ]', '147  Sm  [ No Gas ]',
        '151  Eu  [ No Gas ] oxide_corrected', 
        '157  Gd  [ No Gas ] oxide_corrected', 
       '159  Tb  [ No Gas ]', '163  Dy  [ No Gas ]', '165  Ho  [ No Gas ]',
       '166  Er  [ No Gas ]', '169  Tm  [ No Gas ]', '172  Yb  [ No Gas ]',
       '175  Lu  [ No Gas ]', '55 -> 55  Mn  [ H2 ]', '88  Sr  [ No Gas ]',
       '27  Al  [ No Gas ]', '89  Y  [ No Gas ]', '25 -> 25  Mg  [ H2 ]' ]

# additional info
all_columns = ['Sample Name', 'formation'] + included_elements + ['total_REE', 
    'Ce_anomaly_geom', 'Ce_anomaly_lin', 'Y/Ho', 'd13c_sp', 'd18o_sp', 
    'has_iodine']

# rename columns
new_columns = ['Sample Name', 'formation'] +\
    [get_element(x) for x in included_elements] + ['total_REE', 
    'Ce_anomaly_geom', 'Ce_anomaly_lin', 'Y/Ho', 'd13c_sp', 'd18o_sp', 
    'has_iodine_above_detection_limit']

# save to excel sheet
data_molpermol = data_molpermol[all_columns]
data_molpermol.columns = new_columns
data_molpermol.to_excel(writer, sheet_name='mol_per_mol', index=False)


samples = data_molpermol['Sample Name'].unique()
change_names = lambda x: x.replace('18-2-', '')

# read REE data in ppb units
data_ppb_nitric = pd.read_csv('nitric_acid_runs_reduced_ppb.csv')
# use same samples as mol per mol 
data_ppb_nitric = data_ppb_nitric[data_ppb_nitric['Sample Name'].isin(samples)]
# change Espanola samples names so that they have same names as mol per mol data
data_ppb_nitric['Sample Name'] = data_ppb_nitric['Sample Name'].apply(change_names)

# read I data in ppb units
data_ppb_stabilizer = pd.read_csv('stabilizer_runs_reduced_ppb.csv')
# use same samples as mol per mol 
data_ppb_stabilizer['Sample Name'] = data_ppb_stabilizer['Sample Name'].apply(change_names)
# change Espanola samples names so that they have same names as mol per mol data
data_ppb_stabilizer = data_ppb_stabilizer[
    data_ppb_stabilizer['Sample Name'].isin(samples)]
# use only I data
data_ppb_stabilizer = data_ppb_stabilizer[
    data_ppb_stabilizer.element_description == '127  I  [ No Gas solution ]']

# merge I data to REE data
data_ppb_stabilizer = data_ppb_stabilizer[['Sample Name',
    'inferred_ppb_zero_intercept']]
data_ppb_stabilizer.columns = ['Sample Name', '127  I  [ No Gas solution ]']
data_ppb = pd.merge(data_ppb_nitric, data_ppb_stabilizer, how='outer')


# save to excel sheet
new_columns = ['Sample Name', 'formation'] +\
    [get_element(x) for x in included_elements] 
data_ppb = data_ppb[['Sample Name', 'formation'] + included_elements]
data_ppb.columns = new_columns
data_ppb.to_excel(writer, sheet_name='ppb', index=False)
writer.close()
