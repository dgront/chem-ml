from tdc.utils import retrieve_label_name_list

# ---------- labels for the hERG Central dataset
label_list = retrieve_label_name_list('herg_central')
print(label_list)

# ---------- Download Ames Mutagenicity dataset
from tdc.single_pred import Tox
ames = Tox(name = 'AMES')
data = ames.get_data()
print(data)

# ---------- Download the ADME datasets
# ---------- see https://tdcommons.ai/single_pred_tasks/adme/#metabolism
from tdc.single_pred import ADME
ADME_NAMES = ['CYP2D6_Veith', 'CYP2C19_Veith', 'CYP3A4_Veith','CYP1A2_Veith', 'CYP2C9_Veith']
for data_name in ADME_NAMES:
    data = ADME(name = data_name)
    df = data.get_data()
    print('Dataset: ', data.name)
    print(df)
    active = int(df['Y'].sum())
    all = int(df['Y'].count())
    print('all :', all)
    print('active:', active ,';',round((active / all)*100,2), '%', '\n')
