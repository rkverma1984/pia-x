import pandas as pd
import os

df = pd.read_excel('GPCRdb_complex_structures-2020-10-24.xlsx')
keys = df.keys()
data_length = df.shape[0]

#uniprot = df['Unnamed: 1']
#pdb_name = df['Unnamed: 7']
#pdb_refined_name = df['Unnamed: 8']
#preferred_chain = df['Unnamed: 10']
#species = df['Unnamed: 5']
#family = df['Unnamed: 3']
#H5_note = df['Unnamed: 15']

f = open('extracted_information_from_excel.txt','w')
header = '\t'.join(str(df[key][0]) for key in keys) +'\n'+'\t'.join(str(df[key][1]) for key in keys)  
f.write(header+'\n')

for i in range(2, data_length):
    if str(df[keys[1]][i]) != 'nan':
        tem = []
        for key in keys:
            if key == 'Unnamed: 7': ### PDB ID
                if not len(str(df[key][i])) == 4:
                    df[key][i] = str(df['Unnamed: 8'][i]).split('_')[0]  ### PDB refined ID
            if key == 'Unnamed: 3': ### Family name 
                org_family_name = str(df[key][i])
                family_name = org_family_name.replace(' ','_')
                family_name = family_name.replace('(','')
                family_name = family_name.replace(')','')
                df[key][i] = str(family_name)
                   
        outstring = '\t'.join(str(df[key][i]) for key in keys)
        f.write(outstring+'\n')
f.close()


