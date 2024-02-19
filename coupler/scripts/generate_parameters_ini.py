#f = open('extracted_information_from_excel.txt','r')
f = open('merged_excel_information.txt','r')
lines=f.read().split('\n')
f.close()
info = [line for line in lines if not line=='']
#PDB_ID uniprot_ID rec_chain_ID family_name Ga_name

f = open('heavy_atoms_parameters.ini','w')

f.write('''[DEFAULT]
top_pdb = data/sample.pdb
traj_dcd = None
receptor_chain_ID = R
ligand_chain_IDs = [A,B]
calc_type = heavy_atom
dist_threshold = 5
frequency_threshold = 0.0
stride = 1
bootstrap = False
receptor_resids = []
lig_resids = []
x_axis_label = G protein residues
y_axis_label = Receptor residues

''')

for i in info:
    #family_name = i.split()[1].upper()
    family_name = i.split()[4]
    uni_id = i.split()[1]
    pdb_id = i.split()[0]
    chain_id = i.split()[2]
    species = i.split()[3].lower()
    if 'turkey' in species:
        species = 'melga'
    ga_family_name = i.split()[5].upper()
    if not pdb_id in ['5G53','6GDG','6WHA','6G79']:
        f.write('''
[{0}]
top_pdb = ../data/pdb/{1}/{4}/{0}.pdb
receptor_chain_ID = {2}
family_name = {4}_{3}
ga_family_name = {5}
'''.format(pdb_id, family_name, chain_id, species, uni_id, ga_family_name))
    if pdb_id == '5G53':
        f.write('''
[{0}]
top_pdb = ../data/pdb/{1}/{4}/{0}.pdb
receptor_chain_ID = {2}
family_name = {4}_{3}
ga_family_name = {5}
ligand_chain_IDs = [C]
'''.format(pdb_id, family_name, chain_id, species, uni_id, ga_family_name))

    if pdb_id == '6GDG':
        f.write('''
[{0}]
top_pdb = ../data/pdb/{1}/{4}/{0}.pdb
receptor_chain_ID = R
family_name = {4}_{3}
ga_family_name = {5}
ligand_chain_IDs = [A,B]
'''.format(pdb_id, family_name, chain_id, species, uni_id, ga_family_name))

    if pdb_id == '6WHA':
        f.write('''
[{0}]
top_pdb = ../data/pdb/{1}/{4}/{0}.pdb
receptor_chain_ID = R
family_name = {4}_{3}
ga_family_name = {5}
ligand_chain_IDs = [A,B]
'''.format(pdb_id, family_name, chain_id, species, uni_id, ga_family_name))

    if pdb_id == '6G79':
        f.write('''
[{0}]
top_pdb = ../data/pdb/{1}/{4}/{0}.pdb
receptor_chain_ID = {2}
family_name = {4}_{3}
ga_family_name = {5}
ligand_chain_IDs = [A,B]
'''.format(pdb_id, family_name, chain_id, species, uni_id, ga_family_name))


f.close()


