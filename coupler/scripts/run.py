from vmd_contact import *
from wrap_toolkit import *
import residue_hash_table
import os
import sys
import importlib
from configparser import ConfigParser

#from contact_frequency.GPCR_PDBs.scripts.vmd_contact import wrap_pdb_distance
#from contact_frequency.GPCR_PDBs.scripts.wrap_toolkit import merge_AB_metrics

config = ConfigParser()
resi3_1=residue_hash_table.residue_name_3to1()

def _classify_paramter_type(parameter_section):
    """
    Define parameters type and ensemble them to a dictionary
    """
    dictionary = {}
    for key in parameter_section.keys():
        if 'threshold' in key:
            dictionary[key] = parameter_section.getfloat(key)
        elif 'stride' in key:
            dictionary[key] = parameter_section.getint(key)
        elif 'bootstrap' in key:
            dictionary[key] = parameter_section.getboolean(key)
        elif 'ligand_chain_ids' in key:
        #elif 'ligand_chain_ids' in key or 'receptor_resids' in key or 'lig_resids' in key:
            dictionary[key] = parameter_section[key][1:-1].split(',')
            if parameter_section[key][1:-1].split(',') == ['']:
                dictionary[key] = []
        elif 'receptor_resids' in key:
        #elif 'receptor_resids' in key or 'lig_resids' in key:
            if parameter_section[key][1:-1].split(',') == ['']:
                dictionary[key] = []
            else:
                #dictionary[key] = [int(p) for p in parameter_section[key][1:-1].split(',')]
                dictionary[key] = eval(parameter_section.get(key))
        elif 'lig_resids' in key:
            dictionary[key] = eval(parameter_section.get(key))
        else:
            if parameter_section[key] == 'None' :
                dictionary[key] = None
            else:
                dictionary[key] = parameter_section[key]
    #print (dictionary)
    return dictionary

def _load_parameters(config_path,pdb_id):
    """
    pdb_id: 4 letter code
    parameters file name must be XXX.ini 
    """
    config_file = os.path.join(config_path )
    if not os.path.exists(config_file):
        sys.exit('No parameters file found.')
    config.read(config_file)
    sections = config.sections()## PDB id
    default_parameters = config["DEFAULT"]
    default_dict = _classify_paramter_type(default_parameters)
    parameters_dict = _classify_paramter_type(config[pdb_id])

    return parameters_dict    


def _convert_receptor_label_full(receptor_labels, family_name):
    """
    Convert original residue index to BW number. Need family TM table first. Otherwise will use the original one.
    """
    tm_listfile = family_name+'_table'
    if not os.path.exists(tm_listfile+'.py'):
        print (tm_listfile, os.getcwd())
        print (' TM table could not find. Will use original residue numbers')
    TM_table = importlib.import_module(tm_listfile).convert_Residue2TMs()
    translated_rec_labels = []
    for r in receptor_labels:
        resnum = int(r.split()[0])
        if len(r.split()[1]) == 3:
            letter_code = resi3_1[r.split()[1]]
        else:
            letter_code = r.split()[1]
        if resnum in TM_table.keys():
            translated_rec_labels.append(TM_table[resnum]+' '+letter_code+' '+str(resnum))
        else:
            translated_rec_labels.append(str(resnum)+' '+letter_code+' '+str(resnum))

    return translated_rec_labels

def _convert_receptor_label(receptor_labels, family_name):
    """
    Convert original residue index to BW number. Need family TM table first. Otherwise will use the original one.
    """
    tm_listfile = family_name+'_table'
    if not os.path.exists(tm_listfile+'.py'):
        print (tm_listfile, os.getcwd())
        print (' TM table could not find. Will use original residue numbers')
    TM_table = importlib.import_module(tm_listfile).convert_Residue2TMs()    
    translated_rec_labels = []
    for r in receptor_labels:
        resnum = int(r.split()[0])
        if len(r.split()[1]) == 3:
            letter_code = resi3_1[r.split()[1]]
        else:
            letter_code = r.split()[1]
        if resnum in TM_table.keys():
            translated_rec_labels.append(TM_table[resnum]+' '+letter_code)
        else:
            translated_rec_labels.append(str(resnum)+' '+letter_code)

    return translated_rec_labels
    
def _convert_ligand_label(ligand_labels,ga_family_name):
    """
    Convert Gprotein or arrestin resname to 1 letter code.
    Labels can be : '27 A', or 'Gb 28 D'
    """

    translated_lig_labels = [] 
    ga_listfile = ga_family_name+'_table'
    if not os.path.exists(ga_listfile+'.py'):
        ga_convert = False
        print ('Ga index hash table could not find. Will use original residue numbers')
    else:
        ga_convert = True
        ga_table = importlib.import_module(ga_listfile).convert_Ga_Residues()

    for l in ligand_labels:
        if len(l.split())==2: # Ga
            if not ga_convert:
                resnum = int(l.split()[0])
                if len(l.split()[1]) ==3:
                    letter_code = resi3_1[l.split()[1]]
                else:
                    letter_code = l.split()[1]
                translated_lig_labels.append(str(resnum)+' '+letter_code)
            else:
                # Use aligned Ga index here 
                resnum = int(l.split()[0])
                if len(l.split()[1]) ==3:
                    letter_code = resi3_1[l.split()[1]]
                else:
                    letter_code = l.split()[1]
                translated_lig_labels.append(ga_table[resnum]+' '+letter_code)
                
        else: # label contains 3 items. GB
            resnum = int(l.split()[-2])
            if len(l.split()[-1]) == 3:
                letter_code = resi3_1[l.split()[-1]]
            else:
                letter_code = l.split()[-1]
            translated_lig_labels.append(' '.join(l.split()[:-2])+' '+str(resnum)+' '+letter_code)

    return translated_lig_labels

def _convert_ligand_label_full(ligand_labels,ga_family_name):
    """
    Convert Gprotein or arrestin resname to 1 letter code.
    Labels can be : '27 A', or 'Gb 28 D'
    """

    translated_lig_labels = []
    ga_listfile = ga_family_name+'_table'
    if not os.path.exists(ga_listfile+'.py'):
        ga_convert = False
        print ('Ga index hash table could not find. Will use original residue numbers')
    else:
        ga_convert = True
        ga_table = importlib.import_module(ga_listfile).convert_Ga_Residues()

    for l in ligand_labels:
        if len(l.split())==2: # Ga
            if not ga_convert:
                resnum = int(l.split()[0])
                if len(l.split()[1]) ==3:
                    letter_code = resi3_1[l.split()[1]]
                else:
                    letter_code = l.split()[1]
                translated_lig_labels.append(str(resnum)+' '+letter_code)
            else:
                # Use aligned Ga index here 
                resnum = int(l.split()[0])
                if len(l.split()[1]) ==3:
                    letter_code = resi3_1[l.split()[1]]
                else:
                    letter_code = l.split()[1]
                translated_lig_labels.append(ga_table[resnum]+' '+letter_code+' '+str(resnum))

        else: # label contains 3 items. GB
            resnum = int(l.split()[-2])
            if len(l.split()[-1]) == 3:
                letter_code = resi3_1[l.split()[-1]]
            else:
                letter_code = l.split()[-1]
            translated_lig_labels.append(' '.join(l.split()[:-2])+' '+str(resnum)+' '+letter_code)

    return translated_lig_labels



def _pdb_dist(**parameters):
    """
    A wrapper for pdb distance function. Support merging matrices. 
    ligand_chain_ids type must be list.
    """

    matrices = {}
    receptor_labels = {}
    ligand_labels = {}
    ligand_chain_ids = parameters['ligand_chain_ids']
    top_pdb = parameters['top_pdb']
    receptor_chain_id = parameters['receptor_chain_id']
    dist_threshold = parameters['dist_threshold']
    calc_type = parameters['calc_type']
    traj_dcd = parameters['traj_dcd']
    stride = parameters['stride']
    receptor_resids = parameters['receptor_resids']
    lig_resids = parameters['lig_resids']

    for ligand_chain_id in ligand_chain_ids:
        if lig_resids != {}:
            lig_resids_list = lig_resids[ligand_chain_id]
        else:
            lig_resids_list = []

        matrices[ligand_chain_id],receptor_labels[ligand_chain_id],ligand_labels[ligand_chain_id] = wrap_pdb_distance(top_pdb, receptor_chain_id, ligand_chain_id, dist_threshold, calc_type, traj_dcd, stride, receptor_resids, lig_resids_list)

        if ligand_chain_id == 'B':
            gb_tag_ligand_labels = []
            for l in ligand_labels[ligand_chain_id]:
                gb_tag_ligand_labels.append('Gb '+l)
            ligand_labels[ligand_chain_id] = gb_tag_ligand_labels

    merged_matrix = matrices[ligand_chain_ids[0]]
    merged_rec_labels = receptor_labels[ligand_chain_ids[0]]
    merged_lig_labels = ligand_labels[ligand_chain_ids[0]]

    if len(ligand_chain_ids) > 1:
        for i in range(0,len(ligand_chain_ids)-1):
            merged_matrix, merged_rec_labels, merged_lig_labels = merge_AB_metrics(merged_matrix, matrices[ligand_chain_ids[i+1]],
    merged_rec_labels,receptor_labels[ligand_chain_ids[i+1]],
    merged_lig_labels,ligand_labels[ligand_chain_ids[i+1]],dist_threshold)

    return merged_matrix, merged_rec_labels, merged_lig_labels


def _traj_freq(**parameters):
    matrices = {}
    receptor_labels = {}
    ligand_labels = {}
    ligand_chain_ids = parameters['ligand_chain_ids']
    top_pdb = parameters['top_pdb']
    receptor_chain_id = parameters['receptor_chain_id']
    dist_threshold = parameters['dist_threshold']
    frequency_threshold = parameters['frequency_threshold']
    calc_type = parameters['calc_type']
    traj_dcd = parameters['traj_dcd']
    stride = parameters['stride']
    receptor_resids = parameters['receptor_resids']
    lig_resids = parameters['lig_resids']
    bootstrap = parameters['bootstrap']

    for ligand_chain_id in ligand_chain_ids:
        if lig_resids != {}:
            lig_resids_list = lig_resids[ligand_chain_id]
        else:
            lig_resids_list = []

        matrices[ligand_chain_id],receptor_labels[ligand_chain_id],ligand_labels[ligand_chain_id] = wrap_traj_freq(top_pdb, traj_dcd, receptor_chain_id, ligand_chain_id, calc_type, dist_threshold, frequency_threshold, stride, bootstrap, receptor_resids, lig_resids_list)

        if ligand_chain_id == 'B':
            gb_tag_ligand_labels = []
            for l in ligand_labels[ligand_chain_id]:
                gb_tag_ligand_labels.append('Gb '+l)
            ligand_labels[ligand_chain_id] = gb_tag_ligand_labels

    merged_matrix = matrices[ligand_chain_ids[0]]
    merged_rec_labels = receptor_labels[ligand_chain_ids[0]]
    merged_lig_labels = ligand_labels[ligand_chain_ids[0]]

    if len(ligand_chain_ids) > 1:
        for i in range(0,len(ligand_chain_ids)-1):
            merged_matrix, merged_rec_labels, merged_lig_labels = merge_AB_metrics_freq(merged_matrix, matrices[ligand_chain_ids[i+1]],
    merged_rec_labels,receptor_labels[ligand_chain_ids[i+1]],
    merged_lig_labels,ligand_labels[ligand_chain_ids[i+1]],frequency_threshold)

    return merged_matrix, merged_rec_labels, merged_lig_labels 


def _traj_dist(**parameters):
    matrices = {}
    receptor_labels = {}
    ligand_labels = {}
    ligand_chain_ids = parameters['ligand_chain_ids']
    top_pdb = parameters['top_pdb']
    receptor_chain_id = parameters['receptor_chain_id']
    dist_threshold = parameters['dist_threshold']
    frequency_threshold = parameters['frequency_threshold']
    calc_type = parameters['calc_type']
    traj_dcd = parameters['traj_dcd']
    stride = parameters['stride']
    receptor_resids = parameters['receptor_resids']
    lig_resids = parameters['lig_resids']
    bootstrap = parameters['bootstrap']

    for ligand_chain_id in ligand_chain_ids:
        if lig_resids != {}:
            lig_resids_list = lig_resids[ligand_chain_id]
        else:
            lig_resids_list = []

        matrices[ligand_chain_id],receptor_labels[ligand_chain_id],ligand_labels[ligand_chain_id] = wrap_traj_distance(top_pdb, receptor_chain_id, ligand_chain_id, dist_threshold, calc_type, traj_dcd, stride, receptor_resids, lig_resids_list)

        #if ligand_chain_id == 'B':
        #    gb_tag_ligand_labels = []
        #    for l in ligand_labels[ligand_chain_id]:
        #        gb_tag_ligand_labels.append('Gb '+l)
        #    ligand_labels[ligand_chain_id] = gb_tag_ligand_labels

    #merged_matrix = matrices[ligand_chain_ids[0]]
    #merged_rec_labels = receptor_labels[ligand_chain_ids[0]]
    #merged_lig_labels = ligand_labels[ligand_chain_ids[0]]

    #if len(ligand_chain_ids) > 1:
    #    for i in range(0,len(ligand_chain_ids)-1):
    #        merged_matrix, merged_rec_labels, merged_lig_labels = merge_AB_metrics_freq(merged_matrix, matrices[ligand_chain_ids[i+1]],
    #merged_rec_labels,receptor_labels[ligand_chain_ids[i+1]],
    #merged_lig_labels,ligand_labels[ligand_chain_ids[i+1]],frequency_threshold)

    #return merged_matrix, merged_rec_labels, merged_lig_labels
    return matrices,receptor_labels,ligand_labels

    

def master(run_type, pdb_id, config_path, convert=True):
    """
    run_type: pdb_distance, traj_frequency
    """
    parameters = _load_parameters(config_path,pdb_id)
    x_axis_label = parameters['x_axis_label']
    y_axis_label = parameters['y_axis_label']

    if run_type == 'pdb_distance':
        matrix, xlabels, ylabels = _pdb_dist(**parameters)

    elif run_type == 'traj_frequency':
        matrix, xlabels, ylabels = _traj_freq(**parameters)

    elif run_type == 'traj_distance':
        matrix, xlabels, ylabels = _traj_dist(**parameters)  # matrix will not be the merged one (Ga+Gb), instead, it will be a dictionary.
    else:
        sys.exit('Unknown run type. Terminate.')

    if convert:        
        if xlabels:
            xlabels = _convert_receptor_label_full(xlabels, parameters['family_name'])
            #xlabels = _convert_receptor_label(xlabels, parameters['family_name'])
        if ylabels:
            ylabels = _convert_ligand_label_full(ylabels, parameters['ga_family_name'])
            #ylabels = _convert_ligand_label(ylabels, parameters['ga_family_name'])


    return matrix, xlabels, ylabels, x_axis_label, y_axis_label


def add_labels(pdb_id, config_path,org_xlabels,org_ylabels):
    parameters = _load_parameters(config_path,pdb_id)
    xlabels = _convert_receptor_label(org_xlabels, parameters['family_name'])
    #ylabels = _convert_ligand_label(org_ylabels, parameters['ga_family_name'])
    #return xlabels, ylabels
    return xlabels, org_ylabels
