import glob,os

def _IO_pdb(pdb_file):
    f = open(pdb_file,'r')
    lines = f.read().split('\n')
    f.close()
    return lines

def _grab_missing_residues_remark(pdb_file):
    lines = _IO_pdb(pdb_file)
    remark_lines = []
    for line in lines:
        if ' MISSING ATOM  ' in line:
            remark_lines.append(line)

    remark_num = []
    if remark_lines != []:
        for r in list(set(remark_lines)):
            remark_num.append(r[:10])
    return remark_num

def extract_missing_residues(pdb_file, remark_tag):
    lines = _IO_pdb(pdb_file)
    indexes = []
    for i in range(len(lines)):
        if remark_tag in lines[i]:
            indexes.append(i)
    starting_index = None
    for index in indexes:
        if '   M RES CSSEQI  ATOMS' in lines[index]:
            starting_index = index
    missing_residues = []
    if starting_index != None:
        for i in range(starting_index+1, indexes[-1]+1):
            #print (lines[i])
            missing_residues.append(lines[i].split()[2]+' '+lines[i][19]+' '+lines[i][20:24].split()[0])
            #missing_residues.append(lines[i].split()[2]+' '+lines[i].split()[3]+' '+lines[i].split()[4])
    else:
        print ('Cannot find the header of REMARK MISSING RESIDUES.')

    return missing_residues

current_working_dir = os.getcwd()
os.chdir('/home/bxie/PycharmProjects/pia-x/coupler/data/pdb')

missing_dic = {}        
org_FNs = glob.glob('*/*/*_org.pdb')
_FNs = glob.glob('*/*/*.pdb')

FNs = []
for FN in _FNs:
    if FN[:-4]+'_org.pdb' in org_FNs:
        FNs.append(FN[:-4]+'_org.pdb')
    else:
        FNs.append(FN)

for FN in FNs:
    remark_nums = _grab_missing_residues_remark(FN)
    
    missing_residues = []
    for r in remark_nums:
        missing_residues += extract_missing_residues(FN, r)
     
    missing_dic[os.path.basename(FN)[:4]] = missing_residues

os.chdir(current_working_dir)

import pickle
with open('missing_residues.pkl','wb') as handle:
    pickle.dump(missing_dic, handle, protocol=pickle.HIGHEST_PROTOCOL)

