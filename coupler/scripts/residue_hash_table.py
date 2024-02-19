def residue_name_3to1():
    three_letters = ['ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','SEC','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
    one_letter = ['R','H','K','D','E','S','T','N','Q','C','U','G','P','A','V','I','L','M','F','Y','W']
    AA_shortname = {}
    for i in range(len(three_letters)):
        AA_shortname[three_letters[i]] = one_letter[i]
    AA_shortname['ASH'] = 'D'
    AA_shortname['HIP'] = 'H'
    AA_shortname['HIE'] = 'H'
    AA_shortname['SEP'] = 'S'
    AA_shortname['TPO'] = 'T'
    return AA_shortname
