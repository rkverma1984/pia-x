f = open('../alignment/Galpha_alignment.index','r')
lines=f.read().split('\n')
f.close()
raw_data = {}
for i in range(len(lines)/2):
    raw_data[lines[i*2][1:]] = lines[i*2+1][:-1]
    hash_num = []
    for l in range(len(lines[i*2+1][:-1].split(','))):
        if not lines[i*2+1][:-1].split(',')[l] == '-':
            hash_num.append((int(lines[i*2+1][:-1].split(',')[l]), l)) # tuple type, 1st item is real resnum. 2nd item is aligned index
    
    f = open(lines[i*2][1:]+'_table.py','w')
    f.write('#!/usr/bin/env python\n')
    f.write('"""\n')
    f.write('G'+lines[i*2][1:][3].lower()+' alpha subtype '+lines[i*2][1:].split('_')[1]+'  species aligned sequence table\n')
    f.write('"""\n')
    f.write('def convert_Ga_Residues():\n')
    f.write('    hash_table = {}\n')
    for h in hash_num:
        f.write('    hash_table['+str(h[0])+'] = str('+str(h[1])+')\n')
    f.write('    return hash_table\n')
    f.close()

