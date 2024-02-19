f = open('conserved_50_positions.txt','r')
lines=f.read().split('\n')
f.close()
tm_50 = {}
for line in lines:
    if line!='':
        if not line.startswith('#'):
            tm_50[line.split('\t')[0].replace('-','_')] = []
            for i in range(7):
                tm_50[line.split('\t')[0].replace('-','_')].append(line.split('\t')[1:][i][1:])

f = open('boundary_information.txt','r')
lines=f.read().split('\n')
f.close()
tm_boundary = []
for item in lines[-2].split()[1:]:
    tm_boundary.append(item.split('-')[0][4:]) # last item

uniprot_info = {}
for line in lines:
    if line!='':
        if not line.startswith('#'):
            uniprot_info[line.split()[0].replace('-','_')] = []
            for i in range(7):
                uniprot_info[line.split()[0].replace('-','_')].append([line.split()[i+1].split(':')[0][1:], line.split()[i+1].split(':')[1][1:], tm_boundary[i]])


for key in uniprot_info.keys():
    f = open(key+'_table.py','w')
    f.write('#!/usr/bin/env python\n')
    f.write('"""\n')
    f.write(key+' TMs table\n')
    f.write('"""\n')
    f.write('def convert_Residue2TMs():\n')
    f.write('    '+key+'_BW_dic = {}\n')

    for i in range(7):
        f.write('    for i in range('+uniprot_info[key][i][0]+'-5, '+uniprot_info[key][i][1]+'+2): #TM'+str(i+1)+'\n')
        f.write("        "+key+"_BW_dic[i]='TM"+str(i+1)+".'+str("+uniprot_info[key][i][2]+'-5+i-('+uniprot_info[key][i][0]+'-5))\n\n')

    for i in range(7):
        f.write('    assert '+key+'_BW_dic['+tm_50[key][i]+'] == '+"'TM"+str(i+1)+'.50'+"'\n")
    f.write('\n')
    f.write('    return '+key+'_BW_dic\n') 
    f.close()



