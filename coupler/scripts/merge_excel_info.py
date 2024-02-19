f = open('filtered_information.txt','r')
lines=f.read().split('\n')
f.close()
lines = [line for line in lines[2:] if line!='']
info = {}
for line in lines:
    extracted_information = line.split('\t')[1]+' '+line.split('\t')[10]+' '+line.split('\t')[5].replace(' ','_')+' '+line.split('\t')[3]
    info[line.split('\t')[7]] = extracted_information

f = open('57_structure_Ga_family_information.dat','r')
lines=f.read().split('\n')
f.close()
lines = [line for line in lines if line!='']
ga_info = {}
for line in lines:
    ga_info[line.split()[0]] = line.split()[4]

f = open('merged_excel_information.txt','w')
#PDB_ID Uniprot_name Chain_ID Species Familyname Ga_Familyname
for key in info.keys():
    f.write(key+' '+info[key]+' '+ga_info[key]+'\n')
f.close()


