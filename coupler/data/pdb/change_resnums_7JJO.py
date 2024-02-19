# Need change residue numbers that start from  SER A 191 to SER A 205

f = open('7JJO_org.pdb','r')
lines=f.read().split('\n')
f.close()

for i in range(len(lines)):
    line = lines[i]
    if line.startswith('ATOM'):
        if line.split()[4] == 'A':
            if 'N   SER A 191' in line:
                starting_index = i
                break

ga_region_indexes = []
for i in range(len(lines)):
    line = lines[i]
    if line.startswith('ATOM'):
        if line.split()[4] == 'A':
            ga_region_indexes.append(i)

ending_index = ga_region_indexes[-1]


f = open('7JJO.pdb','w')
for i in range(starting_index):
    f.write(lines[i]+'\n')
for i in range(starting_index, ending_index+2):
    new_line = lines[i][:22]+"%4d"%(int(lines[i][22:26])+14)+lines[i][26:]
    f.write(new_line+'\n')
for i in range(ending_index+2,len(lines)-1):
    f.write(lines[i]+'\n')
f.close()

