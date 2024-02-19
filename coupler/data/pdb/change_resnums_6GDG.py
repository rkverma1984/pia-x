# Need change residue numbers that start from  SER A 191 to SER A 205

f = open('6GDG_org.pdb','r')
lines=f.read().split('\n')
f.close()

for i in range(len(lines)):
    line = lines[i]
    if line.startswith('ATOM'):
        if line.split()[4] == 'D':
            if 'N   ARG D 255' in line:
                starting_index = i
                break

ga_region_indexes = []
for i in range(len(lines)):
    line = lines[i]
    if line.startswith('ATOM'):
        if line.split()[4] == 'D':
            ga_region_indexes.append(i)

ending_index = ga_region_indexes[-1]


f = open('6GDG.pdb','w')
for i in range(starting_index):
    f.write(lines[i]+'\n')
for i in range(starting_index, ending_index+2):
    new_line = lines[i][:22]+"%4d"%(int(lines[i][22:26])+10)+lines[i][26:]
    f.write(new_line+'\n')
for i in range(ending_index+2,len(lines)-1):
    f.write(lines[i]+'\n')
f.close()

