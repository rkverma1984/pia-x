
def change_id(pdb_file,org_chain_id, new_chain_id):
    f = open(pdb_file,'r')
    lines=f.read().split('\n')
    f.close()

    f = open('new_'+pdb_file,'w')
    for line in lines:
        if line.startswith('ATOM'):
            if line[21] == org_chain_id:
                new_line = line[:21]+new_chain_id+line[22:]
            else:
                new_line = line     
            f.write(new_line+'\n')
    f.close()
    return


