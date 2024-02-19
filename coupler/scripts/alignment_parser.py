import os
fasta_filename = '../alignment/receptor_alignment.r2.fasta'
conserved50_list = [] # If .50 positions are known based on the aligned sequences, type them here.

# based on OPRM knowledge
# TM1.50 in MVTAITIMALYSIVCVVGLFGN
# TM2.50 in TLPFQ
# TM3.50 in TSIFTLCTMSVDR
# TM4.50 in IVNVCNW
# TM5.50 in CVFIFAFIMP
# TM6.50 in TRMVLVVVAVFIVCWTP
# TM7.50 in HFCIALGYTNSCLNP

boundary_based_on_paper = {'1':[30, 59],'2':[38, 67], '3':[22, 54], '4':[40, 62], '5':[35, 60], '6':[30, 60], '7':[33, 68]  }

#----------------- Based on the alignment fasta ------------------------------#
def decide_50_position(oprm_seq):
    """
    Based on OPRM. If the .50 positions are known, this step can be skipped.
    Return .50 real index
    """
    tm150_keyword = 'MVTAITIMALYSIVCVVGLFGN'
    tm250_keyword = 'NLALAD'
    tm350_keyword = 'TSIFTLCTMSVDR'
    tm450_keyword = 'IVNVCNW'
    tm550_keyword = 'CVFIFAFIMP'
    tm650_keyword = 'TRMVLVVVAVFIVCWTP'
    tm750_keyword = 'HFCIALGYTNSCLNP'
    
    index_150 = None
    index_250 = None
    index_350 = None
    index_450 = None
    index_550 = None
    index_650 = None
    index_750 = None
    indexes_conserved50 = []
    for i in range(len(oprm_seq)):
        if oprm_seq[i:i+len(tm150_keyword)] == tm150_keyword:
            index_150 = i+len(tm150_keyword)
            indexes_conserved50.append(index_150)
        if oprm_seq[i:i+len(tm250_keyword)] == tm250_keyword:
            index_250 = i+len(tm250_keyword)
            indexes_conserved50.append(index_250)
        if oprm_seq[i:i+len(tm350_keyword)] == tm350_keyword:
            index_350 = i+len(tm350_keyword)
            indexes_conserved50.append(index_350)
        if oprm_seq[i:i+len(tm450_keyword)] == tm450_keyword:
            index_450 = i+len(tm450_keyword)
            indexes_conserved50.append(index_450)
        if oprm_seq[i:i+len(tm550_keyword)] == tm550_keyword:
            index_550 = i+len(tm550_keyword)
            indexes_conserved50.append(index_550)
        if oprm_seq[i:i+len(tm650_keyword)] == tm650_keyword:
            index_650 = i+len(tm650_keyword)
            indexes_conserved50.append(index_650)
        if oprm_seq[i:i+len(tm750_keyword)] == tm750_keyword:
            index_750 = i+len(tm750_keyword)
            indexes_conserved50.append(index_750)
    return indexes_conserved50
            
            
def _find_ahead(ahead, seq, conserved_index):
    while ahead != 0:
        ahead -= 1
        conserved_index -= 1
        while seq[conserved_index] == '-':
            conserved_index -= 1
    return conserved_index
    
def _find_ending(ending, seq, conserved_index):
    while ending != 0:
        ending -= 1
        conserved_index += 1
        while seq[conserved_index] == '-':
            conserved_index += 1
    return conserved_index


def decide_TM_region(seq, indexes_conserved50):
    assert len(indexes_conserved50) == 7
    boundary_list = []
    for t in range(len(indexes_conserved50)):
        boundary_index = boundary_based_on_paper[str(t+1)]
        position_50 = indexes_conserved50[t]
        starting_point = _find_ahead(50 - boundary_index[0], seq, position_50)
        ending_point = _find_ending(boundary_index[1] - 50, seq, position_50)
        boundary_list.append((starting_point, ending_point))
    return boundary_list
 
 
def old_read_fasta(fasta_file):
    f = open(fasta_file,'r')
    lines=f.read().split('\n')
    f.close()
    seq_name_indexes = []
    for i in range(len(lines)):
        if lines[i].startswith('>'):
            seq_name_indexes.append(i)

    sequences_dict = {}
    for s in seq_name_indexes:
        sequences_dict[lines[s]] = lines[s+1]
    return sequences_dict


def read_fasta(fasta_file):
    f = open(fasta_file,'r')
    lines=f.read().split('\n')
    f.close()
    seq_name_indexes = []
    for i in range(len(lines)):
        if lines[i].startswith('>'):
            seq_name_indexes.append(i)

    sequences_dict = {}
    for s in range(len(seq_name_indexes)-1):
        starting_index = seq_name_indexes[s]
        ending_index = seq_name_indexes[s+1]
        sequences_dict[lines[starting_index]] = ''.join( [line for line in lines[starting_index+1:ending_index] if line !=''] )

    starting_index = seq_name_indexes[-1]
    ending_index = -1
    sequences_dict[lines[starting_index]] = ''.join( [line for line in lines[starting_index+1:ending_index] if line !=''] )
    return sequences_dict



#-------------------------------------#
# Runnning starts from below
#-----------------------------------#

sequences_dict = read_fasta(fasta_filename)
oprm_seq = sequences_dict['>OPRM_MOUSE']
if conserved50_list == []:
    conserved50_list = decide_50_position(oprm_seq)
print ('For checking.')
print (conserved50_list)

#----------------------------------------------------------------------#
# This boundary is based on the alignment sequences. May not necessary
#----------------------------------------------------------------------#
boundary = {}
for key in list(sequences_dict.keys()):
    sequence = sequences_dict[key]
    boundary[key] = decide_TM_region(sequence, conserved50_list)


#---------------------------------------------------#
# Based on real sequences, Not aligned, No dash
#---------------------------------------------------#
conserved50_real_positions = {}
real_boundary = {}

for key in list(sequences_dict.keys()):
    conserved50_real_positions[key] = []
    real_boundary[key] = []
    
    sequence = sequences_dict[key]
    no_dash_sequence = []
    empty_positions = []
    for i in range(len(sequence)):
        if sequence[i] == '_' or sequence[i] == '-' :
            empty_positions.append(i)
        else:
            no_dash_sequence.append(sequence[i])
            
    for c in range(len(conserved50_list)):
        empty_length = len([e for e in empty_positions if e < conserved50_list[c]])
        #residue_name = sequence[conserved50_list[c]-1]
        residue_name = no_dash_sequence[conserved50_list[c] - empty_length-1]
        conserved50_real_positions[key].append(residue_name+str(conserved50_list[c] - empty_length))

    for b in range(len(boundary[key])):
        tm_boundary = boundary[key][b]
        starting_point = int(conserved50_real_positions[key][b][1:]) - (50-boundary_based_on_paper[str(b+1)][0])
        #s_residue_name = sequence[tm_boundary[0]-1]
        s_residue_name = no_dash_sequence[starting_point-1]
        ending_point =  int(conserved50_real_positions[key][b][1:]) + (boundary_based_on_paper[str(b+1)][1]-50)
        #e_residue_name = sequence[tm_boundary[1]-1]
        e_residue_name = no_dash_sequence[ending_point-1]
        real_boundary[key].append((s_residue_name+str(starting_point), e_residue_name+str(ending_point)))


f = open('boundary_information.txt','w')
f.write('#Name                TM1           TM2           TM3           TM4           TM5           TM6           TM7\n')
for key in real_boundary.keys():
    tem = [str(item[0])+':'+str(item[1]) for item in real_boundary[key]]
    #f.write("%5s"%key[1:].split('-')[0].upper()+'   '+ ''.join('%14s'%t for t in tem)+'\n')
    #f.write("%11s"%(key[1:].split('-')[0].upper()+'-'+key[1:].split('-')[1])+'   '+ ''.join('%14s'%t for t in tem)+'\n') # For old version
    f.write("%11s"%(key[1:].split('_')[0].upper()+'_'+key[1:].split('_')[1])+'   '+ ''.join('%14s'%t for t in tem)+'\n')
    
tem = []
for key in range(1,8):
    indexes = boundary_based_on_paper[str(key)]
    tem.append('TM'+str(key)+'.'+str(indexes[0])+'-TM'+str(key)+'.'+str(indexes[1]))

f.write('#             '+''.join('%14s'%t for t in tem)+'\n')
f.close()


#print (conserved50_real_positions)
f = open('conserved_50_positions.txt','w')
f.write('#Name\t\tTM1\tTM2\tTM3\tTM4\tTM5\tTM6\tTM7\n')
for key in conserved50_real_positions.keys():
    #f.write(key[1:].split('-')[0].upper()+'\t'+'\t'.join([str(c) for c in conserved50_real_positions[key]])+'\n')
    f.write(key[1:].split('_')[0].upper()+'_'+key[1:].split('_')[1]+'\t'+'\t'.join([str(c) for c in conserved50_real_positions[key]])+'\n')
    #f.write(key[1:].split('-')[0].upper()+'-'+key[1:].split('-')[1]+'\t'+'\t'.join([str(c) for c in conserved50_real_positions[key]])+'\n')
f.close()



