import numpy as np
import os

def merge_AB_metrics(new_matrix,B_new_matrix,filtered_rec_labels,B_filtered_rec_labels,filtered_lig_labels,B_filtered_lig_labels,cutoff):
    union_labels = {}
    for a in filtered_rec_labels:
        union_labels[int(a.split()[0])] = a.split()[1]
    for b in B_filtered_rec_labels:
        union_labels[int(b.split()[0])] = b.split()[1]
    union_labels_keys = sorted(union_labels)
    merged_matrix = np.full((len(union_labels_keys),len(filtered_lig_labels+B_filtered_lig_labels)), cutoff+5, dtype= float)
    filter_y_indexes_Ga = []
    filter_y_indexes_Gb = []

    for r in filtered_rec_labels:
        for i in range(len(union_labels_keys)):
            if int(r.split()[0]) ==union_labels_keys[i]:
                filter_y_indexes_Ga.append(i)
    for r in B_filtered_rec_labels:
        for i in range(len(union_labels_keys)):
            if int(r.split()[0]) ==union_labels_keys[i]:
                filter_y_indexes_Gb.append(i)

    #print (filter_y_indexes_Ga)    
    #print (filter_y_indexes_Gb)
    filter_x_indexes_Ga = [a for a in range(len(filtered_lig_labels))] # G protein
    filter_x_indexes_Gb = [a for a in range(len(filtered_lig_labels), len(filtered_lig_labels)+len(B_filtered_lig_labels))]
    
    for i in range(len(filter_y_indexes_Ga)):
        for j in range(len(filter_x_indexes_Ga)):
            merged_matrix[filter_y_indexes_Ga[i]][filter_x_indexes_Ga[j]] = new_matrix[i][j]

    for i in range(len(filter_y_indexes_Gb)):
        for j in range(len(filter_x_indexes_Gb)):
            merged_matrix[filter_y_indexes_Gb[i]][filter_x_indexes_Gb[j]] = B_new_matrix[i][j]
    union_xlabels = []
    for key in union_labels_keys:
        union_xlabels.append(str(key)+' '+(union_labels[key]))        
    return merged_matrix, union_xlabels, filtered_lig_labels+B_filtered_lig_labels


def merge_AB_metrics_freq(new_matrix,B_new_matrix,filtered_rec_labels,B_filtered_rec_labels,filtered_lig_labels,B_filtered_lig_labels,cutoff):
    union_labels = {}
    for a in filtered_rec_labels:
        union_labels[int(a.split()[0])] = a.split()[1]
    for b in B_filtered_rec_labels:
        union_labels[int(b.split()[0])] = b.split()[1]
    union_labels_keys = sorted(union_labels)
    merged_matrix = np.zeros((len(union_labels_keys),len(filtered_lig_labels+B_filtered_lig_labels)))
    filter_y_indexes_Ga = []
    filter_y_indexes_Gb = []

    for r in filtered_rec_labels:
        for i in range(len(union_labels_keys)):
            if int(r.split()[0]) ==union_labels_keys[i]:
                filter_y_indexes_Ga.append(i)
    for r in B_filtered_rec_labels:
        for i in range(len(union_labels_keys)):
            if int(r.split()[0]) ==union_labels_keys[i]:
                filter_y_indexes_Gb.append(i)

    #print (filter_y_indexes_Ga)    
    #print (filter_y_indexes_Gb)
    filter_x_indexes_Ga = [a for a in range(len(filtered_lig_labels))] # G protein
    filter_x_indexes_Gb = [a for a in range(len(filtered_lig_labels), len(filtered_lig_labels)+len(B_filtered_lig_labels))]

    for i in range(len(filter_y_indexes_Ga)):
        for j in range(len(filter_x_indexes_Ga)):
            merged_matrix[filter_y_indexes_Ga[i]][filter_x_indexes_Ga[j]] = new_matrix[i][j]

    for i in range(len(filter_y_indexes_Gb)):
        for j in range(len(filter_x_indexes_Gb)):
            merged_matrix[filter_y_indexes_Gb[i]][filter_x_indexes_Gb[j]] = B_new_matrix[i][j]
    union_xlabels = []
    for key in union_labels_keys:
        union_xlabels.append(str(key)+' '+(union_labels[key]))
    return merged_matrix, union_xlabels, filtered_lig_labels+B_filtered_lig_labels


def make_union_list(label_list1, label_list2):
    union_list = list(set(label_list1).union(set(label_list2)))
    turpled_union_list = [] # sorting list
    for u in union_list:
        if not 'Gb' in u:
            turpled_union_list.append((int(u.split()[0]), u.split()[1]))
        else:
            turpled_union_list.append((int(u.split()[1]), u.split()[2]+' Gb'))
    turpled_union_list = sorted(turpled_union_list)
    sorted_union_list = []
    for t in turpled_union_list:
        if not 'Gb' in t[1]:
            sorted_union_list.append(str(t[0])+' '+t[1])
        else:
            sorted_union_list.append('Gb '+str(t[0])+' '+t[1].split()[0])
    return sorted_union_list

def expanded_matrix(matrix, union_xlables, union_ylables, xlabels, ylabels):
    """
    expand matrix size to an appointed size. Fill in the blank position with zero.
    Used for union frequency matrix.

    """
    expanded_fen_matrix = np.zeros((len(union_xlables), len(union_ylables)))

    need_added_fen_x_indexes = []
    existed_fen_xindexes = []

    for i in range(len(union_xlables)):
        if not union_xlables[i] in xlabels:
            need_added_fen_x_indexes.append(i)
        else:
            existed_fen_xindexes.append(i)

    need_added_fen_y_indexes = []
    existed_fen_yindexes = []

    for i in range(len(union_ylables)):
        if not union_ylables[i] in ylabels:
            need_added_fen_y_indexes.append(i)
        else:
            existed_fen_yindexes.append(i)        

    for i in range(len(existed_fen_xindexes)):
        for j in range(len(existed_fen_yindexes)):
            new_xindex = existed_fen_xindexes[i]
            new_yindex = existed_fen_yindexes[j]
            expanded_fen_matrix[new_xindex][new_yindex] = matrix[i][j]
            
    return  expanded_fen_matrix

def expanded_matrix_dist_type(matrix, union_xlables, union_ylables, xlabels, ylabels, dist_cutoff):
    """
    expand matrix size to an appointed size. Fill in the blank position with dist_cutoff.
    Used for union distance matrix.

    """
    expanded_fen_matrix = np.full((len(union_xlables), len(union_ylables)),float(dist_cutoff))

    need_added_fen_x_indexes = []
    existed_fen_xindexes = []

    for i in range(len(union_xlables)):
        if not union_xlables[i] in xlabels:
            need_added_fen_x_indexes.append(i)
        else:
            existed_fen_xindexes.append(i)

    need_added_fen_y_indexes = []
    existed_fen_yindexes = []

    for i in range(len(union_ylables)):
        if not union_ylables[i] in ylabels:
            need_added_fen_y_indexes.append(i)
        else:
            existed_fen_yindexes.append(i)

    for i in range(len(existed_fen_xindexes)):
        for j in range(len(existed_fen_yindexes)):
            new_xindex = existed_fen_xindexes[i]
            new_yindex = existed_fen_yindexes[j]
            expanded_fen_matrix[new_xindex][new_yindex] = matrix[i][j]

    return  expanded_fen_matrix



def _look_up_pixel_from_matrix(matrix, look_up_xlabel, look_up_ylabel, xlabels, ylabels):
    """
    given an certain label, search the value from the raw matrix.
    """    
    selected_x = None
    selected_y = None
    for i in range(len(xlabels)):
        if xlabels[i] == look_up_xlabel:
            selected_x = i
    for i in range(len(ylabels)):
        if ylabels[i] == look_up_ylabel:
            selected_y = i
    value = matrix[selected_x][selected_y]
    return value


def expand_matrix_dist_by_added_pixel(matrix, union_xlables, union_ylables, xlabels, ylabels, raw_matrix, raw_xlabels, raw_ylabels, dist_cutoff):
    expanded_fen_matrix = np.zeros((len(union_xlables), len(union_ylables)))

    need_added_fen_x_indexes = []
    existed_fen_xindexes = []
    need_added_x_labels = []

    for i in range(len(union_xlables)):
        if not union_xlables[i] in xlabels:
            need_added_fen_x_indexes.append(i)
            need_added_x_labels.append(union_xlables[i])
        else:
            existed_fen_xindexes.append(i)

    need_added_fen_y_indexes = []
    existed_fen_yindexes = []
    need_added_y_labels = []

    for i in range(len(union_ylables)):
        if not union_ylables[i] in ylabels:
            need_added_fen_y_indexes.append(i)
            need_added_y_labels.append(union_ylables[i])
        else:
            existed_fen_yindexes.append(i)
 

    for i in range(len(existed_fen_xindexes)):
        for j in range(len(existed_fen_yindexes)):

            new_xindex = existed_fen_xindexes[i]
            new_yindex = existed_fen_yindexes[j]


            expanded_fen_matrix[new_xindex][new_yindex] = matrix[i][j]

    #for i in range(len(need_added_x_labels)):
    #    for j in range(len(need_added_y_labels)):
           

    return  expanded_fen_matrix




def shrink_matrix_with_extra_boundary(matrix, xlabel, ylabel, shrink_cutoff, para_cutoff):
    """
    so far this is used for the distance matrix. The shrinked matrix includes cutoff boundary.
    """
    not_empty_pairs = [] # based on the config parameter
    shrinked_not_empty_pairs = []
    
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i][j] < para_cutoff+5:
                if matrix[i][j] <= shrink_cutoff:
                    shrinked_not_empty_pairs.append((i,j))
                else:
                    not_empty_pairs.append((i,j))
    shrink_not_empty_x = []
    shrink_not_empty_y = []
    for s in shrinked_not_empty_pairs:
        shrink_not_empty_x.append(s[0])
        shrink_not_empty_y.append(s[1])
    shrink_not_empty_x = sorted(list(set(shrink_not_empty_x)))
    shrink_not_empty_y = sorted(list(set(shrink_not_empty_y)))

    filtered_matrix = np.full((len(shrink_not_empty_x), len(shrink_not_empty_y)), float(para_cutoff+5))
    for i in range(len(shrink_not_empty_x)):
        for j in range(len(shrink_not_empty_y)):
            filtered_matrix[i][j] = matrix[shrink_not_empty_x[i]][shrink_not_empty_y[j]]
          
    new_xlabel = [xlabel[n] for n in shrink_not_empty_x]
    new_ylabel = [ylabel[n] for n in shrink_not_empty_y]

    return filtered_matrix, new_xlabel, new_ylabel
   

def shrink_matrix(matrix, xlabel, ylabel, shrink_cutoff):
    """
    so far this is used for the distance matrix. The shrinked matrix includes cutoff boundary.
    """
    not_empty_pairs = [] # based on the config parameter
    shrinked_not_empty_pairs = []

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i][j] <= shrink_cutoff:
                shrinked_not_empty_pairs.append((i,j))
    shrink_not_empty_x = []
    shrink_not_empty_y = []
    for s in shrinked_not_empty_pairs:
        shrink_not_empty_x.append(s[0])
        shrink_not_empty_y.append(s[1])
    shrink_not_empty_x = sorted(list(set(shrink_not_empty_x)))
    shrink_not_empty_y = sorted(list(set(shrink_not_empty_y)))
    filtered_matrix = np.zeros((len(shrink_not_empty_x), len(shrink_not_empty_y)))
    for i in range(len(shrink_not_empty_x)):
        for j in range(len(shrink_not_empty_y)):
            filtered_matrix[i][j] = matrix[shrink_not_empty_x[i]][shrink_not_empty_y[j]]
    new_xlabel = [xlabel[n] for n in shrink_not_empty_x]
    new_ylabel = [ylabel[n] for n in shrink_not_empty_y]
    return filtered_matrix, new_xlabel, new_ylabel




def clean_matrix(matrix, xlabel, ylabel, cutoff, mark_out_range_pair= []):
    new_matrix = np.zeros((len(xlabel), len(ylabel)))
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if mark_out_range_pair != []:
                if not (i,j) in mark_out_range_pair:
                    if np.abs(matrix[i][j]) > cutoff:
                        new_matrix[i][j] = matrix[i][j]
            else:
                if np.abs(matrix[i][j]) > cutoff:
                    new_matrix[i][j] = matrix[i][j]

    return new_matrix, xlabel, ylabel
     

def find_different_pair(distance,xlabel, ylabel,threshold):
    """
    report difference between two different matrix
    input: distance is a matrix that is obtained by using one matrix to subtract the other matrix

    """
    
    tem = []
    for i in range(distance.shape[0]):
        for j in range(distance.shape[1]):
            if distance[i][j]>threshold or distance[i][j] < -1.0*threshold:
                tem.append((i,j))
    different_interaction_pairs = []
    for t in tem:
        # Make sure only index number was included into account.
        if len(xlabel[t[0]].split())>1:
            x_num = None
            for x in xlabel[t[0]].split():
                if x.isdigit():
                    x_num = x 
        else:
            x_num = xlabel[t[0]]
        if len(ylabel[t[1]].split())>1:
           y_num = None
           for y in ylabel[t[1]].split():
               if y.isdigit():
                   y_num = y
        else:
           y_num = ylabel[t[1]]
        different_interaction_pairs.append((x_num,y_num,distance[t[0], t[1]]))
        #different_interaction_pairs.append((xlabel[t[0]],ylabel[t[1]], distance[t[0], t[1]]))
    return different_interaction_pairs

def find_different_pair_distance_type(distance,xlabel, ylabel,threshold):
    """
    report difference between two different matrix
    input: distance is a matrix that is obtained by using one matrix to subtract the other matrix

    """

    tem = []
    for i in range(distance.shape[0]):
        for j in range(distance.shape[1]):
            if distance[i][j]<threshold: 
                tem.append((i,j))
    different_interaction_pairs = []
    for t in tem:
        # Make sure only index number was included into account.
        if len(xlabel[t[0]].split())>1:
            x_num = None
            for x in xlabel[t[0]].split():
                if x.isdigit():
                    x_num = x
        else:
            x_num = xlabel[t[0]]
        if len(ylabel[t[1]].split())>1:
           y_num = None
           for y in ylabel[t[1]].split():
               if y.isdigit():
                   y_num = y
        else:
           y_num = ylabel[t[1]]
        different_interaction_pairs.append((x_num,y_num,distance[t[0], t[1]]))
        #different_interaction_pairs.append((xlabel[t[0]],ylabel[t[1]], distance[t[0], t[1]]))
    return different_interaction_pairs




def write_down_difference(diff_pair, outputfile):
    f = open(outputfile,'w')
    receptor_list = []
    ga_list = []
    values = []
    for d in diff_pair:
        receptor_list.append(d[0])
        ga_list.append(d[1])
        values.append(d[2])
    f.write("receptor_list = ["+"'"+"','".join(str(x.split()[0]) for x in receptor_list)+"']\n")
    f.write("galpha_list = ["+"'"+"','".join(str(x.split()[0]) for x in ga_list)+"']\n")
    f.write("radius_list = [" + ",".join(str(round(x,2)) for x in values)+"]\n")
    f.close()
    return

def filter_empty_matrix(matrix, xlabel, ylabel):
    empty_x = []
    empty_y = []
    not_empty_x = []
    not_empty_y = []

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i][j] == 0:
                empty_x.append(i)
                empty_y.append(j)
            else:
                not_empty_x.append(i)
                not_empty_y.append(j)
    filtered_matrix = np.zeros((len(not_empty_x), len(not_empty_y)))
    for i in range(len(not_empty_x)):
        for j in range(len(not_empty_y)):
            filtered_matrix[i][j] = matrix[not_empty_x[i]][not_empty_y[j]]
    new_xlabel = [xlabel[n] for n in not_empty_x]
    new_ylabel = [ylabel[n] for n in not_empty_y]

    return filtered_matrix, new_xlabel, new_ylabel


def filter_duplicate_matrix(matrix,xlabel, ylabel):
    non_repeatable_xlabel = list(set(xlabel))
    non_repeatable_ylabel = list(set(ylabel))

    tem_x = []
    for n in non_repeatable_xlabel:
        if not 'Gb' in n:
            tem_x.append((int(n.split()[0]), n.split()[1]))
        else:
            tem_x.append((int(n.split()[1]), n.split()[2]+' Gb'))
    tem_x = sorted(tem_x)
    sorted_xlabel = []
    for t in tem_x:
        if not 'Gb' in t[1]:
             sorted_xlabel.append(str(t[0])+' '+t[1])
        else:
             sorted_xlabel.append('Gb '+str(t[0])+' '+t[1].split()[0])

    tem_y = []
    for n in non_repeatable_ylabel:
        if not 'Gb' in n:
            tem_y.append((int(n.split()[0]), n.split()[1]))
        else:
            tem_y.append((int(n.split()[1]), n.split()[2]+' Gb'))
    tem_y = sorted(tem_y)
    sorted_ylabel = []
    for t in tem_y:
        if not 'Gb' in t[1]:
            sorted_ylabel.append(str(t[0])+' '+t[1])
        else:
            sorted_ylabel.append('Gb '+str(t[0])+' '+t[1].split()[0])        
        
    non_rep_xindexes = []
    non_rep_yindexes = []
    for nx in sorted_xlabel: 
        for i in range(len(xlabel)):
            if nx == xlabel[i]:
                non_rep_xindexes.append(i)
                break
    for ny in sorted_ylabel:
        for i in range(len(ylabel)):
            if ny == ylabel[i]:
                non_rep_yindexes.append(i)
                break

    new_matrix = np.zeros((len(non_rep_xindexes),  len(non_rep_yindexes)))
    for i in range(len(non_rep_xindexes)):
        for j in range(len(non_rep_yindexes)):
            new_matrix[i][j] = matrix[non_rep_xindexes[i]][non_rep_yindexes[j]]
    return new_matrix, sorted_xlabel, sorted_ylabel


def reorder_matrix(matrix, xlabel, ylabel):
    reorder_xlabel_index = []
    reorder_ylabel_index = []
    for x in range(len(xlabel)): ### reorder G protein index, so that Ga would appear first, and Gb as second.
        if not 'Gb' in xlabel[x]:
            reorder_xlabel_index.append(x)
    for x in range(len(xlabel)):
        if 'Gb' in  xlabel[x]:
            reorder_xlabel_index.append(x)

    for y in range(len(ylabel)): ### reorder G protein index, so that Ga would appear first, and Gb as second.
        if not 'Gb' in ylabel[y]:
            reorder_ylabel_index.append(y)
    for y in range(len(ylabel)):
        if 'Gb' in  ylabel[y]:
            reorder_ylabel_index.append(y)

    new_matrix = np.zeros((len(xlabel), len(ylabel)))
    for i in range(len(reorder_xlabel_index)):
        for j in range(len(reorder_ylabel_index)):
            new_matrix[i][j] = matrix[reorder_xlabel_index[i]][reorder_ylabel_index[j]]
    return new_matrix, [xlabel[r] for r in reorder_xlabel_index], [ylabel[r] for r in reorder_ylabel_index]

