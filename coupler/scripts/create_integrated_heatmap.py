"""
This script is based on the results of pdb distance matrics, and try to summary the similarities and differences.
"""
import numpy as np 

def pick_non_empty_pairs(matrix, xlabels, ylabels, cutoff):
    non_empty_pair = []
    assert len(xlabels) == matrix.shape[0]
    assert len(ylabels) == matrix.shape[1]
    for i in range(len(xlabels)):
        for j in range(len(ylabels)):
            if matrix[i][j] < cutoff:
                non_empty_pair.append((xlabels[i], ylabels[j]))
    return non_empty_pair


def pick_non_empty_pairs_boolean(matrix, xlabels, ylabels):
    non_empty_pair = []
    assert len(xlabels) == matrix.shape[0]
    assert len(ylabels) == matrix.shape[1]
    for i in range(len(xlabels)):
        for j in range(len(ylabels)):
            if matrix[i][j] != 0:
                non_empty_pair.append((xlabels[i], ylabels[j]))
    return non_empty_pair


def _make_union(pair_lists):
    """
    pair_lists = [pair_list1, pair_list2, pair_list3, ...]
    pair_list* should be a list filled with tuple.
    """
    union_list = []
    for i in range(len(pair_lists)):
        pair_list = pair_lists[i]
        for p in pair_list:
            if not p in union_list:
                union_list.append(p)
    return union_list

def _intersection(lst1, lst2):
    """
    Find the common set between two lists
    """
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3

def _split_label(label):
    """
    split the label to number(index), TM/ICL/ECL/Gb, resname
    """
    if 'TM' in label:
        numstr = float(label.split()[0][2:])
        resname = label.split()[1]
        tag = 'TM'
    elif 'ICL' in label:
        numstr = float(label.split()[0][3:])
        tag = 'ICL'
        resname = label.split()[1]
    elif 'ECL' in label:
        numstr = float(label.split()[0][3:])
        tag = 'ECL'
        resname = label.split()[1]
    elif 'Gb' in label:
        numstr = int(label.split()[1])
        tag = 'Gb'
        resname = label.split()[2]
    else:
        numstr = int(label.split()[0])
        tag = ''
        resname = label.split()[1]
    return numstr, tag, resname
   
def _split_label_full(label):
    """
    split the label to number(index), TM/ICL/ECL/Gb, resname
    """
    if 'TM' in label:
        numstr = float(label.split()[0][2:])
        resname = label.split()[1]
        tag = 'TM'
        real_resnum = label.split()[2]

    elif 'ICL' in label:
        numstr = float(label.split()[0][3:])
        tag = 'ICL'
        resname = label.split()[1]
        real_resnum = label.split()[2]

    elif 'ECL' in label:
        numstr = float(label.split()[0][3:])
        tag = 'ECL'
        resname = label.split()[1]
        real_resnum = label.split()[2]

    elif 'Gb' in label:
        numstr = int(label.split()[1])
        tag = 'Gb'
        resname = label.split()[2]
        real_resnum = label.split()[1]

    else:
        numstr = int(label.split()[0])
        tag = ''
        resname = label.split()[1]
        real_resnum = label.split()[2]

    return numstr, tag, resname, real_resnum


 
def _keep_2_digits(single_label):
    """
    This is only used for the receptor TM index.
    Since after converting the index to float, it will probably lose 1 digit. e.g. TM1.60 -->TM1.6
    This function makes sure that the TM index has two digits.
    """
    new_label = single_label
    if 'TM' in single_label:
        new_label = 'TM'+"%3.2f"%float(single_label.split()[0][2:])+' '+single_label.split()[1]
    return new_label

def _keep_2_digits_full(single_label):
    """
    This is only used for the receptor TM index.
    Since after converting the index to float, it will probably lose 1 digit. e.g. TM1.60 -->TM1.6
    This function makes sure that the TM index has two digits.
    """
    #new_label = single_label
    if 'TM' in single_label:
        new_label = 'TM'+"%4.2f"%float(single_label.split()[0][2:])+' '+single_label.split()[1]+' '+single_label.split()[2]
    elif 'ICL' in single_label:
        if len(single_label.split()[0][3:]) == 5: #ICL1.150
            new_label = 'ICL'+"%5.3f"%float(single_label.split()[0][3:])+' '+single_label.split()[1]+' '+single_label.split()[2]
        else: #ICL1.62
            new_label = 'ICL'+"%4.2f"%float(single_label.split()[0][3:])+' '+single_label.split()[1]+' '+single_label.split()[2]
    elif 'ECL' in single_label:
        new_label = 'ECL'+"%5.3f"%float(single_label.split()[0][3:])+' '+single_label.split()[1]+' '+single_label.split()[2]
    else:
        new_label = single_label
    return new_label



def _sort_list(label_list):
    tuple_list = []
    for l in label_list:
        numstr, tag, resname = _split_label(l)
        tuple_list.append((numstr, tag+' '+resname))
    tuple_list = sorted(tuple_list)
    sorted_label_list = []
    for t in tuple_list:
        if len(t[1].split()) == 2:
            sorted_label_list.append(t[1].split()[0]+str(t[0])+' '+t[1].split()[1])
        else:
            sorted_label_list.append(str(t[0])+' '+t[1].split()[0])

    resorted_label_list = []
    label_type = 'ga'
    
    for s in sorted_label_list:
        if 'Gb' in s:
            label_type = 'gb'
            break
    for s in sorted_label_list:
        if 'TM' in s:
            label_type = 'receptor'
            break
        
    if label_type == 'gb':
        for s in sorted_label_list:
            if not 'Gb' in s:
                resorted_label_list.append(s)
        for s in sorted_label_list:
            if 'Gb' in s:
                resorted_label_list.append('Gb '+s[2:])
                
    elif label_type == 'receptor':
        for s in sorted_label_list:
            if 'TM1' in s:
                resorted_label_list.append(_keep_2_digits(s))
        for s in sorted_label_list:
            if 'ICL1' in s:
                resorted_label_list.append(s)
        for s in sorted_label_list:
            if 'TM2' in s:
                resorted_label_list.append(_keep_2_digits(s))
        for s in sorted_label_list:
            if 'ECL1' in s:
                resorted_label_list.append(s)
        for s in sorted_label_list:
            if 'TM3' in s:
                resorted_label_list.append(_keep_2_digits(s))
        for s in sorted_label_list:
            if 'ICL2' in s:
                resorted_label_list.append(s)
        for s in sorted_label_list:
            if 'TM4' in s:
                resorted_label_list.append(_keep_2_digits(s))
        for s in sorted_label_list:
            if 'ECL2' in s:
                resorted_label_list.append(s)
        for s in sorted_label_list:
            if 'TM5' in s:
                resorted_label_list.append(_keep_2_digits(s))
        for s in sorted_label_list:
            if 'ICL3' in s:
                resorted_label_list.append(s)
        for s in sorted_label_list:
            if 'TM6' in s:
                resorted_label_list.append(_keep_2_digits(s))
        for s in sorted_label_list:
            if 'ECL3' in s:
                resorted_label_list.append(s)
        for s in sorted_label_list:
            if 'TM7' in s:
                resorted_label_list.append(_keep_2_digits(s))
        for s in sorted_label_list:
            if 'TM8' in s:
                resorted_label_list.append(_keep_2_digits(s))
        for s in sorted_label_list:
            if not 'TM' in s and not 'CL' in s:
                resorted_label_list.append(s)
    else:
        for s in sorted_label_list:
            resorted_label_list.append(s)

    return resorted_label_list

def _sorted_list_full(label_list):
    tuple_list = []
    for l in label_list:
        numstr, tag, resname, real_resnum = _split_label_full(l)
        tuple_list.append((numstr, tag+' '+resname+' '+real_resnum))
    tuple_list = sorted(tuple_list)
    sorted_label_list = []
    for t in tuple_list:
        if len(t[1].split()) == 3:
            sorted_label_list.append(t[1].split()[0]+str(t[0])+' '+t[1].split()[1]+' '+t[1].split()[2])
        else:
            sorted_label_list.append(str(t[0])+' '+t[1].split()[0]+' '+t[1].split()[1])

    resorted_label_list = []
    label_type = 'ga'
    
    for s in sorted_label_list:
        if 'Gb' in s:
            label_type = 'gb'
            break
    for s in sorted_label_list:
        if 'TM' in s:
            label_type = 'receptor'
            break
        
    if label_type == 'gb':
        for s in sorted_label_list:
            if not 'Gb' in s:
                resorted_label_list.append(s)
        for s in sorted_label_list:
            if 'Gb' in s:
                resorted_label_list.append('Gb '+s[2:].split()[0]+' '+s[2:].split()[1])
                
    elif label_type == 'receptor':
        for s in sorted_label_list:
            if 'TM1' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'ICL1' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'TM2' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'ECL1' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'TM3' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'ICL2' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'TM4' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'ECL2' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'TM5' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'ICL3' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'TM6' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'ECL3' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'TM7' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if 'TM8' in s:
                resorted_label_list.append(_keep_2_digits_full(s))
        for s in sorted_label_list:
            if not 'TM' in s and not 'CL' in s:
                resorted_label_list.append(s)
    else:
        for s in sorted_label_list:
            resorted_label_list.append(s)

    return resorted_label_list



def build_matrix(pair_lists):
    union_list = _make_union(pair_lists)
    xlist = [u[0] for u in union_list]
    ylist = [u[1] for u in union_list]
    #xlist = _sort_list(list(set(xlist)))
    #ylist = _sort_list(list(set(ylist)))
    xlist = _sorted_list_full(list(set(xlist)))
    ylist = _sorted_list_full(list(set(ylist)))
    matrix = np.zeros((len(xlist), len(ylist)))
    indexes = {}

    for i in range(len(pair_lists)):
        pair_list = pair_lists[i] 
        indexes[i] = []   
        for p in pair_list:
            for j in range(len(xlist)):
                if xlist[j] == p[0]:
                    for m in range(len(ylist)):
                        #if ylist[m] == p[1]:
                        if not 'Gb' in ylist[m]:
                            if ylist[m].split()[0] == p[1].split()[0]:
                                indexes[i].append((j,m))
                        else:
                            if ylist[m].split()[1] == p[1].split()[1]:
                                indexes[i].append((j,m))
    for i in range(len(indexes)): ### loop systems
        for index in indexes[i]: ### loop each system pair
            x = index[0]
            y = index[1]
            matrix[x][y] +=1
    return matrix, xlist, ylist


def build_short_name_matrix(pair_lists):
    union_list = _make_union(pair_lists)
    xlist = [u[0] for u in union_list]
    ylist = [u[1] for u in union_list]
    xlist = sorted(list(set(xlist))) 
    ylist = sorted(list(set([int(y) for y in ylist])))

    matrix = np.zeros((len(xlist), len(ylist)))
    indexes = {}

    for i in range(len(pair_lists)):
        pair_list = pair_lists[i]
        indexes[i] = []
        for p in pair_list:
            for j in range(len(xlist)):
                if xlist[j] == p[0]:
                    for m in range(len(ylist)):
                        if str(ylist[m]) == p[1]:
                            indexes[i].append((j,m))

    for i in range(len(indexes)): ### loop systems
        for index in indexes[i]: ### loop each system pair
            x = index[0]
            y = index[1]
            matrix[x][y] +=1
    return matrix, xlist, [str(y) for y in ylist]
           
