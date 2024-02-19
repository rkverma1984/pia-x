#!/usr/bin/env python
"""
VMD-Python calculate contact frequency
Created by Bing, 09/15/2020
Modified on 09/16/2020, add distance calculation function
Modified on 10/27/2020, add wrap up function
Modified on 12/8/2020, move the "if" statement "calc_type" out of the loop;
                       change label from "%3d" to "%d", sort label by lambda calculus;
                       modify the denominator "num_steps" of frequency calculation to "num_steps - 1"
"""
import os
import numpy as np
import residue_hash_table
from tqdm import tqdm
try:
    from vmd import atomsel
    from vmd import molecule
except:
    print('vmd function is not loaded. can only read checkpoint files.')

residue_names_dic = residue_hash_table.residue_name_3to1()


# print (residue_names_dic)

def _chain2resid(chainID):
    resids = atomsel( "not hydrogen and chain " + str( chainID ) ).resid
    tem_resids = list( set( resids ) )
    resids = []
    for t in tem_resids:
        if not 'ACE' in _resid2resname(t, chainID) and not 'NMA' in  _resid2resname(t, chainID):
            resids.append(t)
    return resids

def _resid2resname(index, chainID):
    residue_name = atomsel( 'resid ' + str( index ) + ' and chain ' + chainID ).resname
    # one_letter_resname = residue_names_dic[str(residue_name)]
    return residue_name


def _chain2resid_for_receptor(chainID):
    resids = atomsel( "protein and not hydrogen and chain " + str( chainID ) ).resid
    tem_resids = list( set( resids ) )
    resids = []
    for t in tem_resids:
        if not 'ACE' in _resid2resname(t, chainID) and not 'NMA' in  _resid2resname(t, chainID):
            resids.append(t)
    return resids


def _index2resid(index, chainID):
    residue_id = atomsel( 'index ' + str( index ) + ' and chain ' + chainID ).resid
    return residue_id


def _index2resname(index, chainID):
    residue_name = atomsel( 'index ' + str( index ) + ' and chain ' + chainID ).resname
    return residue_name



class Frequency:
    def __init__(self, pdb_file, dcd_file, receptor_chain_ID, ligand_chain_ID, calc_type, cutoff=5,
                 frequency_threshold=0.005, stride=100, bootstrap=False, rec_resids=[], lig_resids=[]):
        self._pdb_file = pdb_file
        self._dcd_file = dcd_file
        self._rec_resids = rec_resids
        self._lig_resids = lig_resids
        self._receptor_chain_ID = receptor_chain_ID
        self._ligand_chain_ID = ligand_chain_ID
        self._cutoff = cutoff
        self._stride = stride
        self._frequency_threshold = frequency_threshold
        self._bootstrap = bootstrap
        self._calc_type = calc_type

    @property
    def contact_freq_vmd(self):
        """
        pdb_file: top pdb
        dcd_file: trajectory file, dcd format
        rec_resids: receptor residues to be calculated for contact, list type
        lig_resids: ligand residues to be calculated for contact, list type
        """
        top_pdb = molecule.load( 'pdb', self._pdb_file )
        molecule.read( top_pdb, "dcd", self._dcd_file, stride=self._stride, waitfor=-1 )
        num_steps = molecule.numframes( top_pdb )  # pdb frame + dcd frame
        receptor_labels = []
        ligand_labels = []

        if self._lig_resids == []:  # ligand residues not defined, calculate all residues in that chain
            lig_resids = _chain2resid( self._ligand_chain_ID )
        else:
            lig_resids = [int(s) for s in self._lig_resids]

        if self._rec_resids == []:  # ligand residues not defined, calculate all residues in that chain
            print( 'It will search all possible residues and may take a while.' )
            rec_resids = _chain2resid_for_receptor( self._receptor_chain_ID )
        else:
            print( 'range ' + str( self._rec_resids[0] ) + ' to ' + str(
                self._rec_resids[-1] ) + ' residues are used for calculating contact frequency.' )
            rec_resids = [int(s) for s in self._rec_resids]
 
        if self._calc_type == 'heavy_atom':
            keyword = ' and not hydrogen'
        elif self._calc_type == 'backbone':
            keyword = ' and backbone'
        elif self._calc_type == 'CA':
            keyword = ' and name CA'
        else:
            print( 'Unknown calculation type. Will take all atoms into account.' )
            keyword = ''

        frequency_matrix = np.zeros( (len( rec_resids ), len( lig_resids )) )
        not_empty_xindex = []
        not_empty_yindex = []
        if self._bootstrap:
            print( 'Boostraping trajectory...' )
        for r in tqdm( range( len( rec_resids ) ) ):  # control by residue index
            r_id = rec_resids[r]
            for l in range( len( lig_resids ) ):
                l_id = lig_resids[l]

                receptor = atomsel( "resid " + str( r_id ) + keyword + " and chain " + self._receptor_chain_ID )
                ligand = atomsel( "resid " + str( l_id ) + keyword + " and chain " + self._ligand_chain_ID )

                if not self._bootstrap:
                    count = 0
                    for t in range( 1, num_steps ):
                        receptor.frame = t
                        ligand.frame = t
                        receptor.update()
                        ligand.update()

                        if receptor.contacts( ligand, self._cutoff )[0] != []:
                            count += 1
                            
                    freq = float( count ) / float( num_steps - 1 )
                else:
                    freq = 0
                    for _ in range( 100 ):  # repeat 100 times
                        tem_count = 0
                        random_steps = np.random.choice( num_steps - 1, num_steps - 1 )
                        for t in random_steps:
                            receptor.frame = int( t )
                            ligand.frame = int( t )
                            receptor.update()
                            ligand.update()
                            if receptor.contacts( ligand, self._cutoff )[0] != []:
                                tem_count += 1
                        tem_freq = float( tem_count ) / float( len( random_steps ) )
                        freq += tem_freq
                    freq = freq / float( 100.0 )
                if freq > self._frequency_threshold:
                    # print (freq)
                    frequency_matrix[r][l] = freq
                    # print ('keep record r_id', r_id)
                    rec_residue = _resid2resname( r_id, self._receptor_chain_ID )[0]
                    lig_residue = _resid2resname( l_id, self._ligand_chain_ID )[0]

                    receptor_labels.append( "%d" % (r_id) + ' ' + residue_names_dic[rec_residue] )
                    if lig_residue in residue_names_dic.keys():
                        ligand_labels.append( "%d" % (l_id) + ' ' + residue_names_dic[lig_residue] )
                    else:
                        ligand_labels.append( "%d" % (l_id) + ' ' + lig_residue )
                    not_empty_xindex.append( r )
                    not_empty_yindex.append( l )
        not_empty_xindex = sorted( list( set( not_empty_xindex ) ) )
        not_empty_yindex = sorted( list( set( not_empty_yindex ) ) )
        filtered_frequency_matrix = np.zeros( (len( not_empty_xindex ), len( not_empty_yindex )) )
        for i in range( len( not_empty_xindex ) ):
            for j in range( len( not_empty_yindex ) ):
                filtered_frequency_matrix[i][j] = frequency_matrix[not_empty_xindex[i]][not_empty_yindex[j]]
        receptor_labels = sorted( list( set( receptor_labels ) ), key = lambda x: (len (x), x) )
        ligand_labels = sorted( list( set( ligand_labels ) ), key = lambda x: (len (x), x) )
        return filtered_frequency_matrix, receptor_labels, ligand_labels


def calc_dist(sel1, sel2):
    sel1_coords = []
    for s in range( len( sel1.x ) ):
        sel1_coords.append( [sel1.x[s], sel1.y[s], sel1.z[s]] )
    sel1_coords = np.array( sel1_coords )

    sel2_coords = []
    for s in range( len( sel2.x ) ):
        sel2_coords.append( [sel2.x[s], sel2.y[s], sel2.z[s]] )
    sel2_coords = np.array( sel2_coords )

    distance = []
    for i in range( len( sel1_coords ) ):
        for j in range( len( sel2_coords ) ):
            dist = np.sqrt(
                (sel1_coords[i][0] - sel2_coords[j][0]) ** 2 + (sel1_coords[i][1] - sel2_coords[j][1]) ** 2 + (
                            sel1_coords[i][2] - sel2_coords[j][2]) ** 2 )
            distance.append( dist )
    try:
        return round( np.min( distance ), 2 )
    except: 
        print (sel1, sel2, 'has something wrong.')

class Distance:
    def __init__(self, pdb_file, receptor_chain_ID, ligand_chain_ID, cutoff, calc_type, dcd_file=None, stride=100,
                 rec_resids=[], lig_resids=[]):
        self._pdb_file = pdb_file
        self._dcd_file = dcd_file
        self._rec_resids = rec_resids
        self._lig_resids = lig_resids
        self._receptor_chain_ID = receptor_chain_ID
        self._ligand_chain_ID = ligand_chain_ID
        self._stride = stride
        self._cutoff = cutoff
        self._calc_type = calc_type

    def contact_dist_vmd(self):
        """
        pdb_file: top pdb
        dcd_file: trajectory file, dcd format
        rec_resids: receptor residues to be calculated for contact, list type
        lig_resids: ligand residues to be calculated for contact, list type
        """
        receptor_labels = []
        ligand_labels = []
        top_pdb = molecule.load( 'pdb', self._pdb_file )

        if self._lig_resids == []:  # ligand residues not defined, calculate all residues in that chain
            lig_resids = _chain2resid( self._ligand_chain_ID )
        else:
            lig_resids = [int(s) for s in self._lig_resids]
        if self._rec_resids == []:  # ligand residues not defined, calculate all residues in that chain
            print( 'It will search all possible residues and may take a while.' )
            rec_resids = _chain2resid_for_receptor( self._receptor_chain_ID )
        else:
            print( 'range ' + str( self._rec_resids[0] ) + ' to ' + str(
                self._rec_resids[-1] ) + ' residues are used for calculating contact frequency.' )
            rec_resids = [int(s) for s in self._rec_resids]

        if self._dcd_file != None:
            # top_pdb = molecule.load('pdb',self._pdb_file)
            molecule.read( top_pdb, "dcd", self._dcd_file, stride=self._stride, waitfor=-1 )
            num_steps = molecule.numframes( top_pdb )  # pdb frame + dcd frame
            
            if self._calc_type == 'heavy_atom':
                keyword = ' and not hydrogen'
            elif self._calc_type == 'backbone':
                keyword = ' and backbone'
            elif self._calc_type == 'CA':
                keyword = ' and name CA'
            else:
                print( 'Unknown calculation type. Will take all atoms into account.' )
                keyword = ''
            
            dist_matrix = np.full( (num_steps - 1, len( rec_resids ), len( lig_resids )), float( self._cutoff + 5 ) )
            not_empty_xindex = []
            not_empty_yindex = []

            for r in tqdm( range( len( rec_resids ) ) ):  # control by residue index
                r_id = rec_resids[r]
                for l in range( len( lig_resids ) ):
                    l_id = lig_resids[l]

                    receptor = atomsel( "resid " + str( r_id ) + keyword + " and chain " + self._receptor_chain_ID )
                    ligand = atomsel( "resid " + str( l_id ) + keyword + " and chain " + self._ligand_chain_ID )

                    for t in range( 1, num_steps ):
                        receptor.frame = t
                        ligand.frame = t
                        receptor.update()
                        ligand.update()
                        dist = calc_dist( receptor, ligand )
                        #if dist < self._cutoff:
                        dist_matrix[t - 1][r][l] = dist
                        receptor_labels.append(
                                "%d" % (r_id) + ' ' + _resid2resname( r_id, self._receptor_chain_ID )[0] )
                        ligand_labels.append(
                                "%d" % (l_id) + ' ' + _resid2resname( l_id, self._ligand_chain_ID )[0] )
                        not_empty_xindex.append( r )
                        not_empty_yindex.append( l )

            not_empty_xindex = sorted( list( set( not_empty_xindex ) ) )
            not_empty_yindex = sorted( list( set( not_empty_yindex ) ) )
            filtered_dist_matrix = np.full( (num_steps - 1, len( not_empty_xindex ), len( not_empty_yindex )), float( self._cutoff + 5  ))
            
            for i in range( len( not_empty_xindex ) ):
                for j in range( len( not_empty_yindex ) ):
                    for t in range( 1, num_steps ):
                        filtered_dist_matrix[t-1][i][j] = dist_matrix[t-1][not_empty_xindex[i]][not_empty_yindex[j]]


        else:
            print( 'Distance calculate based on pdb frame.' )
            # top_pdb = molecule.load('pdb',_pdb_file)
            
            if self._calc_type == 'heavy_atom':
                keyword = ' and not hydrogen'
            elif self._calc_type == 'backbone':
                keyword = ' and backbone'
            elif self._calc_type == 'CA':
                keyword = ' and name CA'
            else:
                print( 'Unknown calculation type. Will take all atoms into account.' )
                keyword = ''

            dist_matrix = np.full( (len( rec_resids ), len( lig_resids )), float( self._cutoff + 5 ) )
            not_empty_xindex = []
            not_empty_yindex = []

            for r in tqdm( range( len( rec_resids ) ) ):  # control by residue index
                r_id = rec_resids[r]
                for l in range( len( lig_resids ) ):
                    l_id = lig_resids[l]

                    receptor = atomsel( "resid " + str( r_id ) + keyword + " and chain " + self._receptor_chain_ID )
                    ligand = atomsel( "resid " + str( l_id ) + keyword + " and chain " + self._ligand_chain_ID )

                    dist = calc_dist( receptor, ligand )
                    #if dist < self._cutoff:
                    dist_matrix[r][l] = dist
                    receptor_labels.append(
                            "%d" % (r_id) + ' ' + _resid2resname( r_id, self._receptor_chain_ID )[0] )
                    ligand_labels.append( "%d" % (l_id) + ' ' + _resid2resname( l_id, self._ligand_chain_ID )[0] )
                    not_empty_xindex.append( r )
                    not_empty_yindex.append( l )

            #not_empty_xindex = sorted( list( set( not_empty_xindex ) ) )
            #not_empty_yindex = sorted( list( set( not_empty_yindex ) ) )
            #filtered_dist_matrix = np.zeros( (len( not_empty_xindex ), len( not_empty_yindex )) )
            #for i in range( len( not_empty_xindex ) ):
            #    for j in range( len( not_empty_yindex ) ):
            #        filtered_dist_matrix[i][j] = dist_matrix[not_empty_xindex[i]][not_empty_yindex[j]]
        filtered_dist_matrix = dist_matrix
        receptor_labels = sorted( list( set( receptor_labels ) ), key = lambda x: (len (x), x) )
        ligand_labels = sorted( list( set( ligand_labels ) ), key = lambda x: (len (x), x) )
        return filtered_dist_matrix, receptor_labels, ligand_labels

###########################
# Wrap up above functions

def wrap_traj_freq(top_pdb, traj_dcd, receptor_chain_ID, ligand_chain_ID, calc_type, cutoff, frequency_threshold,
                   stride, bootstrap, receptor_resids, lig_resids):
    """
    This part will help to wrap up the trajectory contact frequency function, and save the matrix and labels into files.
    Before calculating, it will automatically detect if there is an existed result.
    
    Input: top_pdb, traj_dcd, and other parameters
    Output: frequency matrix, receptor label,  signaling protein label.
    
    """
    record_keyword = os.path.basename( top_pdb )[:-4] + '_chain' + receptor_chain_ID + 'vs' + ligand_chain_ID
    checkpoint_keyword = record_keyword + '_' + calc_type+'_'+str(cutoff)
    if not os.path.isdir('checkpoint'):
        os.mkdir('checkpoint')
    if not os.path.exists( os.path.join('checkpoint' ,checkpoint_keyword + '_dist.npy') ) and not os.path.exists( os.path.join('checkpoint',
            checkpoint_keyword + '_xlabel.dat') ) and not os.path.exists( os.path.join('checkpoint',checkpoint_keyword + '_ylabel.dat') ):

        step = Frequency( top_pdb, traj_dcd, receptor_chain_ID, ligand_chain_ID, calc_type, cutoff, frequency_threshold,
                          stride, bootstrap, receptor_resids, lig_resids )
        frequency_matrix, receptor_labels, ligand_labels = step.contact_freq_vmd

        np.save( os.path.join('checkpoint', checkpoint_keyword + '_dist.npy'), frequency_matrix )
        f = open( os.path.join('checkpoint',checkpoint_keyword + '_xlabel.dat'), 'w' )
        for r in receptor_labels:
            f.write( r + '\n' )
        f.close()
        f = open( os.path.join('checkpoint', checkpoint_keyword + '_ylabel.dat'), 'w' )
        for r in ligand_labels:
            f.write( r + '\n' )
        f.close()
    else:
        frequency_matrix = np.load( os.path.join('checkpoint', checkpoint_keyword + '_dist.npy'), allow_pickle=True )
        f = open( os.path.join('checkpoint', checkpoint_keyword + '_xlabel.dat'), 'r' )
        lines = f.read().split( '\n' )
        f.close()
        receptor_labels = [line for line in lines if not line == '']
        f = open( os.path.join('checkpoint', checkpoint_keyword + '_ylabel.dat'), 'r' )
        lines = f.read().split( '\n' )
        f.close()
        ligand_labels = [line for line in lines if not line == '']

    return frequency_matrix, receptor_labels, ligand_labels


def wrap_pdb_distance(top_pdb, receptor_chain_ID, ligand_chain_ID, cutoff, calc_type, traj_dcd, stride, receptor_resids,
                      lig_resids):
    """
    This part will help to wrap up the pdb contact distance function, and save the matrix and labels into files.
    Before calculating, it will automatically detect if there is an existed result.

    Input: top_pdb, traj_dcd, and other parameters
    Output: distance matrix, receptor label,  signaling protein label.

    """
    record_keyword = os.path.basename( top_pdb )[:-4] + '_chain' + receptor_chain_ID + 'vs' + ligand_chain_ID
    checkpoint_keyword = record_keyword + '_' + calc_type+'_'+str(cutoff)
    if not os.path.isdir('checkpoint'):
        os.mkdir('checkpoint')
    if not os.path.exists( os.path.join('checkpoint' ,checkpoint_keyword + '_dist.npy') ) and not os.path.exists( os.path.join('checkpoint',
            checkpoint_keyword + '_xlabel.dat') ) and not os.path.exists( os.path.join('checkpoint',checkpoint_keyword + '_ylabel.dat') ):

        step = Distance( top_pdb, receptor_chain_ID, ligand_chain_ID, cutoff, calc_type, traj_dcd, stride, receptor_resids,
                         lig_resids )
        distance_matrix, receptor_labels, ligand_labels = step.contact_dist_vmd()

        np.save( os.path.join('checkpoint', checkpoint_keyword + '_dist.npy'), distance_matrix )
        f = open( os.path.join('checkpoint', checkpoint_keyword + '_xlabel.dat'), 'w' )
        for r in receptor_labels:
            f.write( r + '\n' )
        f.close()
        f = open( os.path.join('checkpoint', checkpoint_keyword + '_ylabel.dat'), 'w' )
        for r in ligand_labels:
            f.write( r + '\n' )
        f.close()
    else:
        distance_matrix = np.load( os.path.join('checkpoint',checkpoint_keyword + '_dist.npy'), allow_pickle=True )
        f = open( os.path.join('checkpoint', checkpoint_keyword + '_xlabel.dat'), 'r' )
        lines = f.read().split( '\n' )
        f.close()
        receptor_labels = [line for line in lines if not line == '']
        f = open( os.path.join('checkpoint', checkpoint_keyword + '_ylabel.dat'), 'r' )
        lines = f.read().split( '\n' )
        f.close()
        ligand_labels = [line for line in lines if not line == '']

    return distance_matrix, receptor_labels, ligand_labels


def wrap_traj_distance(top_pdb, receptor_chain_ID, ligand_chain_ID, cutoff, calc_type, traj_dcd, stride,
                       receptor_resids,
                       lig_resids):
    """
    This part will help to wrap up the pdb contact distance function, and save the matrix and labels into files.
    Before calculating, it will automatically detect if there is an existed result.

    Input: top_pdb, traj_dcd, and other parameters
    Output: distance matrix, receptor label,  signaling protein label.

    """
    record_keyword = os.path.basename( top_pdb )[:-4] + '_chain' + receptor_chain_ID + 'vs' + ligand_chain_ID
    checkpoint_keyword = record_keyword + '_' + calc_type+'_'+str(cutoff)
    if not os.path.isdir('checkpoint'):
        os.mkdir('checkpoint')
    if not os.path.exists( os.path.join('checkpoint' ,checkpoint_keyword + '_dist.npy') ) and not os.path.exists( os.path.join('checkpoint',
            checkpoint_keyword + '_xlabel.dat') ) and not os.path.exists( os.path.join('checkpoint',checkpoint_keyword + '_ylabel.dat') ):
        step = Distance( top_pdb, receptor_chain_ID, ligand_chain_ID, cutoff, calc_type, traj_dcd, stride,
                         receptor_resids,
                         lig_resids )
        distance_matrix, receptor_labels, ligand_labels = step.contact_dist_vmd()

        np.save( os.path.join('checkpoint' ,checkpoint_keyword + '_dist.npy'), distance_matrix )
        f = open( os.path.join('checkpoint', checkpoint_keyword + '_xlabel.dat'), 'w' )
        for r in receptor_labels:
            f.write( r + '\n' )
        f.close()
        f = open( os.path.join('checkpoint', checkpoint_keyword + '_ylabel.dat'), 'w' )
        for r in ligand_labels:
            f.write( r + '\n' )
        f.close()
    else:
        distance_matrix = np.load( os.path.join('checkpoint', checkpoint_keyword + '_dist.npy'), allow_pickle=True )
        f = open( os.path.join('checkpoint', checkpoint_keyword + '_xlabel.dat'), 'r' )
        lines = f.read().split( '\n' )
        f.close()
        receptor_labels = [line for line in lines if not line == '']
        f = open( os.path.join('checkpoint', checkpoint_keyword + '_ylabel.dat'), 'r' )
        lines = f.read().split( '\n' )
        f.close()
        ligand_labels = [line for line in lines if not line == '']

    return distance_matrix, receptor_labels, ligand_labels
