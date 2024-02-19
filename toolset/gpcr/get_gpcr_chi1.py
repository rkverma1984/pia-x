# define enviroment
import sys, os
from pathlib import Path
home = str(Path.home())
toolset_dir = home + '/repositories/pia-x/toolset'
toolset_common = toolset_dir + "/common"
toolset_gpcr = toolset_dir + "/gpcr"
sys.path.insert(0, toolset_dir)
sys.path.insert(0, toolset_common)
sys.path.insert(0, toolset_gpcr)
from common_gpcr import *
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
from MDAnalysis.analysis.dihedrals import Dihedral
from sys import argv
import pickle
import warnings
warnings.filterwarnings('ignore')


def get_gpcr_rec_chi1(u, resid, verbose=False):
    protein = u.select_atoms('protein and segid PROR and resid ' + str(resid))
    if verbose: print(str(resid))
    try:
        sel_resid_resname = list(set(list(protein.resnames)))[0]
        if verbose: print(str(resid), sel_resid_resname)

        if list(set(list(protein.resnames)))[0] == 'SER':
            try:
                # SET ch1 is defined using N,CA,CB and OG, from http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
                protein = u.select_atoms('protein and segid PROR and resid '+ str(resid) +' and (name N or name CA or name CB or name OG)')
                # protein
                ags = [protein]
                # ags
                # ags[0].names
                R = Dihedral(ags).run()
                mda_chi1 = list(R.angles.T[0])
            except:
                # mda_chi1 = "missing_atom"
                print("missing_atom; not calculated; set chi1 as 9999")
                mda_chi1 = 9999.0
        elif list(set(list(protein.resnames)))[0] == 'VAL':
            # VAL ch1 is defined using N-CA-CB-CG1, from http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
            protein = u.select_atoms('protein and segid PROR and resid '+ str(resid) +' and (name N or name CA or name CB or name CG1)')
            ags = [protein]
            try:
                R = Dihedral(ags).run()
                mda_chi1 = list(R.angles.T[0])
            except:
                print("missing_atom; not calculated; set chi1 as 9999")
                mda_chi1 = 9999.0
        elif list(set(list(protein.resnames)))[0] == 'CYS':
            # CYS ch1 is defined using N-CA-CB-SG, from http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
            protein = u.select_atoms('protein and segid PROR and resid ' + str(resid) + ' and (name N or name CA or name CB or name SG)')
            # protein
            ags = [protein]
            # ags
            # ags[0].names
            R = Dihedral(ags).run()
            mda_chi1 = list(R.angles.T[0])
        elif list(set(list(protein.resnames)))[0] == 'THR':
            # THR ch1 is defined using N-CA-CB-OG1, from http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
            protein = u.select_atoms('protein and segid PROR and resid ' + str(resid) + ' and (name N or name CA or name CB or name OG1)')
            # protein
            ags = [protein]
            # ags
            # ags[0].names
            R = Dihedral(ags).run()
            mda_chi1 = list(R.angles.T[0])
        elif list(set(list(protein.resnames)))[0] == 'ILE':
            # ILE ch1 is defined using N-CA-CB-CG1, from http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
            try:
                protein = u.select_atoms('protein and segid PROR and resid ' + str(resid) + ' and (name N or name CA or name CB or name CG1)')
                # protein
                ags = [protein]
                # ags
                # ags[0].names
                R = Dihedral(ags).run()
                mda_chi1 = list(R.angles.T[0])
            except:
                #mda_chi1 = "missing_atom"
                print("missing_atom; not calculated; set chi1 as 9999")
                mda_chi1 = 9999.0
        elif list(set(list(protein.resnames)))[0] == 'GLY':
            #mda_chi1 = "GLY"
            print("GLY; not calculated; set chi1 as 9999")
            mda_chi1 = 9999.0
        elif list(set(list(protein.resnames)))[0] == 'ALA':
            #mda_chi1 = "ALA"
            print("ALA; not calculated; set chi1 as 9999")
            mda_chi1 = 9999.0
        else:
            try:
                omegas = [res.chi1_selection() for res in protein.residues[0:len(protein.residues)]]
                omegas[0].dihedral.value()
                dihs = dihedrals.Dihedral(omegas).run();
                mda_chi1 = dihs.angles.T[0]
            except:
                if list(set(list(protein.resnames)))[0] == 'HIE':
                    # HIE ch1 is defined using N-CA-CB-CG, from http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
                    protein = u.select_atoms(
                        'protein and segid PROR and resid ' + str(resid) + ' and (name N or name CA or name CB or name CG)')
                    # protein
                    ags = [protein]
                    # ags
                    # ags[0].names
                    R = Dihedral(ags).run()
                    mda_chi1 = list(R.angles.T[0])
                else:
                    print("not calculated; set chi1 as 9999")
                    #mda_chi1 = "na"
                    mda_chi1 = 9999.0

    except:
        print(str(resid), "missing residue")
        mda_chi1 = 9999.0

    if verbose:
        print(mda_chi1)
    return mda_chi1


if __name__ == '__main__':

    # argv
    pdbin = argv[1]
    dcdin = argv[2]
    bwidxin = argv[3]
    gpcrin = argv[4]

    # check
    output_dir = "output"
    check_output_dir(output_dir)

    bwtable = convert_Residue2TMs(gpcrin)


    u = mda.Universe(pdbin, dcdin)
    resid = get_key(bwtable, "TM"+bwidxin)
    mda_chi1 = get_gpcr_rec_chi1(u, resid)

    # save to files
    with open(output_dir + '/gpcr_chi1_' + str(bwidxin) + '.p', 'wb') as fp:
        pickle.dump(mda_chi1, fp)