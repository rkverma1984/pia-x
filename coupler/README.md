
# Coupler
We developed a program that can quantitively characterize the interactions between rhodopin-like GPCRs and signaling proteins (G-protein and beta-arrestion). Specifically, this method can  
- report the residue pairs between a receptor and a signaling protein that are 
  - below a distance threshold (a single frame, such as crystal or cryo-EM structures)
    - heavy atom-heavy atom (default 5A)
    - Ca-Ca (default 8A)
  - above contact frequency threshold (MD frame ensembles)
- map these residue pairs in matrices 
- compare the matrices to
  - find the unions of receptor (or signaling protein) residues from two (or more) matices
      - based on the union sets of residues, calculate the differential heatmap between two matrices
  - find both the common and unique receptor (or signaling protein) residues from two (or more) matices
      - show the common and/or unique interactions on a matrix
- map the above on 3D structure using PyMol

- [ ] connect to PIA so to show the allosteric pathways from the ligand binding site to receptor/signaling protein interface 

## Getting Started
 **Make soft links of the following** ```vmd_contact.py```, ```_heatmap.py```, ```residue_hash_table.py```, ```wrap_toolkit.py```, ```run.py``` scripts to your working directory.
 
 **Adapt** ```parameters.ini``` and ```$PROTEIN_table.py``` to your own system. These two files should be in the working directory.

### Prerequisites 
The calculations are based VMD-Python API, which requires 
```
python              3.8.3
matplotlib          3.3.1
mdtraj              1.9.4
msgpack             1.0.0
nose                1.3.7
pandas              1.1.2
scipy               1.5.2
seaborn             0.10.1
sympy               1.6.2
tqdm                4.47.0
vmd-python          3.0.6
```
On Cragger, the VMD-Python is installed under both python3.7.6 and python3.8.3.

We can load this moduler by:
```
module load python/anaconda3-2020.07-py3.8.3 
```

To prepare the parameters for a specific receptor, see parameters.ini file in examples folder. To be noticed that, the parameters


### Running examples

- Calculate distance based on PDB frame
```
See pdb_distance_example/CB1_pdb_example.ipynb and D2_pdb_example.ipynb
```
- Calculate contact frequency based on MD trajectory (when necessary, ssh tunnel into a compute node)
```
See traj_freq_example/Traj_sample_contact_frequency.ipynb
```
- Make union list and calculate the difference between two matrices
```
See union_traj_freq_example/Union_trajectory_sample_contact_frequency.ipynb
```
- Visualize contacts in the structures
```
In pymol_script folder, run view_contact.pml. The input parameter file 'cont_freq_resids_list/diff*.py' can be obtained through the above union step.
```
## Contributors

* Bing Xie  
* Lei Shi*  
