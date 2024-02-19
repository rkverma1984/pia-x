'''
Parameters for dat
'''



## 
TMBB = 'backbone and \
       (resid 65:92 or resid 96:124 or \
       resid 136:172 or resid 238:255 or \
       resid 262:284 or resid 308:335 or \
       resid 343:374 or resid 405:436 or \
       resid 445:465 or resid 470:496 or \
       resid 519:540 or resid 558:584)'



## general atom selection
prot    = 'protein'
protbb  = 'protein and backbone'
protbb  = 'protein and name CA'
lig     = 'not protein and not resname SOD CLA'
lig_noh = 'not type H and not protein and not resname SOD CLA'
ions3   = 'resname SOD CLA'
lig_neighbor = 'protein and not type H and around 7 (not protein and not resname SOD CLA)'

