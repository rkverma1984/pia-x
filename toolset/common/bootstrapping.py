import mdtraj as md
import random

## skip first 300ns (30ns=125frames)
## bootstrapping 5000 frames
nframeskip=1250
nframepick=5000

## loading dcds
dcd1 = md.load_dcd('mdrun1.dcd', top = "init.pdb")
dcd2 = md.load_dcd('mdrun1.1.dcd', top = "init.pdb")

## merge all dcds
dcds = dcd1[nframeskip:].join(dcd2[nframeskip:], check_topology=True, discard_overlapping_frames=True)

## bootstrapping 1
random.seed(1)
bstp_index = random.choices(range(dcds.n_frames), k=nframepick)
print(bstp_index)
dcd_bstp = dcds[bstp_index]
dcd_bstp.save_dcd("bstp1.dcd")

## bootstrapping 1001
random.seed(1001)
bstp_index = random.choices(range(dcds.n_frames), k=nframepick)
print(bstp_index)
dcd_bstp = dcds[bstp_index]
dcd_bstp.save_dcd("bstp2.dcd")

