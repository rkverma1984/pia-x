"""
provide data regarding trajectories in form of dictionaries.
state1_traj = dictionary containing traj names for state1 (mentioned in sample_parameters.py)
state2_traj = dictionary containing traj names for state2 (mentioned in sample_parameters.py)

for both these dictionaries, the traj name is a key while values are list containing clusterids to be included for analysis.

"""



from collections import OrderedDict

"""
provide which cluster (based on PIA calculation)is needed to be used for each trajectory. See below.
"""


"""
state 1: Provide list of trajectories for state 1.
"""

state1_traj = OrderedDict()
state2_traj = OrderedDict()
state1_traj = {
	'd3fl.dop.1.na.charmm36' : ['1', '2'], #use frames from cluster 1 and 2 for this traj
}

"""
state 2. Provide list of trajectories for state 2.
"""
state2_traj = {
	'd3fl.etq.1.na.charmm36' : ['1', '2'], #use frames from cluster 1 and 2 for this traj
}

trajs = [state1_traj, state2_traj]
