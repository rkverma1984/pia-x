Pairwise Interaction Analyzer (PIA) was originally applied in our analysis of [molecular dynamics simulations (MD) of LeuT](http://doi.org/10.1074/jbc.M114.625343). The method was fully described in the context of [detecting allosteric pathways in transmembrane molecular machines](http://doi.org/10.1016/j.bbamem.2016.01.010)
(contributor: Sebastian Stolzenberg).

The PIA-X in this repository is a light version of original PIA, focusing on the differential impacts of various ligand bindings on protein conformations, especially those surrounding ligand binding pocket(s). This program was first applied in our analysis of [MD simulations of dopamine D3 receptor](http://doi.org/10.1021/acs.jmedchem.6b01148), and was then adapted for those of [the dopamine transporter (DAT)](http://doi.org/10.1021/acschemneuro.7b00094) and [the serotonin transporter (SERT)](http://doi.org/10.1016/j.neuropharm.2018.10.040)
(contributors in chronological order: Clare Zhu, Mayako Michino, Ravi Kumar Verma, Ara Abramyan, Kuo Hao Lee). The current package was modularized by Ravi Kumar Verma and is maintained by Kuo Hao Lee.

***Inputs***
    parameters_file          : a file containing all the parameters needed for the Pairwise interaction calculations.
                                   Visit sample_parameters.py to have a look on the format.
    traj_data_file           : a file consisting information about trajectory names and cluster_ids needed to be included in analysis.
                                   Visit sample_traj_data.py to have a look on the format.

***Installation***
	Modify python path in ~/.bashrc file:

	export PYTHONPATH=$PYTHONPATH:/home/user/PycharmProjects/Git-repo/protocols/PIA-GPCR/modular_pia_gpcr/:/home/user/PycharmProjects
	/Git-repo/protocols/PIA-GPCR/modular_pia_gpcr/functions/
	alias pia_gpcr='python /home/user/PycharmProjects/Git-repo/protocols/PIA-GPCR/modular_pia_gpcr/pairwise_interactions_calculation_tool.py'

	Change user.

***Usage format***
	pia_gpcr -p sample_parameters.py -t sample_traj_data.py [Arguments]...
	unit test could be run only from the ~/PycharmProjects/Git-repo/protocols/PIA-GPCR/modular_pia_gpcr/ directory.

	To run unit-test do:
		pia_gpcr -d yes [Arguments]...

***Arguments***
	-p, --pfile                                Parameters file name
	-t, --tfile                                Trajectory file name
	-d, --do=yes|no                            Run unit test
	-w, --write-inter=yes|no, default yes      If "yes" write the metrics into .npy files.
	                                           The .npy files stores traj and frame
	                                           ids for each of the individual matrix.
	-o, --out-plt=yes|no|plots-only            Produce differential graphs
	              , default yes                if "only", graphs are generated using previously
	                                           calculated differential metrics
	-c, --color, default RdBu                  specify which colors to use in graphs

Available colors:

	Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu,
	GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1,
	Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples,
	Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1, Set1_r,
	Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r,
	YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg, brg_r, bwr, bwr_r, cool,
	cool_r, coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r,
	gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern,
	gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r,
	inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean, ocean_r, pink, pink_r, plasma,
	plasma_r, prism, prism_r, rainbow, rainbow_r, seismic, seismic_r, spectral, spectral_r, spring, spring_r, summer, summer_r,
	terrain, terrain_r, viridis, viridis_r, winter, winter_r


Description:
	the pairwise_interactions_calculation_tool.py is a python based tool to perform analyses identical to PIA-GPCR.
	Work flawlessly on trajectories.
	
	Calculate four pairwise metrics.
		Distances between the helical subsegments.
		Angle between the helical subsegments.
		Distances between the CA atoms of the binding site residues.
		Chi1 angles for the binding site residues.
	
Variables needed to be modified:
	Modify relevant parameters in parameters file such as name of the trajectories etc.
	provide a valid cluster id input file (see sample_cluster.dat). The name of the file should be defined in
	parameters file.
	
Unit test:
	compares metrics generated using mdtraj with those generated from pymol.
	unit_test.py could be run using:
		pia_gpcr -d yes -w yes


Speed:
	takes about 4 hours to complete analysis for all four metrics on a dataset of 20,000 frames.


** New **
    1. Added functionality to create a logfile
    2. Now work
    3. Added functionality to perform calculations on pdb files as well.
        use trajectory_suffix='pdb' in parameters file to use this.
        use "PDB_ch_0" format in traj_data_file.py. number denotes the Chain ID index.
    4: include in parameters file
        mass=1 if you want to calculate center off mass, or
        mass=0 if centroid is to be calculated.

        if "atomselection_com_TM" in parameters file is CA both mass=1 or mass=0 gives same results.

	5.  helix_angle_calc.py functions has been modified to speed up helix verctor calculations.