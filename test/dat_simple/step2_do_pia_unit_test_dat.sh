#!/bin/bash

module load python/anaconda-stock-py2.7.12
module switch python/anaconda-stock-py2.7.12


rm -rf pia_average pia_frames pia_plots
python ~/PycharmProjects/Git-repo/protocols/PIA-GPCR/modular_pia_gpcr/pairwise_interactions_calculation_tool.py -d yes



