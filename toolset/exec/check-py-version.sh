#!/bin/bash

## Andy's version
#module load python/anaconda3-5.1-py3.6.4
## khlee version
#module load python/anaconda3-2020.02-py3.7.6


which python

echo "

"

## using pip to check version
pip freeze | grep 'scikit\|MolVS'


## using python to check version
python << EOF

import rdkit
print("rdkit version: ",rdkit.__version__)

from rdkit import rdBase
print("rdBase.rdkitVersion : ",rdBase.rdkitVersion)
print("rdBase.boostVersion : ",rdBase.boostVersion)


import xgboost
print("xgboost version: ",xgboost.__version__)

EOF

echo " 

"
pip freeze
