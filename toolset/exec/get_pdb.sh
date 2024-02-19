#!/bin/bash

pdbid=$1

wget http://www.rcsb.org/pdb/files/${pdbid}.pdb.gz
gunzip ${pdbid}.pdb.gz


