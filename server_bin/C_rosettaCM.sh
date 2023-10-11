#!/bin/bash

cd $1

$ROSETTA3/source/bin/rosetta_scripts.static.linuxgccrelease \
 -database $ROSETTA3/database/ \
 -in:file:fasta seq.fasta \
 -parser:protocol C_rosettaCM.xml \
 -nstruct 5 \
 -relax:jump_move true \
 -relax:dualspace \
 -out::suffix _singletgt \
 -edensity::mapfile inputmap.map \
 -edensity::mapreso 5.0 \
 -edensity::cryoem_scatterers \
 -beta \
 -default_max_cycles 200
