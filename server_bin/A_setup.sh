#!/bin/bash

cd $1

$ROSETTA3/source/bin/partial_thread.static.linuxgccrelease \
 -database $ROSETTA3/database/ \
 -in::file::fasta seq.fasta \
 -in::file::alignment alignment.txt \
 -in::file::template_pdb 1tmpA.pdb
