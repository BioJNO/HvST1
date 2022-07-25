#!/bin/bash

#SBATCH --mail-user=jamie.orr@hutton.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4

iqtree -B 1000 -nt 4 -s ST1_ortho_alignment.fa
