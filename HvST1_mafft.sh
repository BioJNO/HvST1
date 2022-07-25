#!/bin/bash

#SBATCH --mail-user=jamie.orr@hutton.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

mafft --maxiterate 1000 --localpair --thread 8 all_longest_orthos.fasta > mafft_longest_orthos.fasta

