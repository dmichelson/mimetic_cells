#!/bin/bash

#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-00:05
#SBATCH -p short
#SBATCH --mem-per-cpu=1M
#SBATCH -o logs/dam41_%j.out
#SBATCH -e logs/dam41_%j.err

# cd /n/scratch3/users/d/dam41/homer

# perl configureHomer.pl -install homer
# perl configureHomer.pl -install mm10

cd /n/groups/cbdm-db/dam41

sbatch find_motifs_homer.sh IgG_CT mm10 /n/groups/cbdm-db/dam41/FC_06680
sbatch find_motifs_homer.sh H3K27ac_CT mm10 /n/groups/cbdm-db/dam41/FC_06680
sbatch find_motifs_homer.sh Grhl1 mm10 /n/groups/cbdm-db/dam41/FC_06785
sbatch find_motifs_homer.sh Hnf4a mm10 /n/groups/cbdm-db/dam41/FC_06785
sbatch find_motifs_homer.sh Pou2f3_all mm10 /n/groups/cbdm-db/dam41/FC_06921
