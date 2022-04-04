#!/bin/bash

#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-00:05
#SBATCH -p short
#SBATCH --mem-per-cpu=1M
#SBATCH -o logs/dam41_%j.out
#SBATCH -e logs/dam41_%j.err

sbatch merge_bams.sh MEC_IgG_CT mm10 /n/groups/cbdm-db/dam41/FC_06680
sbatch merge_bams.sh MEC_H3K27ac_CT mm10 /n/groups/cbdm-db/dam41/FC_06680
sbatch merge_bams.sh MEC_Grhl1 mm10 /n/groups/cbdm-db/dam41/FC_06785
sbatch merge_bams.sh MEC_Hnf4a mm10 /n/groups/cbdm-db/dam41/FC_06785
sbatch merge_bams.sh MEC_Pou2f3_all mm10 /n/groups/cbdm-db/dam41/FC_06921
