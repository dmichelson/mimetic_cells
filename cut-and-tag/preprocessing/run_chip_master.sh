#!/bin/bash

#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-00:05
#SBATCH -p short
#SBATCH --mem-per-cpu=1M
#SBATCH -o logs/dam41_%j.out
#SBATCH -e logs/dam41_%j.err

sbatch run_chip.sh LIB050654_CHS00207465_S1 MEC_IgG_CT_R1 mm10 /n/groups/cbdm-db/dam41/FC_06680
sbatch run_chip.sh LIB050654_CHS00207472_S8 MEC_IgG_CT_R2 mm10 /n/groups/cbdm-db/dam41/FC_06680
sbatch run_chip.sh LIB050654_CHS00207471_S7 MEC_H3K27ac_CT_R1 mm10 /n/groups/cbdm-db/dam41/FC_06680
sbatch run_chip.sh LIB050654_CHS00207478_S14 MEC_H3K27ac_CT_R2 mm10 /n/groups/cbdm-db/dam41/FC_06680
sbatch run_chip.sh LIB051310_CHS00211592_S2 MEC_Grhl1_R2 mm10 /n/groups/cbdm-db/dam41/FC_06785
sbatch run_chip.sh LIB051310_CHS00211593_S3 MEC_Grhl1_R3 mm10 /n/groups/cbdm-db/dam41/FC_06785
sbatch run_chip.sh LIB051310_CHS00211594_S4 MEC_Grhl1_R4 mm10 /n/groups/cbdm-db/dam41/FC_06785
sbatch run_chip.sh LIB051310_CHS00211595_S5 MEC_Hnf4a_R1 mm10 /n/groups/cbdm-db/dam41/FC_06785
sbatch run_chip.sh LIB051310_CHS00211596_S6 MEC_Hnf4a_R2 mm10 /n/groups/cbdm-db/dam41/FC_06785
sbatch run_chip.sh LIB051310_CHS00211597_S7 MEC_Hnf4a_R3 mm10 /n/groups/cbdm-db/dam41/FC_06785
sbatch run_chip.sh LIB051310_CHS00211598_S8 MEC_Hnf4a_R4 mm10 /n/groups/cbdm-db/dam41/FC_06785
sbatch run_chip.sh LIB052401_CHS00217421_S1 MEC_Pou2f3_fix_R1 mm10 /n/groups/cbdm-db/dam41/FC_06921
sbatch run_chip.sh LIB052401_CHS00217422_S2 MEC_Pou2f3_fix_R2 mm10 /n/groups/cbdm-db/dam41/FC_06921
sbatch run_chip.sh LIB052401_CHS00217423_S3 MEC_Pou2f3_fix_R3 mm10 /n/groups/cbdm-db/dam41/FC_06921
sbatch run_chip.sh LIB052401_CHS00217424_S4 MEC_Pou2f3_fix_R4 mm10 /n/groups/cbdm-db/dam41/FC_06921
sbatch run_chip.sh LIB052401_CHS00217425_S5 MEC_Pou2f3_nofix_R1 mm10 /n/groups/cbdm-db/dam41/FC_06921
sbatch run_chip.sh LIB052401_CHS00217426_S6 MEC_Pou2f3_nofix_R2 mm10 /n/groups/cbdm-db/dam41/FC_06921
sbatch run_chip.sh LIB052401_CHS00217427_S7 MEC_Pou2f3_nofix_R3 mm10 /n/groups/cbdm-db/dam41/FC_06921
sbatch run_chip.sh LIB052401_CHS00217428_S8 MEC_Pou2f3_nofix_R4 mm10 /n/groups/cbdm-db/dam41/FC_06921