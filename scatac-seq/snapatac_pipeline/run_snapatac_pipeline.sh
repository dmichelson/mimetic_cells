#!/bin/bash

#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 0-24:00
#SBATCH -p medium
#SBATCH --mem-per-cpu=24G
#SBATCH -o logs/dam41_%j.out
#SBATCH -e logs/dam41_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=daniel_michelson@hms.harvard.edu

module load gcc/6.2.0
module load python/2.7.12
#module load java/jdk-1.8u112
#module load git/2.9.5
#module load bwa/0.7.8
#module load cuda/9.0
#module load sratoolkit/2.8.1
#module load fastqc/0.11.5
#module load bowtie2/2.2.9
#module load samtools/1.3.1
#module load bedtools/2.26.0
#module load homer/4.9
#module load multiqc/1.5
#module load trimmomatic/0.36
#module load R/3.4.1
#module load deeptools/3.0.2
#module load star/2.5.4a
#module load rsem/1.3.0
#module load picard/2.8.0
module load macs2/2.1.1.20160309
export PATH=/n/scratch2/dam41/homer/bin:$PATH

arg1="/n/groups/cbdm-db/dam41/mec_scatac"
arg2="WT1_WT2_KO1_KO2_merged"
arg3="2020-04-13"
arg4="8"

source activate snapatac_env_plusChromVar

cd /n/groups/cbdm-db/dam41/mec_scatac

Rscript scripts/snapatac_pipeline/preprocess_snap.R $arg1 $arg2 $arg3 $arg4
Rscript scripts/snapatac_pipeline/filter_bmat.R $arg1 $arg2 $arg3 $arg4
Rscript scripts/snapatac_pipeline/cluster_viz_all.R $arg1 $arg2 $arg3 $arg4
Rscript scripts/snapatac_pipeline/remove_clusters.R $arg1 $arg2 $arg3 $arg4
Rscript scripts/snapatac_pipeline/cluster_viz_subset.R $arg1 $arg2 $arg3 $arg4
Rscript scripts/snapatac_pipeline/integrate_scrna.R $arg1 $arg2 $arg3 $arg4
Rscript scripts/snapatac_pipeline/integrate_all_genes.R $arg1 $arg2 $arg3 $arg4
Rscript scripts/snapatac_pipeline/refilter_bmat.R $arg1 $arg2 $arg3 $arg4
Rscript scripts/snapatac_pipeline/cluster_clusters.R $arg1 $arg2 $arg3 $arg4
Rscript scripts/snapatac_pipeline/call_ensemble_peaks.R $arg1 $arg2 $arg3 $arg4

./scripts/run_snaptools_addpmat.sh KO2_reprep_Aire_MEChi1_10xATAC_5K $arg2 $arg3 &
./scripts/run_snaptools_addpmat.sh KO_Aire_MEChi1_10xATAC_5K_REA09522 $arg2 $arg3 &
./scripts/run_snaptools_addpmat.sh WT2_reprep_Aire_MEChi1_10xATAC_5K $arg2 $arg3 &
./scripts/run_snaptools_addpmat.sh WT_Aire_MEChi1_10xATAC_5K_REA09522 $arg2 $arg3

wait

Rscript scripts/snapatac_pipeline/find_dars.R $arg1 $arg2 $arg3 $arg4
Rscript scripts/snapatac_pipeline/run_chromvar.R $arg1 $arg2 $arg3 $arg4
Rscript scripts/snapatac_pipeline/run_homer_motif.R $arg1 $arg2 $arg3 $arg4