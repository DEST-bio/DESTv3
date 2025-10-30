#!/usr/bin/env bash
#
#SBATCH -J dockerMap # A single job name for the array
#SBATCH -c 20 ### 20 cores
#SBATCH -N 1 # on one node
#SBATCH -t 72:00:00
#SBATCH --mem 90G
#SBATCH -o /scratch/aob2x/logs/RunDest.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/RunDest.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### run as: sbatch /home/aob2x/DESTv3/mappingPipeline/misc/testRun.sh
### sacct -j 4841903
### cat /scratch/aob2x/logs/RunDest.4841903*.err


module load apptainer/1.3.4
#apptainer build -F /scratch/aob2x/dest_v3.sif docker://alanbergland/dest_v3:latest

singularity run \
/scratch/aob2x/dest_v3.sif  \
/project/berglandlab/DEST/raw_reads/DrosEU_3_Jan2023/DrosEu-194_1.fastq.gz \
/project/berglandlab/DEST/raw_reads/DrosEU_3_Jan2023/DrosEu-194_2.fastq.gz \
DE_Bad_Bro_1_2020-07-16 \
/scratch/aob2x/dest_v3_output/ \
--threads 4 \
--max-cov 0.95 \
--min-cov 4 \
--base-quality-threshold 25 \
--num-flies 200 \
--reference_genome /scratch/aob2x/tmpRef/holo_dmel_6.12.fa \
--focal_file /scratch/aob2x/tmpRef/focalFile.csv \
--do_snape 1 \
--do_poolsnp 1 \
-prepRef 0
