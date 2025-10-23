#!/usr/bin/env bash
#
#SBATCH -J dockerMap # A single job name for the array
#SBATCH -c 11 ### 20 cores
#SBATCH -N 1 # on one node
#SBATCH -t 72:00:00
#SBATCH --mem 90G
#SBATCH -o /scratch/aob2x/logs/RunDest.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/RunDest.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### run as: sbatch /home/aob2x/DESTv3/mappingPipeline/misc/testRun.sh
### sacct -j 4841893
### cat /scratch/aob2x/logs/RunDest.4841893*.out


### modules
  module load apptainer

  singularity run \
  /scratch/aob2x/dest_v3.sif  \
  /standard/BerglandTeach/data/fastq/SRP002024/SRR036932.fastq.gz \
  SRR036932 \
  /scratch/aob2x/dest_v3_output/ \
  --cores 4 \
  --max-cov 0.95 \
  --min-cov 4 \
  --base-quality-threshold 25 \
  --num-flies 200 \
  --reference_genome /scratch/aob2x/tmpRef/holo_dmel_6.12.fa \
  --focal_file /scratch/aob2x/tmpRef/focalFile.csv \
  --do_snape 1 \
  --do_poolsnp 1 \
  --prep_reference 0 \
  -do_se
