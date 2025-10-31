#!/bin/bash

### Updated by Alan Bergland
## version=3.0
### Oct 20, 2025


check_exit_status () {
  if [ ! "$2" -eq "0" ]; then
    echo "Step $1 failed with exit status $2"
    exit $2
  fi
  echo "Checked step $1"
}

#### DEFAULT parameters
do_single_end=0
read1="test_file_1"
read2="test_file_2"
sample="default_sample_name"
output="."
threads="1"
max_cov=0.95
min_cov=10
theta=0.005
D=0.01
priortype="informative"
fold="unfolded"
maxsnape=0.9
nflies=40
base_quality_threshold=25
illumina_quality_coding=1.8
minIndel=5
###do_prep=1
do_snape=1
do_poolsnp=1
ref_genome="path_to_ref_fasta"
focalFile="path_to_focalFile_csv"
prepRef=1

### Parse positional arguments
# Credit: https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
	-do_se|--single_end)
	  do_single_end=1
    shift # past argument
	  ;;
    -dps|--do_poolsnp)
    do_poolsnp=$2
    shift # past argument
    ;;
    #-dnp|--do_not_prep)
    #do_prep=0
    #shift # past argument
    #;;
    -bq|--base-quality-threshold)
    base_quality_threshold="$2"
    shift # past argument
    shift # past value
    ;;
    -ill|--illumina-quality-coding)
    illumina_quality_coding="$2"
    shift # past argument
    shift # past value
    ;;
    -mindel|--min-indel)
    minIndel="$2"
    shift # past argument
    shift # past value
    ;;
    -tr|--threads)
    threads="$2"
    shift # past argument
    shift # past value
    ;;
    -x|--max-cov)
    max_cov="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--min-cov)
    min_cov="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--help)
    echo "Usage:"
    echo "  singularity run [options] <image> -h          Display this help message."
    echo "  Pair-end reads mode add --sequencing = 'pe'; Single-end mode add --sequencing = 'se'"
    echo "	Pair-end reads mode requires:"
    echo "  singularity run [options] <image> <fastq_file_1_path> <fastq_file_2_path> <sample_name> <output_dir> <num_cores>"
    echo "	Single-end reads mode requires:"
    echo "  singularity run [options] <image> <fastq_file_Single_path> <sample_name> <output_dir> <num_cores>"
    exit 0
    shift # past argument
    ;;
    -t|--theta)
    theta=$2
    shift # past argument
    shift # past value
    ;;
    -D)
    D=$2
    shift # past argument
    shift # past value
    ;;
    -p|--priortype)
    theta=$2
    shift # past argument
    shift # past value
    ;;
    -f|--fold)
    fold=$2
    shift # past argument
    shift # past value
    ;;
    -ms|--maxsnape)
    maxsnape=$2
    shift # past argument
    shift # past value
    ;;
    -nf|--num-flies)
    nflies=$2
    shift # past argument
    shift # past value
    ;;
    -ds|--do_snape)
    do_snape=$2
    shift # past argument
    shift # past value
    ;;
    -ref|--reference_genome)
    ref=$2
    shift # past argument
    shift # past value
    ;;
    -focalFile|--focal_file)
    focalFile=$2
    shift # past argument
    shift # past value
    ;;
    -prepRef|--prep_reference)
    prepRef=$2
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") #save it to an array
    shift
    ;;
esac
done

#################
### prepare reference genome
### module load samtools bwa picard
echo $prepRef
echo $do_poolsnp
echo $do_snape
if [ $prepRef -eq "1" ] && [ $do_poolsnp -eq "1" ]; then
  echo "Cannot prep ref and run mapping at once"
  exit 1
fi

if [ $prepRef -eq "1" ] && [ $do_snape -eq "1" ]; then
  echo "Cannot prep ref and run mapping at once"
  exit 1
fi


if [ $prepRef -eq "1" ]; then
  if [ ! -f ${ref}.amb ]; then bwa index ${ref}; fi
  if [ ! -f ${ref}.fai ]; then samtools faidx ${ref}; fi
  if [ ! -f ${ref}.dict ]; then
    refDict=$( echo ${ref} | sed 's/fa/dict/g' )
    java -jar $PICARD CreateSequenceDictionary R=${ref} O=${refDict}
  fi

  while read p; do
     #prefix=mel
     prefix=$( echo $p | cut -f1 -d',')
     echo $prefix
     chrs=$( echo $p | cut -f2 -d',')
     echo $chrs
     refOut=$( echo ${ref} | sed "s/fa/${prefix}.fa/g" )
     samtools faidx ${ref} ${chrs} > ${refOut}

     python3 /opt/DESTv3/mappingPipeline/scripts/PickleRef.py \
         --ref ${refOut} \
         --output ${refOut}
  done < ${focalFile}
  exit
fi


#####
if [ $do_single_end -eq "1" ]
	then
	echo "Processing ---> single end run"

fi

if [ $do_single_end -eq "0" ]
	then
	echo "Processing ---> paired end run"

fi

#####
#####

set -- "${POSITIONAL[@]}"

#### --> Evaluate whether the pipeline is being used as Paired end or single end

# if [ $do_single_end -eq "0" ] && [ $# != 4 ]
#   then
#     echo "ERROR: For paired end reads (default) you need to supply <fastq_file_1_path> <fastq_file_2_path> <sample_name> <output_dir> as positional arguments, and no others"
#     exit 1
# fi
#
# if [ $do_single_end -eq "1" ] && [ $# != 3 ]
#   then
#     echo "ERROR: For single end reads you need to supply <fastq_file_path> <sample_name> <output_dir> as positional arguments, and no others"
#     exit 1
# fi

######

if [ $do_single_end -eq "0"  ]
	then
read1=$1; shift
read2=$1; shift
sample=$1; shift
output=$1; shift
fi

if [ $do_single_end -eq "1" ]
 	then
read1=$1; shift
sample=$1; shift
output=$1; shift
fi

########

echo -e "This is DEST v. ${version} \n Parameters as interpreted + those assumed by default --> \n"
echo -e \
"pe (0) or se (1)? =" $do_single_end "\n" \
"r1 =" $read1 "\n" \
"r2 =" $read2 "\n" \
"sample name =" $sample "\n" \
"output =" $output "\n" \
"number of flies =" $nflies "\n" \
"cpus =" $threads "\n" \
"max cov =" $max_cov "\n" \
"min cov =" $min_cov "\n" \
"theta =" $theta "\n" \
"D =" $D "\n" \
"prior =" $priortype "\n" \
"folded? =" $fold "\n" \
"max snape ="$maxsnape "\n" \
"base quality threshold =" $base_quality_threshold "\n" \
"illumina quality coding =" $illumina_quality_coding "\n" \
"minIndel =" $minIndel "\n" \
"do Prep =" $do_prep "\n" \
"do snape? (0 = no; 1 = yes) -->" $do_snape "\n" \
"do poolsnp? (0 = no; 1 = yes) -->" $do_poolsnp "\n" \
"reference genome =" $ref "\n" \
"focal chromosome file=" $focalFile "\n" \
"prep ref?=" $prepRef "\n"

########

if [ ! -f "$read1" ] && [ $do_single_end -eq "0"  ]; then
  echo "ERROR: for paired end run"
  echo "ERROR DETAILS: $read1 does not exist"
  exit 1
fi

if [ ! -f "$read2" ] && [ $do_single_end -eq "0"  ]; then
  echo "ERROR: for paired end run"
  echo "ERROR DETAILS: $read2 does not exist"
  exit 1
fi

if [ ! -f "$read1" ] && [ $do_single_end -eq "1"  ]; then
  echo "ERROR: for single end run"
  echo "ERROR DETAILS: $read1 does not exist"
  exit 1
fi

if [ ! -d $output/$sample/ ]; then
  mkdir -p $output/$sample/
fi

if [ ! -d $output/$sample/${sample}_fastqc ]; then
  mkdir $output/$sample/${sample}_fastqc
fi

if [ ! -d $output/$sample/${sample}_fastqc/trimmed ]; then
  mkdir $output/$sample/${sample}_fastqc/trimmed
fi


##### Begin prep
#### Single end case

  if [ $do_single_end -eq "1" ]
	then

  echo "begin read preparation | Single End Mode"

  fastqc $read1 -o $output/$sample/${sample}_fastqc

  cutadapt \
  -q 18 \
  --minimum-length 25 \
  -o $output/$sample/${sample}.trimmed1.fq.gz \
  -b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
  -O 15 \
  -n 3 \
  --cores=$threads \
  $read1

  check_exit_status "cutadapt" $?

  fastqc $output/$sample/${sample}.trimmed1.fq.gz  -o $output/$sample/${sample}_fastqc/trimmed

  check_exit_status "fastqc" $?

  #Automatically uses all available cores
  #bbmerge.sh \
  #in=$output/$sample/${sample}.trimmed1.fq.gz \
  #out=$output/$sample/${sample}.merged.fq.gz \
  #outu=$output/$sample/${sample}.1_un.fq.gz \

  #check_exit_status "bbmerge" $?

  #rm $output/$sample/${sample}.trimmed*

fi

#### Paired end case

  if [ $do_single_end -eq "0" ]
	then

  echo "begin read preparation | Paired End Mode"

  fastqc $read1 $read2 -o $output/$sample/${sample}_fastqc

  cutadapt \
  -q 18 \
  --minimum-length 25 \
  -o $output/$sample/${sample}.trimmed1.fq.gz \
  -p $output/$sample/${sample}.trimmed2.fq.gz \
  -b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
  -B CAAGCAGAAGACGGCATACGAGAT \
  -O 15 \
  -n 3 \
  --cores=$threads \
  $read1 $read2

  check_exit_status "cutadapt" $?

  fastqc $output/$sample/${sample}.trimmed1.fq.gz $output/$sample/${sample}.trimmed2.fq.gz -o $output/$sample/${sample}_fastqc/trimmed

  check_exit_status "fastqc" $?

  #Automatically uses all available cores
  bbmerge.sh in1=$output/$sample/${sample}.trimmed1.fq.gz in2=$output/$sample/${sample}.trimmed2.fq.gz out=$output/$sample/${sample}.merged.fq.gz outu1=$output/$sample/${sample}.1_un.fq.gz outu2=$output/$sample/${sample}.2_un.fq.gz

  check_exit_status "bbmerge" $?

  rm $output/$sample/${sample}.trimmed*

fi

  ##### BWA mem portion
  ##### Map as Paired end

  if [ $do_single_end -eq "0" ]
  then

  bwa mem -t $threads -M -R "@RG\tID:$sample\tSM:sample_name\tPL:illumina\tLB:lib1" ${ref} $output/$sample/${sample}.1_un.fq.gz $output/$sample/${sample}.2_un.fq.gz | samtools view -@ $threads -Sbh -q 20 -F 0x100 - > $output/$sample/${sample}.merged_un.bam

  rm $output/$sample/${sample}.1_un.fq.gz
  rm $output/$sample/${sample}.2_un.fq.gz

  bwa mem -t $threads -M -R "@RG\tID:$sample\tSM:sample_name\tPL:illumina\tLB:lib1" ${ref} $output/$sample/${sample}.merged.fq.gz | samtools view -@ $threads -Sbh -q 20 -F 0x100 - > $output/$sample/${sample}.merged.bam

  check_exit_status "bwa_mem" $?

  rm $output/$sample/${sample}.merged.fq.gz

  java -jar $PICARD MergeSamFiles I=$output/$sample/${sample}.merged.bam I=$output/$sample/${sample}.merged_un.bam SO=coordinate USE_THREADING=true O=$output/$sample/${sample}.sorted_merged.bam

  echo "Mapped as PE done!"

  fi
  ### ^^^ This closes mapping as PE

  ### Begin Mapping as Single end

  if [ $do_single_end -eq "1" ]
  then

  bwa mem -t $threads -M -R "@RG\tID:$sample\tSM:sample_name\tPL:illumina\tLB:lib1" ${ref} $output/$sample/${sample}.trimmed1.fq.gz | samtools view -@ $threads -Sbh -q 20 -F 0x100 - > $output/$sample/${sample}.merged.bam

  java -jar $PICARD SortSam \
 	I=$output/$sample/${sample}.merged.bam \
 	O=$output/$sample/${sample}.sorted_merged.bam \
 	SO=coordinate \
 	VALIDATION_STRINGENCY=SILENT


  rm $output/$sample/${sample}.trimmed*
  echo "Mapped as SE done!"

  fi
  #### ^^^^ Done mapping as single end

#### Continue with Picard and remove duplication

  check_exit_status "Picard_MergeSamFiles" $?

  rm $output/$sample/${sample}.merged.bam
  rm $output/$sample/${sample}.merged_un.bam

  java -jar $PICARD MarkDuplicates \
  REMOVE_DUPLICATES=true \
  I=$output/$sample/${sample}.sorted_merged.bam \
  O=$output/$sample/${sample}.dedup.bam \
  M=$output/$sample/${sample}.mark_duplicates_report.txt \
  VALIDATION_STRINGENCY=SILENT

  check_exit_status "Picard_MarkDuplicates" $?

  rm $output/$sample/${sample}.sorted_merged.bam

  samtools index $output/$sample/${sample}.dedup.bam

  java -jar $GATK -T RealignerTargetCreator \
  -nt $threads \
  -R ${ref} \
  -I $output/$sample/${sample}.dedup.bam \
  -o $output/$sample/${sample}.hologenome.intervals

  check_exit_status "GATK_RealignerTargetCreator" $?

  java -jar $GATK \
  -T IndelRealigner \
  -R ${ref} \
  -I $output/$sample/${sample}.dedup.bam \
  -targetIntervals $output/$sample/${sample}.hologenome.intervals \
  -o $output/$sample/${sample}.contaminated_realigned.bam

  check_exit_status "GATK_IndelRealigner" $?

  rm $output/$sample/${sample}.dedup.bam*

  # samtools index $output/$sample/${sample}.contaminated_realigned.bam

  #Number of reads mapping to simulans and mel
  # grep -v "sim_" $output/$sample/${sample}.${sample}.original_idxstats.txt | awk -F '\t' '{sum+=$3;} END {print sum;}' > $output/$sample/${sample}.num_mel.txt
  # grep "sim_" $output/$sample/${sample}.${sample}.original_idxstats.txt | awk -F '\t' '{sum+=$3;} END {print sum;}' > $output/$sample/${sample}.num_sim.txt

  #Filter out the simulans contaminants
  while read p; do
     #prefix=mel
     prefix=$( echo $p | cut -f1 -d',')
     echo $prefix
     chrs=$( echo $p | cut -f2 -d',')
     echo $chrs
     refOut=$( echo ${ref} | sed "s/fa/${prefix}.fa/g" )

     samtools view -@ $threads $output/$sample/${sample}.contaminated_realigned.bam ${chrs} -b > $output/$sample/${sample}.${prefix}.bam

     #refOut=/scratch/aob2x/tmpRef/holo_dmel_6.12.sim.fa
     samtools sort --reference ${refOut} -@ ${threads} --write-index -O BAM $output/$sample/${sample}.${prefix}.bam -o $output/$sample/${sample}.${prefix}.sort.bam
     mv $output/$sample/${sample}.${prefix}.sort.bam $output/$sample/${sample}.${prefix}.bam

  done < ${focalFile}

  mv $output/$sample/${sample}.contaminated_realigned.bam  $output/$sample/${sample}.original.bam
  rm $output/$sample/${sample}.contaminated_realigned.bai

  samtools index $output/$sample/${sample}.original.bam
  #samtools mpileup $output/$sample/${sample}.mel.bam -B -f /opt/hologenome/raw/D_melanogaster_r6.12.fasta > $output/$sample/${sample}.mel_mpileup.txt

  check_exit_status "mpileup" $?

if [ $do_poolsnp -eq "1" ]; then
  while read p; do
    prefix=$( echo $p | cut -f1 -d',')
    echo $prefix
    chrs=$( echo $p | cut -f2 -d',')
    echo $chrs
    refOut=$( echo ${ref} | sed "s/fa/${prefix}.fa/g" )

    samtools mpileup $output/$sample/${sample}.${prefix}.bam \
    -B \
    -Q ${base_quality_threshold} \
    -f ${refOut} > $output/$sample/${sample}.${prefix}_mpileup.txt
    check_exit_status "mpileup" $?

    python3 /opt/DESTv3/mappingPipeline/scripts/Mpileup2Sync.py \
    --mpileup $output/$sample/${sample}.${prefix}_mpileup.txt \
    --ref ${refOut}.ref \
    --output $output/$sample/${sample}.${prefix} \
    --base-quality-threshold $base_quality_threshold \
    --coding $illumina_quality_coding \
    --minIndel $minIndel

    check_exit_status "Mpileup2Sync" $?

    #For the PoolSNP output
    python3 /opt/DESTv3/mappingPipeline/scripts/MaskSYNC_snape_complete.py \
    --sync $output/$sample/${sample}.${prefix}.sync.gz \
    --output $output/$sample/${sample}.${prefix} \
    --indel $output/$sample/${sample}.${prefix}.indel \
    --coverage $output/$sample/${sample}.${prefix}.cov \
    --mincov $min_cov \
    --maxcov $max_cov \
    --maxsnape $maxsnape

    check_exit_status "MaskSYNC" $?

    mv $output/$sample/${sample}.${prefix}_masked.sync.gz $output/$sample/${sample}.${prefix}.masked.sync.gz
    gunzip $output/$sample/${sample}.${prefix}.masked.sync.gz
    bgzip $output/$sample/${sample}.${prefix}.masked.sync
    tabix -s 1 -b 2 -e 2 $output/$sample/${sample}.${prefix}.masked.sync.gz

    check_exit_status "tabix" $?

    #rm $output/$sample/${sample}.${prefix}_mpileup.txt

    echo "Read 1: $read1" >> $output/$sample/${sample}.parameters.txt
    echo "Read 2: $read2" >> $output/$sample/${sample}.parameters.txt
    echo "Sample name: $sample" >> $output/$sample/${sample}.parameters.txt
    echo "Output directory: $output" >> $output/$sample/${sample}.parameters.txt
    echo "Number of cores used: $threads" >> $output/$sample/${sample}.parameters.txt
    echo "Max cov: $max_cov" >> $output/$sample/${sample}.parameters.txt
    echo "Min cov $min_cov" >> $output/$sample/${sample}.parameters.txt
    echo "base-quality-threshold $base_quality_threshold" >> $output/$sample/${sample}.parameters.txt
    echo "illumina-quality-coding $illumina_quality_coding" >> $output/$sample/${sample}.parameters.txt
    echo "min-indel $minIndel" >> $output/$sample/${sample}.parameters.txt
    echo "species prefix ${prefix}" >> $output/$sample/${sample}.parameters.txt
  done < ${focalFile}
fi

#Generate the SNAPE SYNC files
if [ $do_snape -eq "1" ]; then
  while read p; do
    prefix=$( echo $p | cut -f1 -d',')
    echo $prefix
    chrs=$( echo $p | cut -f2 -d',')
    echo $chrs
    refOut=$( echo ${ref} | sed "s/fa/${prefix}.fa/g" )

    /opt/DESTv2/mappingPipeline/scripts/Mpileup2Snape.sh \
    ${sample}.${prefix}_mpileup.txt \
    $output \
    $sample.${prefix} \
    $theta \
    $D \
    $priortype \
    $fold \
    $nflies

    check_exit_status "Mpileup2SNAPE" $?

    gzip -f $output/$sample/${sample}.${prefix}.SNAPE.output.txt

    python3 /opt/DESTv3/mappingPipeline/scripts/SNAPE2SYNC.py \
    --input $output/$sample/${sample}.${prefix}.SNAPE.output.txt.gz \
    --ref ${refOut} \
    --output $output/$sample/${sample}.${prefix}.SNAPE

    check_exit_status "SNAPE2SYNC" $?

    python3 /opt/DESTv3/mappingPipeline/scripts/MaskSYNC_snape_complete.py \
    --sync $output/$sample/${sample}.${prefix}.SNAPE.sync.gz \
    --output $output/$sample/${sample}.${prefix}.SNAPE.complete \
    --indel $output/$sample/${sample}.indel \
    --coverage $output/$sample/${sample}.cov \
    --mincov $min_cov \
    --maxcov $max_cov \
    --maxsnape $maxsnape \
    --SNAPE

    check_exit_status "MaskSYNC_SNAPE_Complete" $?

    mv $output/$sample/${sample}.${prefix}.SNAPE.complete_masked.sync.gz $output/$sample/${sample}.${prefix}.SNAPE.complete.masked.sync.gz

    python3 /opt/DESTv3/mappingPipeline/scripts/MaskSYNC_snape_monomorphic_filter.py \
    --sync $output/$sample/${sample}.${prefix}.SNAPE.sync.gz \
    --output $output/$sample/${sample}.${prefix}.SNAPE.monomorphic \
    --indel $output/$sample/${sample}.indel \
    --coverage $output/$sample/${sample}.cov \
    --mincov $min_cov \
    --maxcov $max_cov \
    --maxsnape $maxsnape \
    --SNAPE

    check_exit_status "MaskSYNC_SNAPE_Monomporphic_Filter" $?

    mv $output/$sample/${sample}.${prefix}.SNAPE.monomorphic_masked.sync.gz $output/$sample/${sample}.${prefix}.SNAPE.monomorphic.masked.sync.gz

    gunzip $output/$sample/${sample}.${prefix}.SNAPE.complete.masked.sync.gz
    bgzip $output/$sample/${sample}.${prefix}.SNAPE.complete.masked.sync
    tabix -s 1 -b 2 -e 2 $output/$sample/${sample}.${prefix}.SNAPE.complete.masked.sync.gz

    gunzip $output/$sample/${sample}.${prefix}.SNAPE.monomorphic.masked.sync.gz
    bgzip $output/$sample/${sample}.${prefix}.SNAPE.monomorphic.masked.sync
    tabix -s 1 -b 2 -e 2 $output/$sample/${sample}.${prefix}.SNAPE.monomorphic.masked.sync.gz

    check_exit_status "tabix" $?

    #gzip $output/$sample/${sample}.mel_mpileup.txt

    echo "Maxsnape $maxsnape" >> $output/$sample/${sample}.parameters.txt
    echo "theta:  $theta" >> $output/$sample/${sample}.parameters.txt
    echo "D:  $D" >> $output/$sample/${sample}.parameters.txt
    echo "priortype: $priortype" >> $output/$sample/${sample}.parameters.txt
    echo "species prefix ${prefix}" >> $output/$sample/${sample}.parameters.txt

  done < ${focalFile}
fi
