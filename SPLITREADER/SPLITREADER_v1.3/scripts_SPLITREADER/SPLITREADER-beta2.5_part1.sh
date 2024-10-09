#!/bin/bash

#############################################
#                                           #  
#               SPLITREADER                 #
#             beta2.5 - part 1              #
#           Baduel et al. 2020              #
#                                           # 
#############################################

##Questions or comments to quadrana(ar)bio.ens.psl.eu

#######################################################################################################################################################################################
# # ARGUMENTS PART

# Verifying of the arguments
if [ "$#" -ne 13 ] && [ "$#" -ne 12 ] && [ "$#" -ne 11 ]; then
  echo ""
  echo "Name : $0"
  echo ""
  echo "Usage: "
  echo "	$0 <bam file name> <bam file directory> <bam file extension> <cohort name> <working directory> <TE library> <script directory> <tmp directory> <entropy> <dust> <LS: library size> <number of threads>"
  echo ""
  echo "Description:"
  echo ""
  echo "	First Part of SPLITREADER, a pipeline that detects Transposable Elements (TE) insertions, from the alignment of samples on their reference genome, based on their short reads paired-end sequences."
  echo "	This first part consists of extracting the unmapped reads of the data."
  echo ""
  echo "Arguments:"
  echo ""
  echo "	<bam file name>:        Name of the BAM file"
  echo "	<bam file directory>:   Directory containing the BAM file"
  echo "	<bam file extension>:   BAM file extension (e.g., '_dupl_fixed_paired.bam')"
  echo "	<cohort name>:          Name of the cohort of genomes analyzed"
  echo "	<working directory>:    Working directory"
  echo "	<TE library>:           TE library file"
  echo "	<script directory>:     Directory containing scripts"
  echo "	<tmp directory>:		Temporary directory"
  echo "	<entropy>:				Entropy score to filter the reads below this value"
  echo "	<dust>:					Dust score to filter the reads above this value"
  echo "    <Memory>:               Memory for the part1, in Giga"
  echo "	<LS: library size>:     Library size (optional : default 300)"
  echo "	<number of threads>:    Number of threads for processing (optional : default 2)"
  echo ""
  exit 1
fi

# Arguments : 

LS_DEFAULT=300
CORES_DEFAULT=2

in=$1 # bam file name
InputDir=$2 # bam file directory
ext=$3 # bam file extension e.g. "_dupl_fixed_paired.bam"
cohort=$4 # name cohort of genomes analyzed
workdir=$5 # working directory
TEfile=$6 # path to TE library
script_dir=$7 # scripts directory
tmp_dir=$8  # tmp directory
entropy=$9 # entropy value for the LC filter
dust=${10} # dust value for the LC filter
RAM=${11} # Memory for part1
LS=${12:-$LS_DEFAULT} # Library Size
pe=YES # Paired-end ?
mult_mapping_limit=4 # number of multiple mapping allowed before reads are rejected (important when TEs are too similar to each other)
#### If not specified, the program estimate them. To this end, the program calculates the minimum number of reads as 3 standard deviation under the mean whole genome coverage
####This value should be at least 3, if not, it is forced to be 3

CORES=${13:-$CORES_DEFAULT} # Number of threads : By default is 2 threads

# PID of this batch
IDPID=$$

# Temporary directory
mkdir -p $tmp_dir/$IDPID                                                  
TmpDir=$tmp_dir/$IDPID          

#############################


# Verifying of the arguments

if [ ! -d "$InputDir" ]; then
  echo ""
  echo "Error: BAMs directory not found!"
  echo ""
  exit 1
fi

if [ ! -f "$InputDir/$in$ext.bam" ]; then
  echo ""
  echo "Error: BAM file not found!"
  echo ""
  exit 1
fi

if [ -z "$cohort" ]; then
  echo ""
  echo "Error: cohortname not specified!"
  echo ""
  exit 1
fi

if [ ! -d "$workdir" ]; then
  echo ""
  echo "Error: Working directory not found!"
  echo ""
  exit 1
fi

if [ ! -d "$TmpDir" ]; then
  echo ""
  echo "Error: Temporary directory not found!"
  echo ""
  exit 1
fi

if [ ! -f "$workdir/TE_sequence/$TEfile.fasta" ]; then
  echo ""
  echo "Error: TE library file not found!"
  echo ""
  exit 1
fi

if ! [[ "$LS" =~ ^[0-9]+$ ]] || [ "$LS" -lt 1 ]; then
  echo ""
  echo "Error: Library size must be a positive integer greater than or equal to 1."
  echo ""
  exit 1
fi

if ! [[ "$CORES" =~ ^[0-9]+$ ]]; then
  echo ""
  echo "Error: Number of threads must be a positive integer."
  echo ""
  exit 1
fi

if [ ! -d "$script_dir" ]; then
  echo ""
  echo "Error: Script directory not found!"
  echo ""
  exit 1
fi


#######################################################################################################################################################################################

OutputDir=$workdir/$cohort/${in}/part1
if [ -e $OutputDir/log.txt ]
then
  rm $OutputDir/log.txt
fi

#Folder containing the Bowtie2 index for combined TE sequence reference
TESequences=$workdir/TE_sequence/$TEfile


readsDir=$TmpDir
function cleanup {
  rm -r -f $TmpDir
  rm -f $TmpDir/${in}*
  }
trap cleanup EXIT



mkdir -p /$workdir/$cohort/$in/part1/


echo "["$(date +"%y-%m-%d %T")"] Running SPLITREADER part 1 for ${in}"  | tee -a $TmpDir/log.txt  


#Extracting unmapped reads
echo "["$(date +"%y-%m-%d %T")"] Extracting unmapped reads from ${in}${ext}" | tee -a $TmpDir/log.txt

if [ -z "$pe" ]
then
  pe=`$samtoolsDir/samtools view -c -f 1 $InputDir/${in}${ext}.bam | awk '{print $1}'` 
else
  if [[ "$pe" == F* ]]
  then
    pe=""
  fi
fi

#######################################################################################################################################################################################
# # Extract unmapped reads and their mates from BAM, then convert them to fastq


# # extract unmapped reads and their mates from BAM
samtools view -h -f 4 $InputDir/${in}${ext}.bam > $TmpDir/${in}${ext}.sam 2>> $TmpDir/log.txt
samtools view -f 8 $InputDir/${in}${ext}.bam >> $TmpDir/${in}${ext}.sam 2>> $TmpDir/log.txt
# convert to BAM
samtools view -Su $TmpDir/${in}${ext}.sam > $TmpDir/${in}${ext}.bam

rm $TmpDir/${in}${ext}.sam

# # Convert BAM to fastq
if [ -z "$pe" ]
  then
 	_JAVA_OPTIONS="-Xmx${RAM}g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp" picard SamToFastq --INPUT $TmpDir/${in}${ext}.bam --FASTQ $readsDir/$in.fastq --TMP_DIR $TmpDir/javatemp  2>> $TmpDir/log.txt
    
  else
	_JAVA_OPTIONS="-Xmx${RAM}g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp" picard SamToFastq --INPUT $TmpDir/${in}${ext}.bam --FASTQ $readsDir/$in.1 --SECOND_END_FASTQ $readsDir/$in.2 --TMP_DIR $TmpDir/javatemp  2>> $TmpDir/log.txt
  

echo "["$(date +"%y-%m-%d %T")"] Extracting discordant reads from ${in}${ext}: Insert size=$LS" | tee -a $TmpDir/log.txt

	samtools view -hF 4 $InputDir/${in}${ext}.bam | awk -v l=$LS '(($4-$8)>(l*10) || ($8-$4)>(l*10) || $7!="=") || $1~/@HD/ || $1~/@SQ/ || $1~/@PG/' | _JAVA_OPTIONS="-Xmx${RAM}g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp" picard SamToFastq -INPUT /dev/stdin -FASTQ $readsDir/$in.1b -SECOND_END_FASTQ $readsDir/$in.2b -TMP_DIR $TmpDir/javatemp  2>> $TmpDir/log.txt


  cat $readsDir/$in.1 $readsDir/$in.1b > $readsDir/$in.1.fastq
  cat $readsDir/$in.2 $readsDir/$in.2b > $readsDir/$in.2.fastq
  rm -f $readsDir/$in.1
  rm -f $readsDir/$in.2
  rm -f $readsDir/$in.1b
  rm -f $readsDir/$in.2b
  
 
fi

rm $TmpDir/${in}${ext}.bam

#######################################################################################################################################################################################
	# # Remove Low Complexity Reads and duplicates

if [ -z "$pe" ]
then
	prinseq++ -fastq $readsDir/$in.1.fastq -out_format 0 -out_name $readsDir/${in} -lc_entropy=$entropy -lc_dust=$dust -threads $CORES
else
	prinseq++ -fastq $readsDir/$in.1.fastq -fastq2 $readsDir/$in.2.fastq -out_format 0 -out_name $readsDir/${in} -lc_entropy=$entropy -lc_dust=$dust -threads $CORES
fi

#######################################################################################################################################################################################
	# # Map fastq to TE joint-reference fasta 

if [ -z "$pe" ]
then
  bowtie2 -x $TESequences -U $readsDir/$in.fastq -S $TmpDir/$in.sam --local --very-sensitive -k $mult_mapping_limit --threads $CORES  2>> $TmpDir/log.txt
else
  bowtie2 -x $TESequences -1 $readsDir/${in}_good_out_R1.fastq -2 $readsDir/${in}_good_out_R2.fastq -S $TmpDir/$in.sam --local --very-sensitive -k $mult_mapping_limit --threads $CORES 2>> $TmpDir/log.txt
fi
      
#######################################################################################################################################################################################
	# # TE Detection: Processing Unmapped Reads

samtools view -H $TmpDir/$in.sam > $TmpDir/$in-TE.sam 2>> $TmpDir/log.txt
  # # Extract reads with soft-clipping, potentially indicating transposable element (TE) insertions
samtools view -F 4 -S $TmpDir/$in.sam | awk '$6~/^[2-8][0-9]S/ || $6~/[2-8][0-9]S$/ || $6~/^1[0-9][0-9]S$/ || $6~/1[0-9][0-9]S$/ {print $0}'>> $TmpDir/$in-TE.sam 2>> $TmpDir/log.txt # -f 1 
  
  # # extract unmapped pair of mapped reads
samtools view -f 4 -F 8 -S $TmpDir/$in.sam | awk '{print $0}' >> $TmpDir/$in-TE.sam 2>> $TmpDir/log.txt # -f 5
  # # extract mapped pair ends of unmapped reads
samtools view -f 8 -F 4 -S $TmpDir/$in.sam | awk '{print $0}' >> $TmpDir/$in-TE.sam 2>> $TmpDir/log.txt # -f 5

  # # convert sam to bam
samtools view -Sb $TmpDir/$in-TE.sam > $OutputDir/$in-TE.bam


echo "["$(date +"%y-%m-%d %T")"] Finished running SPLITREADER part 1 for ${in}"  | tee -a $TmpDir/log.txt

cp $TmpDir/log.txt $OutputDir/

rm -f $readsDir/$in*.fastq
