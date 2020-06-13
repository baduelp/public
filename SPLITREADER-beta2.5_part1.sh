#!/bin/bash 

#############################################
#                                           #  
#               SPLITREADER                 #
#             beta2.5 - part 1              #
#           Baduel et al. 2020              #
#                                           # 
#############################################

##Questions or comments to quadrana(ar)bio.ens.psl.eu

in=$1 # bam file name
TE_fam=$2 # HERE PART1
InputDir=$3 # bam file directory
ext=$4 # bam file extension e.g. "_dupl_fixed_paired.bam"
cohort=$5 # name cohort of genomes analyzed
workdir=$6 # working directory
TEfile=$7 # directory with TE annotations files

#############################################################
#IN THE FOLLOWING VARIABLE YOU CAN PROVIDE THE MINIMUM READ LEGTH. By default this is 100nt
LENGTH=100

#### If not specified, the program will calculate the longest read


#############################################################
#IN THE FOLLOWING VARIABLE YOU CAN PROVIDE THE MINIMUM NUMBER OF SPLIT-READS IN EACH EXTRIMITY OF THE TE. By default is 5 reads
READS=2
maxcov=100
LS=200
pe=YES
mult_mapping_limit=4 # number of multiple mapping allowed before reads are rejected (important when TEs are too similar to each other)
#### If not specified, the program estimate them. To this end, the program calculates the minimum number of reads as 3 standard deviation under the mean whole genome coverage
####This value should be at least 3, if not, it is forced to be 3
 
#############################################################

#############################################################
#IN THE FOLLOWING VARIABLE YOU CAN EXPLICITE THE NUMBER OF THREADS YOU WANT TO USE FOR MAPPING. By default is 2 threads
CORES=2
#############################################################


########################### edit the following paths !!! ####################


OutputDir=/$workdir/$cohort/${in}/$TE_fam
if [ -e $OutputDir/log.txt ]
then
  rm $OutputDir/log.txt
fi

#Folder containing the Bowtie2 index for combined TE sequence reference
TESequences=/$workdir/TE_sequence/$TEfile

#Bowtie2 executable path
Bowtie2Dir=/usr/local/bin
# #samtools executable path
samtoolsDir=/usr/local/bin
#bedtools executable path
bedtoolsdir=/usr/local/bin
#picard tools executable path
picardDir=/usr/share/java/picard.jar
# PID of this batch
IDPID=$$

# Temporary directory
# TmpDir=/$workdir/localtmp/QD-$in
TmpDir=/localtmp/PB-$IDPID
mkdir -p $TmpDir
readsDir=$TmpDir
function cleanup {
  rm -r -f $TmpDir
  rm -f $TmpDir/${in}*
  }
 trap cleanup EXIT


echo "["$(date +"%y-%m-%d %T")"] Running SPLITREADER 1.5 part 1"  | tee -a $TmpDir/log.txt  


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

# # extract unmapped reads and their mates from BAM
samtools view -h -f 4 $InputDir/${in}${ext}.bam > $TmpDir/${in}${ext}.sam 2>> $TmpDir/log.txt
samtools view -f 8 $InputDir/${in}${ext}.bam >> $TmpDir/${in}${ext}.sam 2>> $TmpDir/log.txt
# convert to BAM
samtools view -Su $TmpDir/${in}${ext}.sam > $TmpDir/${in}${ext}.bam

# # Convert BAM to FASTQ
if [ -z "$pe" ]
  then
  # java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp -jar $picardDir/SamToFastq.jar  INPUT=$TmpDir/${in}${ext}.bam FASTQ=$readsDir/$in.fastq TMP_DIR=$TmpDir/javatemp  2>> $TmpDir/log.txt
  java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp -cp $picardDir net.sf.picard.sam.SamToFastq  INPUT=$TmpDir/${in}${ext}.bam FASTQ=$readsDir/$in.fastq TMP_DIR=$TmpDir/javatemp  2>> $TmpDir/log.txt

    
  else
  java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp -cp $picardDir net.sf.picard.sam.SamToFastq INPUT=$TmpDir/${in}${ext}.bam FASTQ=$readsDir/$in.1 SECOND_END_FASTQ=$readsDir/$in.2 TMP_DIR=$TmpDir/javatemp  2>> $TmpDir/log.txt
   
  if [ -z "$LS" ]
    then
    LS=300
    echo "["$(date +"%y-%m-%d %T")"] Extracting discordant reads from ${in}${ext}: Insert size=$LS [Default]"  | tee -a $TmpDir/log.txt
    else
    echo "["$(date +"%y-%m-%d %T")"] Extracting discordant reads from ${in}${ext}: Insert size=$LS [User Defined]" | tee -a $TmpDir/log.txt
  fi
  samtools view -hF 4 $InputDir/${in}${ext}.bam | awk -v l=$LS '(($4-$8)>(l*10) || ($8-$4)>(l*10) || $7!="=") || $1~/@HD/ || $1~/@SQ/ || $1~/@PG/' | java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp -cp $picardDir net.sf.picard.sam.SamToFastq INPUT=/dev/stdin FASTQ=$readsDir/$in.1b SECOND_END_FASTQ=$readsDir/$in.2b TMP_DIR=$TmpDir/javatemp  2>> $TmpDir/log.txt

  cat $readsDir/$in.1 $readsDir/$in.1b > $readsDir/$in.1.fastq
  cat $readsDir/$in.2 $readsDir/$in.2b > $readsDir/$in.2.fastq
  rm -f $readsDir/$in.1
  rm -f $readsDir/$in.2
  rm -f $readsDir/$in.1b
  rm -f $readsDir/$in.2b
  
 
fi
 
# # Map fastq to TE joint-reference fasta 
if [ -z "$pe" ]
then
  $Bowtie2Dir/bowtie2 -x $TESequences -U $readsDir/$in.fastq -S $TmpDir/$in.sam --local --very-sensitive -k $mult_mapping_limit --threads $CORES  2>> $TmpDir/log.txt
else
  $Bowtie2Dir/bowtie2 -x $TESequences -1 $readsDir/$in.1.fastq -2 $readsDir/$in.2.fastq -S $TmpDir/$in.sam --local --very-sensitive -k $mult_mapping_limit --threads $CORES   2>> $TmpDir/log.txt
fi
      
#############################################################
  ###extract soft-clipped reads with at least 20nt softclipped at 5' or 3' read's end (S>20 or S<$length-20??)
samtools view -H $TmpDir/$in.sam > $TmpDir/$in-TE.sam 2>> $TmpDir/log.txt
samtools view -F 4 -S $TmpDir/$in.sam | awk '$6~/^[2-8][0-9]S/ || $6~/[2-8][0-9]S$/ || $6~/^1[0-9][0-9]S$/ || $6~/1[0-9][0-9]S$/ {print $0}'>> $TmpDir/$in-TE.sam 2>> $TmpDir/log.txt # -f 1 
  
  ### extract unmapped pair of mapped reads
samtools view -f 4 -F 8 -S $TmpDir/$in.sam | awk '{print $0}' >> $TmpDir/$in-TE.sam 2>> $TmpDir/log.txt # -f 5
  ### extract mapped pair ends of unmapped reads
samtools view -f 8 -F 4 -S $TmpDir/$in.sam | awk '{print $0}' >> $TmpDir/$in-TE.sam 2>> $TmpDir/log.txt # -f 5

### convert sam to bam
# cp $TmpDir/$in-TE.sam $OutputDir/$in-TE.sam
samtools view -Sb $TmpDir/$in-TE.sam > $OutputDir/$in-TE.bam

echo "["$(date +"%y-%m-%d %T")"] Finished running SPLITREADER beta 1.5 part 1"  | tee -a $TmpDir/log.txt  

cp $TmpDir/log.txt $OutputDir/

rm -f $readsDir/$in*.fastq

