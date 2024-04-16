#!/bin/bash 

#############################################
#                                           #  
#               SPLITREADER                 #
#             beta1.5 - part 2              #
#                                           # 
#############################################

##Questions or comments to quadrana(ar)biologie.ens.fr

#######################################################################################################################################################################################
# # ARGUMENTS PART

# Verifying of the arguments
if [ "$#" -ne 10 ] && [ "$#" -ne 9 ] && [ "$#" -ne 8 ] && [ "$#" -ne 7 ]; then
  echo ""
  echo "Name : $0"
  echo ""
  echo "Usage: "
  echo "	$0 <bam file name> <cohort name> <workdir> <genome file> <TE annotation file> <max coverage> <TmpDir> <LS: library size> <number of threads> <Length>"
  echo ""
  echo "Description:"
  echo ""
  echo "	Second Part of SPLITREADER, a pipeline that detects Transposable Elements (TE) insertions, from the alignment of samples on their reference genome, based on their short reads paired-end sequences."
  echo "	This second part consists of mapping the unmapped reads from Part1 and map it to the given TE library."
  echo ""
  echo "Arguments:"
  echo ""
  echo "	<bam file name>:        Name of the BAM file"
  echo "	<cohort name>:   	Name of the cohort of genomes analyzed"
  echo "	<workdir>:   		Working directory"
  echo "	<genome file>:          Genome bowtie2 index"
  echo "	<TE annotation file>:   File of the annotation of the TEs"
  echo "	<max coverage>:         Maximum coverage"
  echo "	<tmp_dir>:				Tmp dir"
  echo "	<LS: library size>:     Library size (optional : default 300)"
  echo "	<number of threads>:    Number of threads for processing (optional : default 2)"
  echo "	<Length>:           	Minimum read length (If not specified, the program will calculate the longest read)"
  echo ""
  exit 1
fi

# Arguments : 

LS_DEFAULT=300
CORES_DEFAULT=2
pe=TRUE # Paired-end ?


in=$1 # bam file name
cohort=$2 # name of cohort of genomes analyzed
workdir=$3 # working directory
GenomeFile=$4 # genome
TEannotfile=$5 # gff or gff3 file
maxcov=$6 # max coverage
tmp_dir=$7 # tmp dir
LS=${8:-$LS_DEFAULT} # Library Size
mult_mapping_limit=4 # number of multiple mapping allowed before reads are rejected (important when TEs are too similar to each other)
#### If not specified, the program estimate them. To this end, the program calculates the minimum number of reads as 3 standard deviation under the mean whole genome coverage
####This value should be at least 3, if not, it is forced to be 3

CORES=${9:-$CORES_DEFAULT} # Number of threads : By default is 2 threads

LENGTH=${10}
        
# PID of this batch
IDPID=$$

# Temporary directory
mkdir -p $tmp_dir/$IDPID                                                  
TmpDir=$tmp_dir/$IDPID 
  
#############################


# Verifying of the arguments

if [ ! -d "$workdir/$cohort/${in}/part1" ]; then
  echo ""
  echo "Error: Directory of the output of Part1 not found!"
  echo ""
  exit 1
fi

if [ ! -d "$TmpDir" ]; then
  echo ""
  echo "Error: Workspace directory not found!"
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

if [ ! -f "$workdir/Reference/$GenomeFile.rev.2.bt2" ]; then
  echo ""
  echo "Error: Genome bowtie2 index not found!"
  echo "Indexing ..."
  echo ""
  if [ -f "$workdir/Reference/$GenomeFile.fasta" ]; then
    bowtie2-build "$workdir/Reference/$GenomeFile.fasta" "$workdir/Reference/$GenomeFile"
  else
    echo "Fasta file of the Genome not found!"
    exit 1
  fi
  echo ""
fi


if ! [[ "$LS" =~ ^[0-9]+$ ]]; then
  echo ""
  echo "Error: Library size must be a positive integer."
  echo ""
  exit 1
fi

if ! [[ "$CORES" =~ ^[0-9]+$ ]]; then
  echo ""
  echo "Error: Number of threads must be a positive integer."
  echo ""
  exit 1
fi


#######################################################################################################################################################################################
# # Verifications


InputDir=$workdir/$cohort/${in}/part1 #TE-sam from part1 

OutputDir=$workdir/$cohort/${in}/part2

if [ ! -e $OutputDir/$in-insertion-sites.bed ]
then
  echo '' > $OutputDir/$in-insertion-sites.bed
fi
#path to list of TEs to analyze (TE-information-all.txt)
#TE-> should be indicated in the first column of the TE-information.txt file located in listDir
#TSD -> should be indicated in the second column of the TE-information.txt file located in listDir
listDir=$workdir/TE_sequence
TE_info='TE-information-all'

TEannot=$workdir/TE_sequence/$TEannotfile.gff

#Bowtie2 index for reference genome
GenomeIndexFile=$workdir/Reference/$GenomeFile

# Temporary directory
##
if [ ! -d "$OutputDir" ]; then
  mkdir -p $OutputDir
fi

function cleanup {
  rm -rf $TmpDir/*
  rmdir $TmpDir
  # rm -f $TmpResultsDir/*
  # rm -r -f $TmpResultsDir
  }
 trap cleanup EXIT
  
############################################################################################################

echo "["$(date +"%y-%m-%d %T")"] Running SPLITREADER part 2 for ${in}" | tee -a $OutputDir/log.txt

if [ ! -e $InputDir/$in-TE.bam ]
then  
  if [ -e $InputDir/$in-TE.sam ]
  then
    echo "convert sam to bam" | tee -a $OutputDir/log.txt
     samtools view -Sb $InputDir/$in-TE.sam > $InputDir/$in-TE.bam
     rm $InputDir/$in-TE.sam
  fi
fi
  
# # Extracting superfamily reads from remapped TE sam
echo "["$(date +"%y-%m-%d %T")"] Extracting superfamily reads from ${in}" | tee -a $OutputDir/log.txt

end=`wc -l $listDir/$TE_info.txt | awk '{print $1}'`

############################################################################################################

# # Starting the SPLITREADER pipeline for each TE in the TE-information.txt file

for ((l=1; $l<=$end; l=$l+1)); do

  STARTTIME=$(date +%s)

  TE=`sed -n "${l}p" $listDir/$TE_info.txt | awk '{print $1}'`
  TSD=`sed -n "${l}p" $listDir/$TE_info.txt | awk '{print $2}'`
  echo -e "\n"

  if [ ! -s $OutputDir/$in-$TE-split.bam ]; then
    TmpResultsDir=$OutputDir/$TE
    mkdir -p $TmpResultsDir
    echo "["$(date +"%y-%m-%d %T")"] ##### RUNNING SPLIT-READ ANALYSIS ON $TE (TSD size = $TSD bp)######"  | tee -a $OutputDir/log.txt  
    echo ""
    ############# 

    # Selecting split-reads by mapping the unmapped reads over TE extremities
      
    #############################################################
      ###filter soft-clipped reads with at least 20nt softclipped at 5' or 3' read's end 
    echo "["$(date +"%y-%m-%d %T")"] Selecting split-reads"  | tee -a $OutputDir/log.txt    
    samtools view -F 4 $InputDir/$in-TE.bam | grep -w $TE | awk '$6~/^[2-8][0-9]S/ || $6~/[2-8][0-9]S$/ || $6~/^1[0-9][0-9]S/ || $6~/1[0-9][0-9]S$/ {print $1"\t"$10"\t"$11}' | sort -k1,1 -k2,2 -u | awk '{print "@"$1"\n"$2"\n+\n"$3}' > $TmpResultsDir/$in-$TE-split.fastq 2>> $OutputDir/log.txt

    #############################################################
      ###filtering pair ends of unmapped reads  
    echo "["$(date +"%y-%m-%d %T")"] Selecting discordant-reads"  | tee -a $OutputDir/log.txt

    samtools view -F 8 $InputDir/$in-TE.bam | grep -w $TE | awk '{print $1}' >> $TmpResultsDir/reads.name.disc  2>> $OutputDir/log.txt
    
    Nsplitreads=`wc -l $TmpResultsDir/$in-$TE-split.fastq | awk '{print $1/4}'` 
    Ndiscreads=`sort -k1,1 $TmpResultsDir/reads.name.disc | uniq | wc -l | awk '{print $1}'`  
    
    echo "$Nsplitreads splitreads and $Ndiscreads discordant reads identified on $TE" | tee -a $OutputDir/log.txt  

    if [ $(($Nsplitreads+$Ndiscreads)) -eq 0 ]
    then
      echo "No reads identified: end analysis on $TE" | tee -a $OutputDir/log.txt 

    else
      
      if [ $Ndiscreads -gt 0 ]
      then

        # grab first in pair
        samtools view -f 64 -u $InputDir/$in-TE.bam > $TmpResultsDir/$in-TE-disc-first.bam 2>> $OutputDir/log.txt
        # grab second in pair
        samtools view -f 128 -u $InputDir/$in-TE.bam > $TmpResultsDir/$in-TE-disc-second.bam 2>> $OutputDir/log.txt
        

        _JAVA_OPTIONS=" -Xmx15g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp" picard FilterSamReads -INPUT $TmpResultsDir/$in-TE-disc-first.bam -FILTER includeReadList -READ_LIST_FILE $TmpResultsDir/reads.name.disc -OUTPUT $TmpResultsDir/$in-$TE-selected-disc-first.sam -TMP_DIR $TmpDir/javatemp  2>> $OutputDir/log.txt


        _JAVA_OPTIONS=" -Xmx15g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp" picard FilterSamReads -INPUT $TmpResultsDir/$in-TE-disc-second.bam -FILTER includeReadList -READ_LIST_FILE $TmpResultsDir/reads.name.disc -OUTPUT $TmpResultsDir/$in-$TE-selected-disc-second.sam -TMP_DIR $TmpDir/javatemp  2>> $OutputDir/log.txt
        
        
        cat $TmpResultsDir/$in-$TE-selected-disc-first.sam | awk '$1!~/^@/ {print $1"\t"$10"\t"$11}' | sort -u -k1,1 -k2,2 |  awk '{print "@"$1"|1\n"$2"\n+\n"$3}' > $TmpResultsDir/$in-$TE-disc.fastq 2>> $OutputDir/log.txt
        
        cat $TmpResultsDir/$in-$TE-selected-disc-second.sam | awk '$1!~/^@/ {print $1"\t"$10"\t"$11}' | sort -u -k1,1 -k2,2 | awk '{print "@"$1"|2\n"$2"\n+\n"$3}' >> $TmpResultsDir/$in-$TE-disc.fastq 2>> $OutputDir/log.txt
      fi
      

      
      rm -f $TmpResultsDir/reads.name
        
      ################################
    
    
      ###Estimating max read size (If necessary)
      
      if [ -z "$LENGTH" ]
      then
        LENGTH=`head -5000 $TmpResultsDir/$in-$TE-split.fastq | awk 'NR%4 == 2 {print length($0)}' | sort | tail -1 `  
        length=$((LENGTH-20))
        echo "["$(date +"%y-%m-%d %T")"] Maximum Read length: $LENGTH [Estimated] " | tee -a $OutputDir/log.txt
        else
        length=$((LENGTH-20))
        echo "["$(date +"%y-%m-%d %T")"] Maximum Read length: $LENGTH [User defined] " | tee -a $OutputDir/log.txt
      fi
      
          
      ###Recursive split-reads mapping
      # step 1 for 3' read extremity: begining the loop.
      
        
      echo "["$(date +"%y-%m-%d %T")"] Analyzing split-reads" | tee -a $OutputDir/log.txt
      
      
      if [ $Nsplitreads -gt 0 ]
      then
        bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$in-$TE-split.fastq -S $TmpResultsDir/$in-$TE-local.sam --local --very-sensitive --threads $CORES --quiet 

        samtools view -H -S $TmpResultsDir/$in-$TE-local.sam > $TmpResultsDir/$in-$TE-split-local-up.sam 
        cat $TmpResultsDir/$in-$TE-split-local-up.sam > $TmpResultsDir/$in-$TE-split-local-down.sam 
        
        ##extracting split reads upstream insertion
        samtools view -F 4 -q 5 -S $TmpResultsDir/$in-$TE-local.sam | awk '$6~/^[0-9][0-9]S/ || $6~/^1[0-9][0-9]S/ {print $0}' >> $TmpResultsDir/$in-$TE-split-local-down.sam 
        samtools view -F 4 -q 5 -S $TmpResultsDir/$in-$TE-local.sam | awk '$6~/[0-9][0-9]S$/ || $6~/1[0-9][0-9]S$/ {print $0}' >> $TmpResultsDir/$in-$TE-split-local-up.sam 

        samtools view -Sbu $TmpResultsDir/$in-$TE-split-local-down.sam | samtools sort - -o $TmpResultsDir/$in-$TE-splitjunction-down.bam 
        samtools view -Sbu $TmpResultsDir/$in-$TE-split-local-up.sam | samtools sort - -o $TmpResultsDir/$in-$TE-splitjunction-up.bam 
          
        ############################################
        
        ##Refining insertion sites

        length=$(($((LENGTH/2))-1))
        echo "["$(date +"%y-%m-%d %T")"] Refining insertion sites" | tee -a $OutputDir/log.txt
        echo -e "\n"
        echo -n "Progression: ["
        
        bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$in-$TE-split.fastq -S $TmpResultsDir/$in-$TE-splitjunction-5-$length.sam --un $TmpResultsDir/$in-$TE-split-5-$length -5 $length --very-sensitive --threads $CORES --local --quiet 2>> $OutputDir/log.txt
      
        samtools view -H $TmpResultsDir/$in-$TE-splitjunction-5-$length.sam > $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam 2>> $OutputDir/log.txt
        cat $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam > $TmpResultsDir/$in-$TE-splitjunction-5-$length-up.sam 2>> $OutputDir/log.txt
        cat $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam > $TmpResultsDir/$in-$TE-splitjunction-3-$length-down.sam 2>> $OutputDir/log.txt
        cat $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam > $TmpResultsDir/$in-$TE-splitjunction-3-$length-up.sam 2>> $OutputDir/log.txt
        
        # reads S at beginning and M at end
        samtools view -F 4 -q 5 -S $TmpResultsDir/$in-$TE-splitjunction-5-$length.sam | awk '($6~/^[0-9][0-9]S/ || $6~/^1[0-9][0-9]S/) && ($6~/[2-9][0-9]M$/) {print $0}' >> $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam 2>> $OutputDir/log.txt

        # reads M at beginning and S at end
        samtools view -F 4 -q 5 -S $TmpResultsDir/$in-$TE-splitjunction-5-$length.sam | awk '($6~/[0-9][0-9]S$/ || $6~/1[0-9][0-9]S$/) && ($6~/^[2-9][0-9]M/)  {print $0}' >> $TmpResultsDir/$in-$TE-splitjunction-5-$length-up.sam 2>> $OutputDir/log.txt

        samtools view -Sbu $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam  | samtools sort - -o $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.bam 2>> $OutputDir/log.txt
        samtools view -Sbu $TmpResultsDir/$in-$TE-splitjunction-5-$length-up.sam  | samtools sort - -o $TmpResultsDir/$in-$TE-splitjunction-5-$length-up.bam 2>> $OutputDir/log.txt
      
        bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$in-$TE-split.fastq -S $TmpResultsDir/$in-$TE-splitjunction-3-$length.sam --un $TmpResultsDir/$in-$TE-split-3-$length -3 $length --very-sensitive --threads $CORES --local --quiet 2>> $OutputDir/log.txt

        samtools view -F 4 -q 5 -S $TmpResultsDir/$in-$TE-splitjunction-3-$length.sam | awk '$6~/^[0-9][0-9]S/ || $6~/^1[0-9][0-9]S/ {print $0}' >> $TmpResultsDir/$in-$TE-splitjunction-3-$length-down.sam 2>> $OutputDir/log.txt
        samtools view -F 4 -q 5 -S $TmpResultsDir/$in-$TE-splitjunction-3-$length.sam | awk '$6~/[0-9][0-9]S$/ || $6~/^1[0-9][0-9]S/ {print $0}' >> $TmpResultsDir/$in-$TE-splitjunction-3-$length-up.sam 2>> $OutputDir/log.txt

        samtools view -Sbu $TmpResultsDir/$in-$TE-splitjunction-3-$length-down.sam  | samtools sort - -o $TmpResultsDir/$in-$TE-splitjunction-3-$length-down.bam 2>> $OutputDir/log.txt
        samtools view -Sbu $TmpResultsDir/$in-$TE-splitjunction-3-$length-up.sam  | samtools sort - -o $TmpResultsDir/$in-$TE-splitjunction-3-$length-up.bam 2>> $OutputDir/log.txt 
        
        echo -n "]"
        echo -e "\n"
      fi

      ###########################
      if [ -n "$pe" ]
      then
        if [ $Ndiscreads -gt 0 ]
        then
          echo "["$(date +"%y-%m-%d %T")"] Analyzing discordant reads" | tee -a $OutputDir/log.txt
          
          # increase mapping penalities for mismatches and indels to prevent mismapped discordant reads to break-down the clusters
          bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$in-$TE-disc.fastq -S $TmpResultsDir/$in-$TE-endtoend.sam --very-sensitive --threads $CORES --quiet  --mp 13 --rdg 8,5 --rfg 8,5
        
          samtools view -Sbu -f 16 -q 5 $TmpResultsDir/$in-$TE-endtoend.sam | samtools sort - -o $TmpResultsDir/$in-$TE-disc-down.bam 
          samtools view -Sbu -F 16 -q 5 $TmpResultsDir/$in-$TE-endtoend.sam | samtools sort - -o $TmpResultsDir/$in-$TE-disc-up.bam 
        
          rm -f $TmpResultsDir/$in-$TE-local.sam
          rm -f $TmpResultsDir/$in-$TE-endtoend.sam
        fi
      fi


      ############################################
      # Post-treatment:
      
      echo "Merge the 5' and 3' clusters to create the downstream and upstream cluster" | tee -a $OutputDir/log.txt

      samtools merge -f -u $TmpResultsDir/$in-$TE-up.bam $TmpResultsDir/$in-$TE-splitjunction-5-$length-up.bam $TmpResultsDir/$in-$TE-splitjunction-3-$length-up.bam $TmpResultsDir/$in-$TE-splitjunction-up.bam  2>> $OutputDir/log.txt

      samtools merge -f -u $TmpResultsDir/$in-$TE-down.bam $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.bam $TmpResultsDir/$in-$TE-splitjunction-3-$length-down.bam $TmpResultsDir/$in-$TE-splitjunction-down.bam  2>> $OutputDir/log.txt

      samtools sort $TmpResultsDir/$in-$TE-down.bam -o $TmpResultsDir/$in-$TE-down 2>> $OutputDir/log.txt
      samtools sort $TmpResultsDir/$in-$TE-up.bam -o $TmpResultsDir/$in-$TE-up 2>> $OutputDir/log.txt
    
    

      echo "Calculate the coverage over mapped regions - filter regions according to minimum and maximum read-depth" | tee -a $OutputDir/log.txt

      samtools depth $TmpResultsDir/$in-$TE-up.bam | awk -v M=$maxcov '$3<(M) {print $1 "\t" $2 "\t"$2"\t"$3}' | sort -k 1,1 -k2,2n > $TmpResultsDir/$in-$TE-up.bed 
      if [ -s $TmpResultsDir/$in-$TE-up.bed ]; then
        mergeBed -i $TmpResultsDir/$in-$TE-up.bed  -c 4 -o max | awk -v l=$LENGTH '($3-$2)<=l {print $0}' > $TmpResultsDir/$in-$TE-up-merge.bed 2>> $OutputDir/log.txt
      else
        echo '' > $TmpResultsDir/$in-$TE-up-merge.bed
      fi
      
      samtools depth $TmpResultsDir/$in-$TE-down.bam | awk -v M=$maxcov ' $3<(M) {print $1 "\t" $2 "\t"$2"\t"$3}' | sort -k 1,1 -k2,2n > $TmpResultsDir/$in-$TE-down.bed 
      if [ -s $TmpResultsDir/$in-$TE-down.bed ]; then
        mergeBed -i $TmpResultsDir/$in-$TE-down.bed -c 4 -o max | awk -v l=$LENGTH '($3-$2)<=l {print $0}' > $TmpResultsDir/$in-$TE-down-merge.bed 2>> $OutputDir/log.txt
      else
        echo '' > $TmpResultsDir/$in-$TE-down-merge.bed
      fi
    
      teid=`echo $TE | sed 's/\@/\t/' | awk '{print $1}'`
      grep -w $teid $TEannot | awk '{print $1"\t"$4"\t"$5}' > $TmpResultsDir/disc-excluding.tmp

      echo "merge cluster of discordant-reads" | tee -a $OutputDir/log.txt
      if [ -n "$pe" ]
      then
        echo "["$(date +"%y-%m-%d %T")"] Searching for discordant-reads clusters..." | tee -a $OutputDir/log.txt
        
        
        samtools depth $TmpResultsDir/$in-$TE-disc-down.bam | awk -v M=$maxcov '$3<(M) {print $1 "\t" $2 "\t"$2 "\t"$3}' | sort -k1,1 -k2,2n > $TmpResultsDir/$in-$TE-disc-down.bed 
        if [ -s $TmpResultsDir/$in-$TE-disc-down.bed ]; then
          mergeBed -i $TmpResultsDir/$in-$TE-disc-down.bed  -d $LS -c 4 -o max | sort -k1,1 -k2,2n > $TmpResultsDir/$in-$TE-disc-down.tmp  2>> $OutputDir/log.txt
        else
          echo '' > $TmpResultsDir/$in-$TE-disc-down.tmp
        fi

        samtools depth $TmpResultsDir/$in-$TE-disc-up.bam | awk -v M=$maxcov '$3<(M) {print $1 "\t" $2 "\t"$2 "\t"$3}' | sort -k1,1 -k2,2n > $TmpResultsDir/$in-$TE-disc-up.bed 
        if [ -s $TmpResultsDir/$in-$TE-disc-up.bed ]; then
          mergeBed -i $TmpResultsDir/$in-$TE-disc-up.bed  -d $LS -c 4 -o max | sort -k1,1 -k2,2n > $TmpResultsDir/$in-$TE-disc-up.tmp 2>> $OutputDir/log.txt
        else
          echo '' > $TmpResultsDir/$in-$TE-disc-up.tmp
        fi
        
        #filtering discordat-read clusters that overlap by a legth >than TSD+15bp
        intersectBed -a $TmpResultsDir/$in-$TE-disc-up.tmp -b $TmpResultsDir/$in-$TE-disc-down.tmp -wa -wb | awk -v l=$LENGTH '($3-$6)>l {print $1"\t"$2"\t"$7}' >> $TmpResultsDir/disc-excluding.tmp  2>> $OutputDir/log.txt

        sort -o $TmpResultsDir/disc-excluding.tmp -k1,1 -k2,2n $TmpResultsDir/disc-excluding.tmp  2>> $OutputDir/log.txt

        # $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-disc-up.tmp -b $TmpResultsDir/$in-$TE-disc-down.tmp -wa -wb | awk -v l=$LENGTH '($3-$6)>l {print $1"\t"$2"\t"$7}' >> $TmpResultsDir/disc-excluding.tmp 2>> $OutputDir/log.txt
        
        intersectBed -a $TmpResultsDir/$in-$TE-disc-up.tmp -b $TmpResultsDir/disc-excluding.tmp	 -sorted -v | awk '{print $1 "\t" $2 "\t"$3}' > $TmpResultsDir/$in-$TE-disc-up.tmp2 2>> $OutputDir/log.txt
        
        intersectBed -a $TmpResultsDir/$in-$TE-disc-down.tmp -b $TmpResultsDir/disc-excluding.tmp -sorted -v | awk '{print $1 "\t" $2 "\t"$3}' > $TmpResultsDir/$in-$TE-disc-down.tmp2 2>> $OutputDir/log.txt
        
        intersectBed -b $TmpResultsDir/$in-$TE-disc-up.bam -a $TmpResultsDir/$in-$TE-disc-up.tmp2 -bed -c | sort -k1,1 -k2,2n > $TmpResultsDir/$in-$TE-disc-up.bed 2>> $OutputDir/log.txt

        intersectBed -b $TmpResultsDir/$in-$TE-disc-down.bam -a $TmpResultsDir/$in-$TE-disc-down.tmp2 -bed -c  | sort -k1,1 -k2,2n > $TmpResultsDir/$in-$TE-disc-down.bed 2>> $OutputDir/log.txt
        
        rm -f $TmpResultsDir/$in-$TE-disc-down.tmp*
        rm -f $TmpResultsDir/$in-$TE-disc-up.tmp*
      fi 

      # echo "Defining TE insertions based on splitreads only: searching for overlapping clusters meeting the expected TSD size and number of supporting reads
      echo "["$(date +"%y-%m-%d %T")"] Defining insertions..." | tee -a $OutputDir/log.txt

      intersectBed -a $TmpResultsDir/$in-$TE-up-merge.bed -b $TmpResultsDir/$in-$TE-down-merge.bed -wo | awk -v tsd=$TSD -v te=$TE '($6-$2)>10 && ($7-$3)>10 && $9>=tsd && $9<(tsd+15) {print $1 "\t" $6 "\t" $3 "\t"te"\t" ($9-1)"\t"$4"\t"$8}' | intersectBed -a stdin -b $TmpResultsDir/disc-excluding.tmp -v >> $TmpResultsDir/$in-insertion-sites.1 2>> $OutputDir/log.txt
      
      if [ -s $TmpResultsDir/$in-insertion-sites.1 ]
      then
      awk -v l=$LS '{print $1"\t"$2"\t"$3+l"\t"$4}' $TmpResultsDir/$in-$TE-disc-up.bed | intersectBed -a $TmpResultsDir/$in-insertion-sites.1 -b stdin -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11}' > $TmpResultsDir/$in-insertion-sites.2  2>> $OutputDir/log.txt
      awk -v l=$LS '{print $1"\t"$2"\t"$3+l"\t"$4}' $TmpResultsDir/$in-$TE-disc-up.bed | intersectBed -a $TmpResultsDir/$in-insertion-sites.1 -b stdin -v | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""0"}' >> $TmpResultsDir/$in-insertion-sites.2  2>> $OutputDir/log.txt
      
      awk -v l=$LS ' $2-l>0 {print $1"\t"$2-l"\t"$3"\t"$4}' $TmpResultsDir/$in-$TE-disc-down.bed | intersectBed -a $TmpResultsDir/$in-insertion-sites.2 -b stdin -wa -wb | awk '($6+$7)>=r {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$12}' >> $OutputDir/$in-insertion-sites.bed  2>> $OutputDir/log.txt
      
      # -v r=$READS '($6+$7)>=r (non informative filter)
      awk -v l=$LS '$2-l>0 {print $1"\t"$2-l"\t"$3"\t"$4}' $TmpResultsDir/$in-$TE-disc-down.bed | intersectBed -a $TmpResultsDir/$in-insertion-sites.2 -b stdin -v | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t""0"}' >> $OutputDir/$in-insertion-sites.bed  2>> $OutputDir/log.txt
      
      else
      echo "No split-reads clusters identified... skipping this step" | tee -a $OutputDir/log.txt
      fi
      
      ################################
        
      #Defining TE insertions based in discordant reads and supported for at least one split-read:
      if [ -n "$pe" ]
        then
        echo "selecting cluster of discordant reads that were not called by split-reads"  | tee -a $OutputDir/log.txt
        echo -n "."  >> $OutputDir/log.txt
        awk -v l=$LS '{print $1"\t"$2"\t"$3+l"\t"$4}' $TmpResultsDir/$in-$TE-disc-up.bed | intersectBed -b $TmpResultsDir/$in-insertion-sites.1 -a stdin -v > $TmpResultsDir/$in-$TE-disc-up-outer.bed 2>> $OutputDir/log.txt
        echo -n "."  >> $OutputDir/log.txt
        awk -v l=$LS '$2>l {print $1"\t"$2-l"\t"$3"\t"$4}' $TmpResultsDir/$in-$TE-disc-down.bed | intersectBed -b $TmpResultsDir/$in-insertion-sites.1 -a stdin -v > $TmpResultsDir/$in-$TE-disc-down-outer.bed 2>> $OutputDir/log.txt
        echo -n "."  >> $OutputDir/log.txt
        #selecting overlapping upstream and downstream clusters of discordant reads   
        intersectBed -a $TmpResultsDir/$in-$TE-disc-up-outer.bed -b $TmpResultsDir/$in-$TE-disc-down-outer.bed -wa -wb | awk -v l=$LS '$3>l {print $1"\t"$2"\t"$7"\t"$4"\t"$8"\t"$3-l"\t"$6+l}' > $TmpResultsDir/$in-$TE-disc-cluster.0 2>> $OutputDir/log.txt
        echo -n "."  >> $OutputDir/log.txt
        awk '$6<=($7+8) {print $1"\t"$6-4"\t"$7+4"\t"$4"\t"$5"\t"$2"\t"$3}' $TmpResultsDir/$in-$TE-disc-cluster.0 > $TmpResultsDir/$in-$TE-disc-cluster.1 2>> $OutputDir/log.txt
        echo -n "."  >> $OutputDir/log.txt
        awk '$6>($7+8) {print $1"\t"$7-4"\t"$6+4"\t"$4"\t"$5"\t"$2"\t"$3}' $TmpResultsDir/$in-$TE-disc-cluster.0 >> $TmpResultsDir/$in-$TE-disc-cluster.1 2>> $OutputDir/log.txt
        echo -n "."  >> $OutputDir/log.txt
        if [[ -s $TmpResultsDir/$in-$TE-disc-cluster.1 ]] ; then
          echo -n "."  >> $OutputDir/log.txt
          intersectBed -a $TmpResultsDir/$in-$TE-disc-cluster.1 -b $TmpResultsDir/$in-$TE-up-merge.bed -wa -wb | awk '{print $1"\t"$2"\t"$10"\t"$4"\t"$5"\t"$11"\t"$6"\t"$7}' | sort -k1,1 -k2,2n | mergeBed -i stdin -c 4,5,6,7,8 -o max,max,max,max,max > $TmpResultsDir/$in-$TE-disc-cluster.3 2>> $OutputDir/log.txt
          echo -n "."  >> $OutputDir/log.txt
          intersectBed -a $TmpResultsDir/$in-$TE-disc-cluster.1 -b $TmpResultsDir/$in-$TE-up-merge.bed -v | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""0""\t"$6"\t"$7}'  >> $TmpResultsDir/$in-$TE-disc-cluster.3 2>> $OutputDir/log.txt
          echo -n "."  >> $OutputDir/log.txt
          intersectBed -a $TmpResultsDir/$in-$TE-disc-cluster.3 -b $TmpResultsDir/$in-$TE-down-merge.bed -wa -wb | awk -v te=$TE '{print $1"\t"$10"\t"$3"\t"te"\t""NA""\t"$6"\t"$12"\t"$4"\t"$5"\t"$7"\t"$8}' | sort -k1,1 -k2,2n | mergeBed -i stdin -c 4,5,6,7,8,9,10,11 -o distinct,distinct,distinct,max,max,max,min,max >> $OutputDir/$in-insertion-sites.bed 2>> $OutputDir/log.txt
          echo -n "."  >> $OutputDir/log.txt
          intersectBed -a $TmpResultsDir/$in-$TE-disc-cluster.3 -b $TmpResultsDir/$in-$TE-down-merge.bed -v | awk -v te=$TE '{print $1"\t"$2"\t"$3"\t"te"\t""NA""\t"$6"\t""0""\t"$4"\t"$5"\t"$7"\t"$8}' | sort -k1,1 -k2,2n | mergeBed -i stdin -c 4,5,6,7,8,9,10,11 -o distinct,distinct,distinct,max,max,max,min,max >> $OutputDir/$in-insertion-sites.bed 2>> $OutputDir/log.txt
          echo -n "."  >> $OutputDir/log.txt
        else

          echo "No discording-reads clusters identified... skiping this step" | tee -a $OutputDir/log.txt
        fi
      fi

      INSERTIONS=`grep -w $TE $OutputDir/$in-insertion-sites.bed | wc -l | awk '{print $1}'` 
      
      echo "["$(date +"%y-%m-%d %T")"] Split-read analyis done: $INSERTIONS putative insertions identified..." | tee -a $OutputDir/log.txt

      
      # if [ $INSERTIONS -gt 0 ]
      # then
        # # merging bam files and moving them to the output folder 
        samtools merge -f $TmpResultsDir/$in-$TE-split.bam $TmpResultsDir/$in-$TE-up.bam $TmpResultsDir/$in-$TE-down.bam $TmpResultsDir/$in-$TE-disc-down.bam $TmpResultsDir/$in-$TE-disc-up.bam 
        samtools sort $TmpResultsDir/$in-$TE-split.bam -o $TmpResultsDir/$in-$TE-split-sorted.bam 
        samtools view -H $TmpResultsDir/$in-$TE-split-sorted.bam | grep ^@SQ > $TmpResultsDir/$in-$TE-split-sorted.sam 
        samtools view -H $TmpResultsDir/$in-$TE-split-sorted.bam | grep ^@PG | head -1 >> $TmpResultsDir/$in-$TE-split-sorted.sam 
        samtools view $TmpResultsDir/$in-$TE-split-sorted.bam >> $TmpResultsDir/$in-$TE-split-sorted.sam 
        samtools view -hb $TmpResultsDir/$in-$TE-split-sorted.sam > $OutputDir/$in-$TE-split.bam 
        samtools index $OutputDir/$in-$TE-split.bam 
      # fi
      ENDTIME=$(date +%s)
      echo "["$(date +"%y-%m-%d %T")"] It takes $((ENDTIME-STARTTIME)) seconds to analyse $TE" | tee -a $OutputDir/log.txt
    fi
      
    rm -rf $TmpResultsDir/*
    rmdir $TmpResultsDir
  
  else
  
    echo -e "["$(date +"%y-%m-%d %T")"] ##### FOUND EXISTING SPLIT-BAM ON $TE (TSD size = $TSD bp)######"  | tee -a $OutputDir/log.txt  
    echo ""
  fi
done
echo "["$(date +"%y-%m-%d %T")"] Finished running SPLITREADER part 2 for ${in}" | tee -a $OutputDir/log.txt
# cp $OutputDir/log.txt $OutputDir/
