#!/bin/bash 

#############################################
#                                           #  
#                                           #
#            SPLITREADER-sort               #
#                                           #
#                                           # 
#############################################

##Questions or comments to pierre.baduel(at)ens.psl.eu

#######################################################################################################################################################################################
# # ARGUMENTS PART

# Verifying of the arguments
if [ "$#" -ne 6 ]; then
  echo ""
  echo "Name : $0"
  echo ""
  echo "Usage: "
  echo "	$0 <in> <workdir> <OutputDir> <popname> <cohort> <TmpDir>"
  echo ""
  echo "Description:"
  echo ""
  echo "	Sorts all the output bed file of SPLITREADER-beta1.5_part2.sh to BEDfiles/SPLITREADER/\$cohortname/ALL"
  echo ""
  echo "Arguments:"
  echo ""
  echo "	<in>:  Basename of the BAM file"
  echo "	<workdir>:   	Workspace Directory"
  echo "	<OutputDir>:   	path/to/output"
  echo "	<popname>:   	Genome name (when different from bam file name)"
  echo "	<cohort>:       Name of Cohort"
  echo "	<TmpDir>:		Tmp Dir"
  echo ""
  exit 1
fi

in=$1
workdir=$2 
OutputDir=$3 
popname=$4 
cohort=$5
tmp_dir=$6

# PID of this batch
IDPID=$$

# Temporary directory
mkdir -p $tmp_dir/$IDPID
TmpDir=$tmp_dir/$IDPID

#############################


# Verifying of the arguments

if [ ! -d "$workdir" ]; then
  echo ""
  echo "Error: Workspace directory not found!"
  echo ""
  exit 1
fi

if [ ! -d "$workdir/$cohort/$in" ]; then
  echo ""
  echo "Error: sample $in not analysed!"
  echo ""
  exit 1
fi

if [ ! -d "$OutputDir" ]; then
  echo "Error: OutputDir not found! (should be: \$workdir/BEDfiles/SPLITREADER/\$cohortname)"
  exit 1
fi

if [ ! -d "$workdir/$cohort/$popname" ]; then
  echo ""
  echo "Error: sample $popname not analysed!"
  echo ""
  exit 1
fi

if [ -z "$cohort" ]; then
  echo "Error: Cohort name not specified!"
  exit 1
fi


#######################################################################################################################################################################################

InputDir=$workdir/$cohort/${in}/part2


function cleanup {
  rm -r -f $TmpDir
  rm -f $TmpDir/${in}*
  }
 trap cleanup EXIT

#######################################################################################################################################################################################
	# # Sorts the bed file from SPLITREADER-beta1.5_part2.sh and put it sorted in $OutputDir/ALL

echo "["$(date +"%y-%m-%d %T")"] Sorting SPLITREADER output for ${in}" | tee -a $InputDir/log.txt  

cp $InputDir/$in-insertion-sites.bed $OutputDir/ALL/$popname-insertion-sites.bed
cut -f 1-4 $OutputDir/ALL/$popname-insertion-sites.bed > $TmpDir/$popname-insertion-sites.cut.bed
sortBed -i $TmpDir/$popname-insertion-sites.cut.bed > $TmpDir/$popname-insertion-sites.sort1.bed
mergeBed -i $TmpDir/$popname-insertion-sites.sort1.bed > $TmpDir/$popname-insertion-sites.mrg.bed
sortBed -i $TmpDir/$popname-insertion-sites.mrg.bed > $OutputDir/ALL/$popname-insertion-sites.sort.bed 

echo "["$(date +"%y-%m-%d %T")"] Finished sorting SPLITREADER output for ${in}"  | tee -a $InputDir/log.txt  

