#!/bin/bash 

#############################################
#                                           #  
#                                           #
#            SPLITREADER-sort               #
#                                           #
#                                           # 
#############################################

# usage SPLITREADER-sort.sh $filename $workspace_dir $workspace_dir/BEDfiles/SPLITREADER/$cohortname $filename $cohortname
##Questions or comments to pierre.baduel(at)ens.psl.eu

in=$1
workdir=$2 
OutputDir=$3 
popname=$4 
cohort=$5

InputDir=$workdir/$cohort/${in}/part2
# PID of this batch
IDPID=$$

# Temporary directory
TmpDir=/localtmp/PB-$IDPID
mkdir -p $TmpDir

function cleanup {
  rm -r -f $TmpDir
  rm -f $TmpDir/${in}*
  }
 trap cleanup EXIT


echo "["$(date +"%y-%m-%d %T")"] Sorting SPLITREADER output" | tee -a $InputDir/log.txt  

cp $InputDir/$in-insertion-sites.bed $OutputDir/$popname-insertion-sites.bed
cut -f 1-4 $OutputDir/$popname-insertion-sites.bed > $TmpDir/$popname-insertion-sites.cut.bed
sortBed -i $TmpDir/$popname-insertion-sites.cut.bed > $TmpDir/$popname-insertion-sites.sort1.bed
mergeBed -i $TmpDir/$popname-insertion-sites.sort1.bed > $TmpDir/$popname-insertion-sites.mrg.bed
sortBed -i $TmpDir/$popname-insertion-sites.mrg.bed > $OutputDir/$popname-insertion-sites.sort.bed 

echo "["$(date +"%y-%m-%d %T")"] Finished sorting SPLITREADER output"  | tee -a $InputDir/log.txt  

