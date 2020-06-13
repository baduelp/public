#!/bin/bash 

#############################################
#                                           #  
#               SPLITREADER                 #
#  calculate negative coverage by genome    #
#           Baduel et al. 2020              #
#                                           # 
#############################################

##Questions or comments to pbaduel(ar)bio.ens.psl.eu

in=$1 # basename of the bam file 
TE_fam=$2 # HERE BAMrc
bamdir=$3 # path to bam file
bamext=$4 # bam file extension 
popname=$5 # genome name (when different from bam file name)
depth=$6 # number of reads required to call an insertion on 1st pass
cohort=$7 # name of cohort
OutputDir=$8 # path/to/output
project_dir=$9 # path/to/files
GenomeIndexDir=${10} # path/to/genome/reference
genomefile=${11} # name of reference genome file

updownbool=1 ## calculates coverage 100bp up and down of insertion sites

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

if [ -e $OutputDir/log.txt ]
then
  echo '' > $OutputDir/log.txt  
fi

echo "["$(date +"%y-%m-%d %T")"] Calculating negative coverage over insertions " | tee -a $OutputDir/log.txt  

echo $bamdir/${in}$bamext.bam

cp $project_dir/$cohort-insertion-sites.0.bed  $OutputDir/

if [ -e $OutputDir/$popname.$cohort-insertion.0.rc ]
then
    ##remove all intervals already calculated
    bedtools subtract -a $OutputDir/$cohort-insertion-sites.0.bed -b $OutputDir/$popname.$cohort-insertion.0.rc -f 1 -r > $OutputDir/$cohort-insertion-sites.0.v2.bed 
    mv $OutputDir/$popname.$cohort-insertion.0.rc $OutputDir/$popname.$cohort-insertion.0.v1.rc 
    
    ## calculate coverage over insertions
    bam-readcount -w 1 -f $GenomeIndexDir/$genomefile $bamdir/${in}$bamext.bam -l $OutputDir/$cohort-insertion-sites.0.v2.bed > $OutputDir/$popname.$cohort-insertion.tmp.0.v2.rc
    # # process rc counts to keep only min DP over each insertion
    perl /users/a2e/pbaduel/myScripts/Process_BAMrc_splitreader.pl $cohort $depth $popname $OutputDir 0.v2
    rm $OutputDir/$popname.$cohort-insertion.tmp.0.v2.rc
    rm $OutputDir/$cohort-insertion-sites.0.v2.bed
    
    ## merge read counts with already calculated
    cat $OutputDir/$popname.$cohort-insertion.0.v*.rc > $OutputDir/$popname.$cohort-insertion.0.rc
    rm $OutputDir/$popname.$cohort-insertion.0.v*
else

  ## calculate coverage over insertions
  bam-readcount -w 1 -f $GenomeIndexDir/$genomefile $bamdir/${in}$bamext.bam -l $OutputDir/$cohort-insertion-sites.0.bed > $OutputDir/$popname.$cohort-insertion.tmp.0.rc
  # # process rc counts to keep only min DP over each insertion
  perl /users/a2e/pbaduel/myScripts/Process_BAMrc_splitreader.pl $cohort $depth $popname $OutputDir 0
  rm $OutputDir/$popname.$cohort-insertion.tmp.0.rc

fi
echo "["$(date +"%y-%m-%d %T")"] Finished local BAMrc "  | tee -a $OutputDir/log.txt  

if [ $updownbool -gt 0 ]
then
  ## calculates coverage 100bp up and down
  cp $project_dir/$cohort-insertion-sites.100up.bed  $OutputDir/
  if [ -e $OutputDir/$popname.$cohort-insertion.100up.rc ]
  then
      bedtools subtract -a $OutputDir/$cohort-insertion-sites.100up.bed -b $OutputDir/$popname.$cohort-insertion.100up.rc -f 1 -r > $OutputDir/$cohort-insertion-sites.100up.v2.bed ##remove all intervals already calculated
      mv $OutputDir/$popname.$cohort-insertion.100up.rc $OutputDir/$popname.$cohort-insertion.100up.v1.rc 
      bam-readcount -w 1 -f $GenomeIndexDir/$genomefile $bamdir/${in}$bamext.bam -l $OutputDir/$cohort-insertion-sites.100up.v2.bed > $OutputDir/$popname.$cohort-insertion.tmp.100up.v2.rc
      perl /users/a2e/pbaduel/myScripts/Process_BAMrc_splitreader.pl $cohort $depth $popname $OutputDir 100up.v2
      rm $OutputDir/$popname.$cohort-insertion.tmp.100up.v2.rc
      rm $OutputDir/$cohort-insertion-sites.100up.v2.bed
      cat $OutputDir/$popname.$cohort-insertion.100up.v*.rc > $OutputDir/$popname.$cohort-insertion.100up.rc
      rm $OutputDir/$popname.$cohort-insertion.100up.v*
  else
    bam-readcount -w 1 -f $GenomeIndexDir/$genomefile $bamdir/${in}$bamext.bam -l $OutputDir/$cohort-insertion-sites.100up.bed > $OutputDir/$popname.$cohort-insertion.tmp.100up.rc
    # # # process rc counts to keep only min DP over each insertion
    perl /users/a2e/pbaduel/myScripts/Process_BAMrc_splitreader.pl $cohort $depth $popname $OutputDir 100up
    rm $OutputDir/$popname.$cohort-insertion.tmp.100up.rc
  fi
  echo "["$(date +"%y-%m-%d %T")"] Finished 100up BAMrc "  | tee -a $OutputDir/log.txt  

  cp $project_dir/$cohort-insertion-sites.100down.bed  $OutputDir/
  if [ -e $OutputDir/$popname.$cohort-insertion.100down.rc ]
  then
      bedtools subtract -a $OutputDir/$cohort-insertion-sites.100down.bed -b $OutputDir/$popname.$cohort-insertion.100down.rc -f 1 -r > $OutputDir/$cohort-insertion-sites.100down.v2.bed ##remove all intervals already calculated
      mv $OutputDir/$popname.$cohort-insertion.100down.rc $OutputDir/$popname.$cohort-insertion.100down.v1.rc 
      bam-readcount -w 1 -f $GenomeIndexDir/$genomefile $bamdir/${in}$bamext.bam -l $OutputDir/$cohort-insertion-sites.100down.v2.bed > $OutputDir/$popname.$cohort-insertion.tmp.100down.v2.rc
      perl /users/a2e/pbaduel/myScripts/Process_BAMrc_splitreader.pl $cohort $depth $popname $OutputDir 100down.v2
      rm $OutputDir/$popname.$cohort-insertion.tmp.100down.v2.rc
      rm $OutputDir/$cohort-insertion-sites.100down.v2.bed
      cat $OutputDir/$popname.$cohort-insertion.100down.v*.rc > $OutputDir/$popname.$cohort-insertion.100down.rc
      rm $OutputDir/$popname.$cohort-insertion.100down.v*
  else
    bam-readcount -w 1 -f $GenomeIndexDir/$genomefile $bamdir/${in}$bamext.bam -l $OutputDir/$cohort-insertion-sites.100down.bed > $OutputDir/$popname.$cohort-insertion.tmp.100down.rc
    # # process rc counts to keep only min DP over each insertion
    perl /users/a2e/pbaduel/myScripts/Process_BAMrc_splitreader.pl $cohort $depth $popname $OutputDir 100down
    rm $OutputDir/$popname.$cohort-insertion.tmp.100down.rc
  fi
  echo "["$(date +"%y-%m-%d %T")"] Finished 100down BAMrc "  | tee -a $OutputDir/log.txt  
fi


echo "["$(date +"%y-%m-%d %T")"] Finished BAMrc "  | tee -a $OutputDir/log.txt  
