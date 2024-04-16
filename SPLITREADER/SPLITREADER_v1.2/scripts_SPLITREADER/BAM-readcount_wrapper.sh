#!/bin/bash 

#############################################
#                                           #  
#               SPLITREADER                 #
#  calculate negative coverage by genome    #
#           Baduel et al. 2020              #
#                                           # 
#############################################

##Questions or comments to pbaduel(ar)bio.ens.psl.eu

#######################################################################################################################################################################################
# # ARGUMENTS PART

# Verifying of the arguments
if [ "$#" -ne 12 ]; then
  echo ""
  echo "Name : $0"
  echo ""
  echo "Usage: "
  echo "	$0 {wildcards.sample} <in> <bamdir> <bamext> <popname> <depth> <cohort> <OutputDir> <project_dir> <GenomeIndexDir> <genomefile> <script_dir> <TmpDir>"
  echo ""
  echo "Description:"
  echo ""
  echo "	This script calculate the negative coverage (RC) for all sites within \$cohort-insertion-sites.[region].bed in the whole-genome alignment bam file (\$bamdir/\${in}\$bamext.bam) and stores the minimum over each interval processed by Process_BAMrc_splitreader.pl for each \$popname in output file: \$popname.\$cohort-insertion.[region].rc"
  echo "	If the read coverage was already calculated over a subset of the sites in \$cohort-insertion-sites.[region].bed it does not re-analyze them to reduce computational time."
  echo "	The \$cohort-insertion-sites.[region].bed input files can be generated directly from the combined list of putative insertion sites (all the \$fam.\$subsetname-insertions.filt.DP\$depth.bed from the Filter_insertions_splitreader.pl) and reformated directly as a bed-file (\$cohort-insertion-sites.0.bed) or shifted 100bp upstream (-100 on both start and stop columns in \$cohort-insertion-sites.100up.bed) or 100bp downstream (\$cohort-insertion-sites.100down.bed)."
  echo ""
  echo "Arguments:"
  echo ""
  echo "	<in>:        		Basename of the BAM file"
  echo "	<bamdir>:   		Directory containing the BAM file"
  echo "	<bamext>:   		BAM file extension (e.g., '_dupl_fixed_paired.bam')"
  echo "	<popname>:          	Genome name (when different from bam file name)"
  echo "	<depth>:    		Number of reads required to call an insertion on 1st pass"
  echo "	<cohort>:           	Name of Cohort"
  echo "	<OutputDir>:     	path/to/output"
  echo "	<project_dir>:     	path/to/files"
  echo "	<GenomeIndexDir>:    	path/to/genomeReference"
  echo "	<genomefile>:     	Name of reference genome file"
  echo "	<script_dir>:     	path/to/scripts"
  echo "	<TmpDir>:			TmpDir"
  echo ""
  exit 1
fi


in=$1 # basename of the bam file 
bamdir=$2 # path to bam file
bamext=$3 # bam file extension 
popname=$4 # genome name (when different from bam file name)
depth=$5 # number of reads required to call an insertion on 1st pass
cohort=$6 # name of cohort
OutputDir=$7 # path/to/output
project_dir=$8 # path/to/files
GenomeIndexDir=$9 # path/to/genome/reference
genomefile=${10} # name of reference genome file
script_dir=${11} # path/to/scripts
tmp_dir=${12}

updownbool=1 ## calculates coverage 100bp up and down of insertion sites

# PID of this batch
IDPID=$$

# Temporary directory
mkdir -p $tmp_dir/$IDPID                                                  
TmpDir=$tmp_dir/$IDPID 

#############################


# Verifying of the arguments

if [ ! -f "$bamdir/$in$bamext.bam" ]; then
  echo ""
  echo "Error: BAM file not found!"
  echo ""
  exit 1
fi

if ! [[ "$depth" =~ ^[0-9]+$ ]] || [ "$depth" -lt 1 ]; then
  echo ""
  echo "Error: Depth must be a positive integer greater than or equal to 1."
  echo ""
  exit 1
fi

if [ -z "$cohort" ]; then
  echo ""
  echo "Error: Cohort name not specified!"
  echo ""
  exit 1
fi

if [ ! -d "$OutputDir" ]; then
  echo ""
  echo "Error: Output directory not found!"
  echo ""
  exit 1
fi

if [ ! -d "$project_dir" ]; then
  echo ""
  echo "Error: Project directory not found!"
  echo ""
  exit 1
fi

if [ ! -d "$GenomeIndexDir" ]; then
  echo ""
  echo "Error: Genome index directory not found!"
  echo ""
  exit 1
fi

if [ ! -f "$GenomeIndexDir/$genomefile" ]; then
  echo ""
  echo "Error: Genome file not found!"
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
# PID of this batch
IDPID=$$

# Temporary directory
#TmpDir=/localtmp/PB-$IDPID
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

#######################################################################################################################################################################################
	# # Calculation of local BAMrc

echo "["$(date +"%y-%m-%d %T")"] Calculating negative coverage over insertions " | tee -a $OutputDir/log.txt  

echo $bamdir/${in}$bamext.bam

cp $project_dir/$cohort-insertion-sites.0.bed  $OutputDir/

if [ -e $OutputDir/$popname.$cohort-insertion.0.rc ]
then
    # #remove all intervals already calculated
    bedtools subtract -a $OutputDir/$cohort-insertion-sites.0.bed -b $OutputDir/$popname.$cohort-insertion.0.rc -f 1 -r > $OutputDir/$cohort-insertion-sites.0.v2.bed 
    mv $OutputDir/$popname.$cohort-insertion.0.rc $OutputDir/$popname.$cohort-insertion.0.v1.rc 
    
    # # calculate coverage over insertions
    bam-readcount -w 1 -f $GenomeIndexDir/$genomefile $bamdir/${in}$bamext.bam -l $OutputDir/$cohort-insertion-sites.0.v2.bed > $OutputDir/$popname.$cohort-insertion.tmp.0.v2.rc
    # # process rc counts to keep only min DP over each insertion
    perl $script_dir/Process_BAMrc_splitreader.pl $cohort $depth $popname $OutputDir 0.v2
    rm $OutputDir/$popname.$cohort-insertion.tmp.0.v2.rc
    rm $OutputDir/$cohort-insertion-sites.0.v2.bed
    
    # # merge read counts with already calculated
    cat $OutputDir/$popname.$cohort-insertion.0.v*.rc > $OutputDir/$popname.$cohort-insertion.0.rc
    rm $OutputDir/$popname.$cohort-insertion.0.v*
else

  # # calculate coverage over insertions
  bam-readcount -w 1 -f $GenomeIndexDir/$genomefile $bamdir/${in}$bamext.bam -l $OutputDir/$cohort-insertion-sites.0.bed > $OutputDir/$popname.$cohort-insertion.tmp.0.rc
  # # process rc counts to keep only min DP over each insertion
  perl $script_dir/Process_BAMrc_splitreader.pl $cohort $depth $popname $OutputDir 0
  rm $OutputDir/$popname.$cohort-insertion.tmp.0.rc

fi
echo "["$(date +"%y-%m-%d %T")"] Finished local BAMrc "  | tee -a $OutputDir/log.txt  


#######################################################################################################################################################################################
	# # Calculation of 100up and 100down BAMrc

if [ $updownbool -gt 0 ]
then
  ## calculates coverage 100bp up and down
  cp $project_dir/$cohort-insertion-sites.100up.bed  $OutputDir/
  if [ -e $OutputDir/$popname.$cohort-insertion.100up.rc ]
  then
      bedtools subtract -a $OutputDir/$cohort-insertion-sites.100up.bed -b $OutputDir/$popname.$cohort-insertion.100up.rc -f 1 -r > $OutputDir/$cohort-insertion-sites.100up.v2.bed ##remove all intervals already calculated
      mv $OutputDir/$popname.$cohort-insertion.100up.rc $OutputDir/$popname.$cohort-insertion.100up.v1.rc 
      bam-readcount -w 1 -f $GenomeIndexDir/$genomefile $bamdir/${in}$bamext.bam -l $OutputDir/$cohort-insertion-sites.100up.v2.bed > $OutputDir/$popname.$cohort-insertion.tmp.100up.v2.rc
      perl $script_dir/Process_BAMrc_splitreader.pl $cohort $depth $popname $OutputDir 100up.v2
      rm $OutputDir/$popname.$cohort-insertion.tmp.100up.v2.rc
      rm $OutputDir/$cohort-insertion-sites.100up.v2.bed
      cat $OutputDir/$popname.$cohort-insertion.100up.v*.rc > $OutputDir/$popname.$cohort-insertion.100up.rc
      rm $OutputDir/$popname.$cohort-insertion.100up.v*
  else
    bam-readcount -w 1 -f $GenomeIndexDir/$genomefile $bamdir/${in}$bamext.bam -l $OutputDir/$cohort-insertion-sites.100up.bed > $OutputDir/$popname.$cohort-insertion.tmp.100up.rc
    # # # process rc counts to keep only min DP over each insertion
    perl $script_dir/Process_BAMrc_splitreader.pl $cohort $depth $popname $OutputDir 100up
    rm $OutputDir/$popname.$cohort-insertion.tmp.100up.rc
  fi
  echo "["$(date +"%y-%m-%d %T")"] Finished 100up BAMrc "  | tee -a $OutputDir/log.txt  

  cp $project_dir/$cohort-insertion-sites.100down.bed  $OutputDir/
  if [ -e $OutputDir/$popname.$cohort-insertion.100down.rc ]
  then
      bedtools subtract -a $OutputDir/$cohort-insertion-sites.100down.bed -b $OutputDir/$popname.$cohort-insertion.100down.rc -f 1 -r > $OutputDir/$cohort-insertion-sites.100down.v2.bed ##remove all intervals already calculated
      mv $OutputDir/$popname.$cohort-insertion.100down.rc $OutputDir/$popname.$cohort-insertion.100down.v1.rc 
      bam-readcount -w 1 -f $GenomeIndexDir/$genomefile $bamdir/${in}$bamext.bam -l $OutputDir/$cohort-insertion-sites.100down.v2.bed > $OutputDir/$popname.$cohort-insertion.tmp.100down.v2.rc
      perl $script_dir/Process_BAMrc_splitreader.pl $cohort $depth $popname $OutputDir 100down.v2
      rm $OutputDir/$popname.$cohort-insertion.tmp.100down.v2.rc
      rm $OutputDir/$cohort-insertion-sites.100down.v2.bed
      cat $OutputDir/$popname.$cohort-insertion.100down.v*.rc > $OutputDir/$popname.$cohort-insertion.100down.rc
      rm $OutputDir/$popname.$cohort-insertion.100down.v*
  else
    bam-readcount -w 1 -f $GenomeIndexDir/$genomefile $bamdir/${in}$bamext.bam -l $OutputDir/$cohort-insertion-sites.100down.bed > $OutputDir/$popname.$cohort-insertion.tmp.100down.rc
    # # process rc counts to keep only min DP over each insertion
    perl $script_dir/Process_BAMrc_splitreader.pl $cohort $depth $popname $OutputDir 100down
    rm $OutputDir/$popname.$cohort-insertion.tmp.100down.rc
  fi
  echo "["$(date +"%y-%m-%d %T")"] Finished 100down BAMrc "  | tee -a $OutputDir/log.txt  
fi


echo "["$(date +"%y-%m-%d %T")"] Finished BAMrc " | tee -a $OutputDir/log.txt 
