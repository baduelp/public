#!/bin/bash 

#############################################
#                                           #  
#               SPLITREADER                 #
#              SRfam_wrapper                #
#           Baduel et al. 2020              #
#                                           # 
#############################################

##Questions or comments to pbaduel(ar)bio.ens.psl.eu

#######################################################################################################################################################################################
# # ARGUMENTS PART

# Verifying of the arguments
if [ "$#" -ne 5 ]; then
  echo ""
  echo "Name : $0"
  echo ""
  echo "Usage: "
  echo "	$0 <project_dir> <fam> <popname> <currdir> <TmpDir>"
  echo ""
  echo "Description:"
  echo ""
  echo "	Splits all the insertions found with SPLITREADER-beta1.5_part2.sh into families and superfamilies"
  echo ""
  echo "Arguments:"
  echo ""
  echo "	<project_dir>:  Workspace Directory"
  echo "	<fam>:   	TE superfamily"
  echo "	<popname>:   	Sample name"
  echo "	<currdir>:      \$project_dir/BEDfiles/SPLITREADER/\$cohortname"
  echo "	<TmpDir>:       TmpDir"
  echo ""
  exit 1
fi

project_dir=$1
fam=$2 
popname=$3 
currdir=$4
tmp_dir=$5

# PID of this batch
IDPID=$$

# Temporary directory
mkdir -p $tmp_dir/$IDPID                                                  
TmpDir=$tmp_dir/$IDPID 

#############################


# Verifying of the arguments

if [ ! -d "$project_dir" ]; then
  echo ""
  echo "Error: Workspace directory not found!"
  echo ""
  exit 1
fi

if [ ! -d "$TmpDir" ]; then
  echo ""
  echo "Error: Workspace directory not found!"
  echo ""
  exit 1
fi

if [ ! -f "$project_dir/TE_sequence/$fam.TElist" ]; then
  echo ""
  echo "Error: TE list file for family $fam not found!"
  echo ""
  exit 1
fi

if [ -z "$popname" ]; then
  echo ""
  echo "Error: Sample name not specified!"
  echo ""
  exit 1
fi

if [ ! -d "$currdir" ]; then
  echo ""
  echo "Error: \$currdir not found! (should be: \$project_dir/BEDfiles/SPLITREADER/\$cohortname)"
  echo ""
  exit 1
fi

#######################################################################################################################################################################################

function cleanup {
  rm -r -f $TmpDir
  rm -f $TmpDir/*
  }
 trap cleanup EXIT

#######################################################################################################################################################################################
	# # Seek TE names for each superfamily in $project_dir/TE_sequence/$fam.TElist, then splits the insertions found by SPLITREADER-beta1.5_part2.sh and sorted by SPLITREADER-sort.sh

cd $currdir

if [ ! -d ./$fam ]
  then 
    mkdir $currdir/$fam
  fi

end=`wc -l $project_dir/TE_sequence/$fam.TElist | awk '{ print $1}'`
for ((l=1; $l<=$end; l=$l+1)); do
	TE=`sed -n "${l}p" $project_dir/TE_sequence/$fam.TElist | awk '{print $1}'`

  	if [ ! -d ./$fam/$TE ]
  	then 
    	mkdir $currdir/$fam/$TE
  	fi

    grep -w $TE $currdir/ALL/$popname-insertion-sites.bed > $currdir/$fam/$TE/$TE.$popname-insertion-sites.bed
    cut -f 1-4 $currdir/$fam/$TE/$TE.$popname-insertion-sites.bed > $TmpDir/$TE.$popname-insertion-sites.cut.bed
    sortBed -i $TmpDir/$TE.$popname-insertion-sites.cut.bed > $TmpDir/$TE.$popname-insertion-sites.sort1.bed
    mergeBed -i $TmpDir/$TE.$popname-insertion-sites.sort1.bed > $TmpDir/$TE.$popname-insertion-sites.mrg.bed
    sortBed -i $TmpDir/$TE.$popname-insertion-sites.mrg.bed > $currdir/$fam/$TE/$TE.$popname-insertion-sites.sort.bed 
done

echo 'SRfam completed' > $currdir/$fam/$fam.$popname-completed.txt

