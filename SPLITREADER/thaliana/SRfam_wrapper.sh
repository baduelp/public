#!/bin/bash 

project_dir=$1
fam=$2 
popname=$3 
currdir=$4

# make local temporary directory
TmpDir=/localtmp/PB-$IDPID
mkdir -p $TmpDir

function cleanup {
  rm -r -f $TmpDir
  rm -f $TmpDir/*
  }
 trap cleanup EXIT

# # by TE
cd $currdir
end=`wc -l $project_dir/TE_sequence/$fam.TElist | awk '{ print $1}'`
for ((l=1; $l<=$end; l=$l+1)); do
  TE=`sed -n "${l}p" $project_dir/TE_sequence/$fam.TElist | awk '{print $1}'`
  # indNB=1
  if [ ! -d ./$fam/$TE ]
  then 
    mkdir ./$fam/$TE
  fi
    rm ./$fam/$TE/$TE.$popname-insertion-sites*.bed
    grep -w $TE ./ALL/$popname-insertion-sites.bed > ./$fam/$TE/$TE.$popname-insertion-sites.bed
    cut -f 1-4 ./$fam/$TE/$TE.$popname-insertion-sites.bed > $TmpDir/$TE.$popname-insertion-sites.cut.bed
    sortBed -i $TmpDir/$TE.$popname-insertion-sites.cut.bed > $TmpDir/$TE.$popname-insertion-sites.sort1.bed
    mergeBed -i $TmpDir/$TE.$popname-insertion-sites.sort1.bed > $TmpDir/$TE.$popname-insertion-sites.mrg.bed
    sortBed -i $TmpDir/$TE.$popname-insertion-sites.mrg.bed > ./$fam/$TE/$TE.$popname-insertion-sites.sort.bed 
done
echo 'SRfam completed' > ./$fam/$fam.$popname-completed.txt