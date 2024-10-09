#!/bin/bash


#######################################################################################################################################################################################
# # ARGUMENTS PART

# Verifying of the arguments

if [ "$#" -ne 3 ]; then
	  echo ""
  echo "Name : $0"
  echo ""
  echo "Description:"
  echo ""
  echo "	Prerequisites for the SPLITREADER pipeline."
  echo ""
  exit 1
fi

# Arguments :

workspace_dir=$1
ref_dir=$2
genome=$3

#############################

# Verifying of the arguments

if [ ! -d "$workspace_dir" ]; then
  echo ""
  echo "Error: Working directory not found!"
  echo ""
  exit 1
fi

if [ ! -d "$workspace_dir/Reference" ]; then
  echo ""
  echo "Error: Reference directory not found!"
  echo ""
  exit 1
fi

if [ ! -f "$workspace_dir/Reference/$genome.fasta" ]; then
  echo ""
  echo "Error: Genome fasta file not found!"
  echo ""
  exit 1
fi

#######################################################################################################################################################################################

cd $workspace_dir

cd Reference/

echo -e -n "normalizing REF genome... \t"

if [ ! -e $genome.norm.fasta ]; then
    picard NormalizeFasta I=$genome.fasta O=$genome.norm.fasta 
    cp $genome.norm.fasta ./$genome.fasta
fi

echo -e "done"
echo -e -n "building REF genome... \t"

if [ ! -e $genome.1.bt2 ]; then
    bowtie2-build $genome.fasta $genome
fi

echo -e "done"
echo -e -n "indexing REF genome... \t"

if [ ! -e $genome.fasta.fai ]; then
    samtools faidx $genome.fasta
fi

echo -e "done"

