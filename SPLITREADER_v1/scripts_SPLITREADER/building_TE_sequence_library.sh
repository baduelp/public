#!/bin/bash


#######################################################################################################################################################################################
# # ARGUMENTS PART

# Verifying of the arguments

if [ "$#" -ne 6 ]; then
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
TE_lib=$2
TE_seq=$3
TE_list=$4
superfamily_TSD=$5
TEfamily_superfamily=$6

#############################

# Verifying of the arguments

if [ ! -d "$workspace_dir" ]; then
  echo ""
  echo "Error: Working directory not found!"
  echo ""
  exit 1
fi

if [ ! -f "$workspace_dir/TE_sequence/$TE_lib.fasta" ]; then
  echo ""
  echo "Error: TE library file not found!"
  echo ""
  exit 1
fi

#######################################################################################################################################################################################

cd $workspace_dir

cd TE_sequence/

echo -e -n "building TE library... \t"

if [ ! -e $TE_lib.1.bt2 ]; then
    if [ ! -e $TE_lib.fasta ]; then
        if [ -e $workspace_dir/TE_sequence/$TE_lib.fa ]; then
            mv $workspace_dir/TE_sequence/$TE_lib.fa $workspace_dir/TE_sequence/$TE_lib.fasta
        fi
        cp $workspace_dir/TE_sequence/$TE_lib.fasta ./
    fi
    bowtie2-build $TE_lib.fasta $TE_lib
fi

echo -e "done"
cd $workspace_dir/TE_sequence
echo -e -n "building TE list... \t"

if [ ! -s TE-information-all.txt ]; then
    if [ ! -e superfamily_TSD.txt ]; then
        cp $workspace_dir/TE_sequence/superfamily_TSD.txt ./
    fi
    if [ ! -e $TE_lib.fasta ]; then
        cp $workspace_dir/TE_sequence/$TE_lib.fasta ./
    fi

    echo -n '' > TE-information-all.txt
    TElistarr=($(awk '/>/{print}' $TE_lib.fasta | cut -c2- ))
    TENB=`wc -l TE-list.txt | awk '{print $1}'`

    for ((t=1; $t<=$TENB; t++)); do
        TEname=$(sed -n "${t}p" TE-list.txt | awk '{print $1}' | sed -e 's/[\r\n]//g' )
        TEinfo=`grep -w $TEname TE-information-all.txt | wc -l | awk '{print $1}' `
        if [ $TEinfo -eq 0 ]; then
            superfam=$(grep -w $TEname $workspace_dir/TE_sequence/TEfamily-superfamily.txt | cut -f 2 | cut -d"/" -f2)
            TSDlength=`grep -w $superfam superfamily_TSD.txt | awk '{print $2}' `
            TSDlength=${TSDlength//$'\r'}
            echo -e "$TEname\t$TSDlength\t$superfam"
            echo -e "$TEname\t$TSDlength\t$superfam" >> TE-information-all.txt
        fi
    done
fi

echo -e "done"

