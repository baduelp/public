#!/bin/bash

#######################################################################################################################################################################################
# # ARGUMENTS PART

# Verifying of the arguments

if [ "$#" -ne 1 ]; then
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

#############################

# Verifying of the arguments

if [ ! -d "$workspace_dir" ]; then
  echo ""
  echo "Error: Working directory not found!"
  echo ""
  exit 1
fi


#######################################################################################################################################################################################

cd $workspace_dir/TE_sequence
    echo -e -n "building fam list... \t"
    
    famName=( $(awk 'NR>1{print $1}' $workspace_dir/TE_sequence/superfamily_TSD.txt) )
    famNB=${#famName[@]}

    for (( famI=0; famI<$famNB; famI++ ))  
    do
        fam=${famName[$famI]}
        
        if [ ! -e $fam.TElist ]
        then
        	echo -e -n "\n\t$fam... \t"
        	grep -w $fam TE-information-all.txt | awk '{print $1}' > $fam.TElist
        	echo -e -n "done"
        fi

        if [ $famI -eq $famNB ]; then
            echo '\n'
        fi
    done
    echo -e "done"

