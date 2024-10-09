#!/bin/bash

workspace_dir=$1
cohortname=$2
bamext=$3
bam_dir=$4

famNames=$(awk 'NR>1{{print $1}}' $workspace_dir/TE_sequence/superfamily_TSD.txt)
SampleNames=($(ls "$bam_dir"/*.bam | sed -e "s/${bamext}.bam$//" | xargs -n1 basename)) 
cd $workspace_dir  

all_files_exist=true

for fam in $famNames 
do
    currdir=$workspace_dir/BEDfiles/SPLITREADER/$cohortname/$fam 
    TEs=$(awk '{{print $1}}' $workspace_dir/TE_sequence/$fam.TElist) 
    
    for TE in $TEs 
    do 
        cd $currdir/$TE 

		
        bed_files=()
        for file in $TE.*-insertion-sites.sort.bed; do
            bed_files+=("$file")
        done
				#bed_files=($(ls *.sort.bed)) # Start of command to execute multiIntersectBed 
				#echo "${{bed_files[@]}}"
        cmd="multiIntersectBed -i " # Add each .bed file to command 
        anyinsertion=0

        for bed_file in "${bed_files[@]}" 
        do 
            if [[ -s $bed_file ]]; then
                anyinsertion=1
                cmd+="$bed_file "
            fi 
        done 
        
        if [ $anyinsertion -gt 0 ]; then
            cmd+="-names " 
            for sample in "${SampleNames[@]}" 
            do 
                if [[ -s $TE.$sample-insertion-sites.bed ]]; then                    
                    cmd+="$sample " 
                fi                
            done 

            $cmd> $TE.$cohortname-intersect.bed 2> $currdir/$TE/$TE.intersect.err
        else
            echo '' > $TE.$cohortname-intersect.bed
        fi

        if [ ! -f $TE.$cohortname-intersect.bed ]; then 
            all_files_exist=false 
        fi
    done 

    if $all_files_exist ; then 
        touch $currdir/fam_intersect_completed.txt 
    fi 
done

