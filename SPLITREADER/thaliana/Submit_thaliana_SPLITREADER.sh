#!/bin/bash

echo -e "\n\n"
echo -e "##############################"
echo -e "       SPLITREADER V2.5"
echo -e "   Baduel et al. MMB 2021"
echo -e "       implemented by"
echo -e "        Pierre Baduel"
echo -e "##############################\n"

# # # List of variable names
workspace_dir='/path/to/working/directory'
ref_dir="$workspace_dir/Reference/" # path to reference genome annotation files
genome='TAIR10'
TE_lib='TE_all_Athaliana' 
TE_annotation='TAIR10_Quesneville_GFF3_transposons_only'
bam_dir="$workspace_dir/BAMs/"
bamext='_dupl_fixed_paired' # extension of bam files
script_dir='/path/to/scripts/' # location of SPLITREADER scripts

# # # Required files
# in $ref_dir a $genome.fasta file with the reference genome sequence
# in $workspace_dir/TE_sequence a tab-delimited superfamily_TSD.txt file with superfamilies in the 1st column and their respective TSD lengths in the 2nd column 
# in $workspace_dir/TE_sequence a TE-list.txt file with a list of the names of the TEs annotated in the reference genome
# in $workspace_dir/TE_sequence a $TE_lib.fa file with the sequence of all the TEs annotated in the reference genome
# in $workspace_dir/TE_sequence a tab-delimited TEfamily-superfamily.txt file with TE names in the 1st column and their respective superfamily in the 2nd column

depth=3 # Minimum number of reads (split+discordant)

# # # Shortcuts
rerunSR1=0
rerunSR2=0
rerunSRsort=0
run_part=10
rerunSRfam=0
refilt1=0
refilt2=0
refilt3=0

# # # Uncomment lines to launch corresponding steps
run_part=1 # go to run SPLITREADER PARTS 1 & 2
# # # rerunSR1=1 # force rerun of SPLITREADER part 1
# # # rerunSR2=1 # force rerun of SPLITREADER part 2
# # # rerunSRsort=1 # force resort of SPLITREADER output
# run_part=2 # go to process SPLITREADER output
# # rerunSRfam=1 # force rerun of processing SPLITREADER output by TE family
# # refilt1=1 # force rerun of filtering SPLITREADER output based on positive coverage (split+disc)
# # refilt2=1 # force restart of counting negative coverage over putative insertion sites (BAM read count)
# # refilt3 # force restart of filtering based on negative coverage drop


# # # List of samples
    cd $bam_dir/
    SampleNames=($(ls *.bam | sed -e 's/\_dupl_fixed_paired.bam$//')) # copy there the bam files extension ($bamext)
    indNB=${#SampleNames[@]}
    cohortname="TE_capture_accessions"

    cd $workspace_dir
    if [ ! -d $cohortname ]
    then
        mkdir $cohortname
        mkdir BAMs
        mkdir BEDfiles
        cd BEDfiles
        if [ ! -d SPLITREADER ]
        then
            mkdir SPLITREADER
        fi
    fi
    if [ ! -e $workspace_dir/BAMs/$cohortname.fastq_names.txt ]; then
        # create list of filenames for BAMs (column 1) and fastq (column 2) in case they differ
        echo '' > $workspace_dir/BAMs/$cohortname.fastq_names.txt
        echo '' > $workspace_dir/BAMs/$cohortname.missing_bams.txt   
                
    fi

cd $workspace_dir

echo -e "------  Project $cohortname  ------\n"
echo -e "### Checking requirements  ###\n"

# # # Building reference genome bowtie2 index
    cd $workspace_dir
    if [ ! -d Reference ]
    then
        mkdir Reference
    fi
    cd Reference/

    if [ ! -e $genome.fasta ]
    then
        cp $ref_dir/$genome.fasta ./
    fi
    echo -e -n "normalizing REF genome... \t"
    
    if [ ! -e $genome.norm.fasta ]
    then
        # java -jar /import/bc_users/bioinfo/vincens/pkg/te-tracker/picard-tools-1.119/NormalizeFasta.jar I=$genome.fasta O=$genome.norm.fasta 
        picard NormalizeFasta I=$genome.fasta O=$genome.norm.fasta 
        cp $genome.norm.fasta ./$genome.fasta
    fi
    echo -e "done"

    echo -e -n "building REF genome... \t"
    if [ ! -e $genome.1.bt2 ]
    then
        # gunzip $genome.fasta.gz
        bowtie2-build $genome.fasta $genome
    fi
    echo -e "done"

    echo -e -n "indexing REF genome... \t"
    if [ ! -e $genome.fasta.fai ]
    then
        samtools faidx $genome.fasta
    fi
    echo -e "done"

# # # Building TE sequence library
    cd $workspace_dir
    if [ ! -d TE_sequence ]
    then
        mkdir TE_sequence
    fi
    cd TE_sequence/

    echo -e -n "building TE library... \t"
    if [ ! -e $TE_lib.1.bt2 ]
    then
        if [ ! -e $TE_lib.fasta ]
        then
            if [ -e $workspace_dir/TE_sequence/$TE_lib.fa ]; then
                mv $workspace_dir/TE_sequence/$TE_lib.fa $workspace_dir/TE_sequence/$TE_lib.fasta
            fi
            cp $workspace_dir/TE_sequence/$TE_lib.fasta ./
        fi
        bowtie2-build $TE_lib.fasta $TE_lib
    fi
    echo -e "done"

# # # Building TE sequence library
    cd $workspace_dir/TE_sequence
    echo -e -n "building TE list... \t"

    if [ ! -s TE-information-all.txt ]
    then
        if [ ! -e superfamily_TSD.txt ]
        then
            cp $workspace_dir/TE_sequence/superfamily_TSD.txt ./
        fi

        if [ ! -e $TE_lib.fasta ]
        then
            cp $workspace_dir/TE_sequence/$TE_lib.fasta ./
        fi

        echo -n '' > TE-information-all.txt
        TElistarr=($(awk '/>/{print}' $TE_lib.fasta | cut -c2- )) #
        TENB=`wc -l TE-list.txt | awk '{print $1}'`
        # # TENB=1
        for ((t=1; $t<=$TENB; t++)); do
            TEname=$(sed -n "${t}p" TE-list.txt | awk '{print $1}' | sed -e 's/[\r\n]//g' )
            TEinfo=`grep -w $TEname TE-information-all.txt | wc -l | awk '{print $1}' `
            # echo -e $TEname "\t" $TEinfo
            if [ $TEinfo -eq 0 ]; then
                superfam=$(grep -w $TEname $workspace_dir/TE_sequence/TEfamily-superfamily.txt | cut -f 2 | cut -d"/" -f2)
                TSDlength=`grep -w $superfam superfamily_TSD.txt | awk '{print $2}' `
                TSDlength=${TSDlength//$'\r'}
                echo -e "$TEname\t$TSDlength\t$superfam" 
                echo -e "$TEname\t$TSDlength\t$superfam" >> TE-information-all.txt # assemble list of TEs for SR2 | tr -d '\n' 
            fi
        done
    fi
    echo -e "done"

# # # Building family TE list
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
            echo ''
        fi
    done
    echo -e "done"

# # # Running SPLITREADER pipeline

if [ $run_part -eq 1 ]
then
    echo -e "###    Running SPLITREADER pipeline    ###\n"

    echo "executable = $script_dir/SPLITREADER-beta1.5_part1.sh" > $workspace_dir/tmp.SR1.sub
    echo "executable = $script_dir/SPLITREADER-beta1.5_part2.sh" > $workspace_dir/tmp.SR2.sub
    echo "executable = $script_dir/SPLITREADER-sort.sh" > $workspace_dir/tmp.SRsort.sub
    
    if [ ! -s $workspace_dir/BAMs/$cohortname.cov ]
    then
        # create list to store genome wide coverages
        echo -e "sample\tmedcov\tmeancov" > $workspace_dir/BAMs/$cohortname.cov
    fi

    for (( sample=0; sample<$indNB; sample++ ))  
    do

        cd $workspace_dir
        filename=${SampleNames[$sample]}
        missing_bam=1
        echo -e "\n#### $filename ####"
        
        cd $bam_dir
        if [ -e $bam_dir/${filename}$bamext.bam ]
        then
            echo "BAM found"
            echo -e "$filename\t$filename\n" >> $workspace_dir/BAMs/$cohortname.fastq_names.txt
        else
            echo "missing bam >> SKIP"
            echo $filename >> $workspace_dir/BAMs/$cohortname.missing_bams.txt
            missing_bam=0
        fi
            
        if [ $missing_bam -gt 0 ]  
        then

            if [ ! -d "/$workspace_dir/$cohortname/$filename" ]
            then
                mkdir /$workspace_dir/$cohortname/$filename
            fi

            

            cd /$workspace_dir/$cohortname/$filename/
            if [ ! -d "part1" ]
            then
                mkdir ./part1
            fi
            cd /$workspace_dir/$cohortname/$filename/part1
            
            curr_dir=$(pwd)

            if [ $rerunSR1 -gt 0 ]; then
                rm log.txt
            fi

            if [ -e log.txt ]
            then
                job_status=$(grep "Finished running" log.txt | wc -l)
            else 
                job_status=0
            fi
            
            if [ -e ${filename}-TE.bam ]
            then
                sam_status=$(ls -l ${filename}-TE.bam | awk '{print $5}')
                if [ "$sam_status" -lt 500000 ]
                then
                    echo "misformed TE-bam: $filename"
                    sam_status=0
                fi
            else 
                sam_status=0
            fi
            
            if [[ ("$job_status" -eq 0) || ("$sam_status" -lt 10 ) ]]
            then
            
                if [ -e $filename.SR1.e ]
                then
                    rm *.e
                    rm *.bam
                fi

                if [ -d "/$workspace_dir/$cohortname/$filename/part2" ]
                then
                    if [ -e "/$workspace_dir/$cohortname/$filename/part2/log.txt" ]
                    then
                        rm /$workspace_dir/$cohortname/$filename/part2/log.txt
                    fi
                fi

                echo "SR1 submitted..." 
                $script_dir/SPLITREADER-beta1.5_part1.sh $filename part1 $bam_dir $bamext $cohortname $workspace_dir $TE_lib > $filename.SR1.e
                
            elif [ "$sam_status" -gt 100000 ]
            then
                echo "SR1 completed"
                cd /$workspace_dir/$cohortname/$filename/
                if [ ! -d "part2" ]
                then
                    mkdir ./part2
                fi
                cd /$workspace_dir/$cohortname/$filename/part2
                curr_dir=$(pwd)
                
                if [ $rerunSR2 -gt 0 ]; then
                    rm log.txt
                    rm $filename-insertion-sites.bed
                fi

                job_status=0
                if [[ ( -e log.txt ) && ( -e $filename-insertion-sites.bed ) ]]
                then
                    job_status=$(grep "SPLITREADER PART2 COMPLETED" log.txt | wc -l)
                fi

                if [ "$job_status" -gt 0 ]
                then
                    echo "SR2 completed"

                    sort_job_status=$(grep "Finished sorting SPLITREADER output" log.txt | wc -l)
                    if [ $rerunSRsort -gt 0 ]; then
                        sort_job_status=0
                        sed -i '/Finished sorting SPLITREADER output/d' log.txt
                    fi
                    
                    if [ ! -d $workspace_dir/BEDfiles/SPLITREADER/$cohortname ]; then
                        mkdir $workspace_dir/BEDfiles/SPLITREADER/$cohortname
                    fi


                    if [ "$sort_job_status" -gt 0 ]
                    then
                        echo "SRsort completed"
                        if [ ! -s $workspace_dir/BEDfiles/SPLITREADER/$cohortname/$filename-insertion-sites.sort.bed ]
                        then
                                echo "no insertions detected"
                                echo "" > $workspace_dir/BEDfiles/SPLITREADER/$cohortname/$filename-insertion-sites.sort.bed
                        fi
                    else
                        echo "SRsort submitted..." 
                        $script_dir/SPLITREADER-sort.sh $filename $workspace_dir $workspace_dir/BEDfiles/SPLITREADER/$cohortname $filename $cohortname
                    fi
                else
                    if [ -e $filename.SR2.e ]
                    then
                        rm *.e
                        rm *.bam*
                    fi
                    
                    echo "SR2 submitted..." 
                    $script_dir/SPLITREADER-beta1.5_part2.sh $filename $cohortname $workspace_dir Reference/$genome $TE_annotation > $filename.SR1.e
                    
                    
                fi
            fi

        fi
        echo ''
    done

# # # Process SPLITREADER output
elif [ $run_part -eq 2 ]
then


    cd $workspace_dir/BEDfiles/SPLITREADER/
    curr_dir=$(pwd)

    filtname="filt$depth"

    all_ind_finished=0
    if [ $run_part -le 2 ]
    then
        
        echo -e "Calling $cohortname insertions by family...\n"
        cd $workspace_dir/BEDfiles/SPLITREADER/$cohortname
        rm *.log
        rm *.e

        currdir=$(pwd)
        if [ ! -d $currdir/ALL ]; then
            mkdir $currdir/ALL
        fi

        if [ -e ${SampleNames[0]}-insertion-sites.sort.bed ]
        then
            mv $currdir/*.bed $currdir/ALL/
        fi

        fam_finished=0
        for (( famI=0; famI<$famNB; famI++ ))  
        do
            fam=${famName[$famI]}
            echo -e "\n\t#### $fam ####"
            if [ ! -d $currdir/$fam ]
            then 
                mkdir $currdir/$fam
            fi

            if [ -e $fam.${SampleNames[0]}.SRfam.e ]
            then
                rm $fam.*.e
            fi

            if [ $rerunSRfam -gt 0 ]; then
                rm $currdir/$fam.$cohortname-insertions.$filtname.DP$depth.bed
                rm $currdir/$fam/$fam.*-completed.txt

            fi

            if [ ! -e $currdir/$fam.$cohortname-insertions.$filtname.DP$depth.bed ]
            then
                if [ ! -e $currdir/$fam/$fam.ALL-completed.txt ]
                then
                    echo -e "\textracting family insertions..."
                    
                    echo "executable = $script_dir/SRfam_wrapper.sh" > $workspace_dir/tmp.SRfam.sub
                    # indNB=1
                    indDONE=0
                    for (( sample=0; sample<$indNB; sample++ ))  
                    do
                        filename=${SampleNames[$sample]}
                        echo -e -n "\t## $filename... \t"
                        
                        if [ ! -e $currdir/$fam/$fam.$filename-completed.txt ]
                        then
                            $script_dir/SRfam_wrapper.sh $workspace_dir $fam $filename $currdir > $fam.$filename.SRfam.e
                            echo -e "SRfam... \t submitted"
                        else
                            echo -e "SRfam... \t done"
                            indDONE=$((indDONE+1))
                        fi
                    done
                    
                    if [ $indDONE -eq $indNB ]; then
                        echo 'SRfam completed' > $currdir/$fam/$fam.ALL-completed.txt 
                    fi 
                    echo '' 

                
                else
                    cd $currdir/$fam
                    end=`wc -l $workspace_dir/TE_sequence/$fam.TElist | awk '{print $1}'`
                    all_TE_finished=1
                    for ((l=1; $l<=$end; l=$l+1)); do
                        TE=`sed -n "${l}p" $workspace_dir/TE_sequence/$fam.TElist | awk '{print $1}'`
                        cd $currdir/$fam/$TE
                        echo -e -n "\t$TE...\t"
                        
                        if [ $refilt1 -gt 0 ]; then
                            rm $TE.$cohortname-insertions.$filtname.DP$depth.bed
                            rm $TE.$cohortname-intersect.*bed
                        fi

                        if [ ! -e $TE.$cohortname-intersect.bed ]
                        then
                            all_TE_finished=0
                            echo -e -n "intersecting individual BED files...\t"

                            # # COMPUTE ALL OVERLAPS BETWEEN PUTATIVE INSERTIONS
                            $script_dir/Intersect_insertions_splitreader.sh $currdir/$fam/$TE $currdir/$fam/$TE $cohortname $TE "${SampleNames[@]}"
                                                            
                        elif [ ! -e $TE.$cohortname-insertions.$filtname.DP$depth.bed ]
                        then
                            all_TE_finished=0    
                            if [ ! -e $TE.$cohortname-intersect.sort.bed ]
                            then
                                cut -f 1-4 $TE.$cohortname-intersect.bed > $TE.$cohortname-intersect.cut.bed
                                sortBed -i $TE.$cohortname-intersect.cut.bed > $TE.$cohortname-intersect.sort.bed
                                mergeBed -d 5 -i $TE.$cohortname-intersect.sort.bed > $TE.$cohortname-intersect.sort1.bed ## MERGE CONTIGUOUS OR NEAR-CONTIGUOUS (5bp) INSERTIONS (BY TE!) TO AVOID DUPLICATES
                                sortBed -i $TE.$cohortname-intersect.sort1.bed > $TE.$cohortname-intersect.mrg.bed
                                rm $TE.$cohortname-intersect.sort1.bed
                                rm $TE.$cohortname-intersect.cut.bed
                            fi

                            if [ -s $TE.$cohortname-intersect.mrg.bed ]
                            then
                                # # FILTER PUTATIVE INSERTIONS BASED ON OVERLAPS
                                echo -e -n "\tfilter DP$depth...\t"
                                $script_dir/Filter_insertions_splitreader.pl $cohortname $depth $workspace_dir $TE $currdir/$fam/$TE 10 50 ${SampleNames[*]} > $TE.$cohortname.SRfilt.e
                            else
                                echo -e -n "\tno insertions detected...\tdone"
                                echo '' > $TE.$cohortname-insertions.$filtname.DP$depth.bed
                            fi
                        else 
                        
                            if [ -e $TE.$cohortname-intersect.mrg.bed ]
                            then
                                rm $TE.$cohortname-intersect.mrg.bed
                                rm $TE.$cohortname-intersect.sort.bed
                                rm *.e
                            fi

                            if [ -e $TE.$cohortname.SRfilt.e ]
                            then
                                rm $TE.$cohortname.SRfilt.e
                            fi

                            cp $TE.$cohortname-insertions.$filtname.DP$depth.bed $currdir/$fam/
                            echo -e -n "\tfilter DP$depth COMPLETED!"                
                        fi
                        echo ''
                    done
                    if [ $all_TE_finished -gt 0 ]; then
                        echo -e -n "\n\tall individual $fam filter DP$depth COMPLETED!\n"        
                        
                        cat $currdir/$fam/*.$cohortname-insertions.$filtname.DP$depth.bed > $currdir/$fam.$cohortname-insertions.$filtname.DP$depth.bed

                        # # #
                        header=$(head -n1 $currdir/$fam.$cohortname-insertions.$filtname.DP$depth.bed) #
                        sed -i '/start/d' $currdir/$fam.$cohortname-insertions.$filtname.DP$depth.bed
                        sed -i '/^$/d' $currdir/$fam.$cohortname-insertions.$filtname.DP$depth.bed
                        echo $header | cat - $currdir/$fam.$cohortname-insertions.$filtname.DP$depth.bed > $currdir/temp && mv $currdir/temp $currdir/$fam.$cohortname-insertions.$filtname.DP$depth.bed
                    else
                        echo -e -n "\n\tindividual $fam filter DP$depth MISSING!\n" 
                    fi
                fi
            else
                fam_finished=$(($fam_finished+1))
                echo -e -n "\t$fam filter DP$depth COMPLETED!\n"  

            fi
        done

        
        if [ $fam_finished -eq $famNB ]
        then
            cd $workspace_dir/BEDfiles/SPLITREADER/     
            echo -e "\nALL FAMS SPLITREADER CALLS DP$depth-filter COMPLETED!\n"      

            echo -e "Filtering $cohortname insertions: NC-filter...\n"      
            
            if [[ $refilt2 -gt 0 || ! -e $cohortname-insertions.$filtname.DP$depth.bed ]]; then
                cat ./$cohortname/*.$cohortname-insertions.$filtname.DP$depth.bed > ./$cohortname-insertions.$filtname.DP$depth.bed
            fi
            
            sed -i '/start/ d' ./$cohortname-insertions.$filtname.DP$depth.bed
            # # make list of insertions for bam-readcount
            cut -f 1-3 $cohortname-insertions.$filtname.DP$depth.bed > $cohortname-insertion-sites.0.bed
            sed -i '/start/ d' $cohortname-insertion-sites.0.bed 
            cut -f 1-3 $cohortname-insertions.$filtname.DP$depth.bed | awk -v OFS='\t' '{$2 = $2 - 100; $3 = $3 - 100 ; if ( $3 > 0 ) print $1, $2, $3}' > $cohortname-insertion-sites.100up.bed
            sed -i '/start/ d' $cohortname-insertion-sites.100up.bed 
            cut -f 1-3 $cohortname-insertions.$filtname.DP$depth.bed | awk -v OFS='\t' '{$2 = $2 + 100; $3 = $3 + 100 ; print $1, $2, $3}' > $cohortname-insertion-sites.100down.bed
            sed -i '/start/ d' $cohortname-insertion-sites.100down.bed 
            
            echo 'executable = $script_dir/BAM-readcount_wrapper.sh' > $workspace_dir/tmp.BAMrc.sub

            all_ind_finished=1
            # indNB=1
            for (( sample=0; sample<$indNB; sample++ ))  
            do

                filename=${SampleNames[$sample]}
                
                echo -e "\n\t#### $filename ####"
                missing_bam=1
                cd $bam_dir
                if [ -e $bam_dir/${filename}$bamext.bam ]
                then
                    echo -e "\tBAM found"

                    if [ ! -e $bam_dir/${filename}$bamext.bam.bai ]
                    then
                        if [ ! -e $bam_dir/${filename}$bamext.bai ]
                        then
                            echo -e "\tbuilding BAM index..."
                            samtools index $bam_dir/${filename}$bamext.bam
                        else
                            cp $bam_dir/${filename}$bamext.bai $bam_dir/${filename}$bamext.bam.bai
                        fi
                    fi
                    
                    DP_info=$(grep -w $filename $workspace_dir/BAMs/$cohortname.cov | wc -l)
                    if [ $DP_info -lt 2 ]
                    then
                        bedtools genomecov -ibam $bam_dir/${filename}$bamext.bam | awk -v printbool1=1 -v name="$filename" -v file1="$workspace_dir/BAMs/$cohortname.cov" ' { sumFRAC+=$5; sumDP+=$3*$2; WG_length+=$3; if (sumFRAC >= .5 && printbool1 > 0) { medDP=$2; } } END { sumDP=sumDP/WG_length; print name,"\011",medDP,"\011",sumDP >> file1; }'
                    fi

                else
                    echo -e "\tmissing bam >> SKIP"
                    echo $filename >> $workspace_dir/BAMs/$cohortname.missing_bams.txt
                    missing_bam=0
                fi

                cd $workspace_dir

                if [ $missing_bam -gt 0 ]
                then
                    cd $workspace_dir/$cohortname/$filename/
                    if [ ! -d "BAMrc" ]
                    then
                        mkdir ./BAMrc
                    fi
                    cd $workspace_dir/$cohortname/$filename/BAMrc
                    curr_dir=$(pwd)
                    job_status=0
                    if [[ -e $filename.$cohortname-insertion.100down.rc && -e $filename.$cohortname-insertion.100up.rc && -e $filename.$cohortname-insertion.0.rc ]]
                    then
                        rc_size=$(ls -l $filename.$cohortname-insertion.0.rc | awk '{print $5}')
                        rc_size_up=$(ls -l $filename.$cohortname-insertion.100up.rc | awk '{print $5}')
                        rc_size_down=$(ls -l $filename.$cohortname-insertion.100down.rc | awk '{print $5}')
                        if [[ $rc_size -ge 10000 && $rc_size_up -ge 10000 && $rc_size_down -ge 10000 ]]
                        then 
                            job_status=$(grep "Finished BAMrc" $filename.BAMrc.e | wc -l)     
                        else
                            echo -e "\tmissing readcount"
                        fi     
                    else
                        echo -e "\tmissing readcount"
                    fi
                    
                    if [ $refilt2 -gt 0 ]; then
                        job_status=0 ## force restart of BAM readcount
                        # rm *
                    fi

                    if [ "$job_status" -gt 0 ]
                    then
                        echo -e "\tBAM readcount completed\n"
                    else
                        all_ind_finished=0
                        # # CALCULATE NEGATIVE COVERAGE by chromosome 
                        echo -e "\tBAM readcount submitted...\n"
                        $script_dir/BAM-readcount_wrapper.sh $filename BAMrc $bam_dir $bamext $filename $depth $cohortname $curr_dir $workspace_dir/BEDfiles/SPLITREADER/ $ref_dir/ $genome.fasta > $filename.BAMrc.e
                    fi
                fi
            done
            
        fi
    else
        all_ind_finished=1

    fi
    
    if [ $all_ind_finished -eq 1 ]; then
        echo 'all BAMrc completed!'

        echo -e "\n######"
        cd $workspace_dir/BEDfiles/SPLITREADER/
        curr_dir=$(pwd)

        maskname='nomask' # at this step problematic regions of the genome can be masked 

        if [ $refilt3 -gt 0 ]; then
            rm  $cohortname-insertions.ratioNC$filtname*.DP$depth.$maskname.bed
            rm $cohortname-insertions.$filtname.DP$depth.header.bed
        fi
        
        if [ ! -e $cohortname-insertions.ratioNC$filtname.NConly.DP$depth.$maskname.bed ]
        then
            if [ ! -e $cohortname-insertions.$filtname.DP$depth.header.bed ]
            then
                header=$(head -n1 ./$cohortname/${famName[0]}.$cohortname-insertions.$filtname.DP$depth.bed) 
                echo  $header | tr " " "\t" | cat - $cohortname-insertions.$filtname.DP$depth.bed > $cohortname-insertions.$filtname.DP$depth.header.bed
            fi

            cut -f 1-4 $cohortname-insertions.$filtname.DP$depth.bed > $cohortname-insertion-sites.bed
            sed -i '/start/ d' $cohortname-insertion-sites.bed 
            bedtools sort -i $cohortname-insertion-sites.bed > $cohortname-insertion-sites.sort.bed
            rm $cohortname-insertion-sites.bed
            
            # if nomask simply copy all putative insertions otherwise add step to remove desired regions 
            cp $cohortname-insertion-sites.sort.bed ./$cohortname-insertions.$filtname.DP$depth.$maskname.bed

            echo "Merge and filter NC BED file"
            # # # ADD NEGATIVE COVERAGE
            $script_dir/Filter_negative_calls_splitreader.pl $cohortname $depth $workspace_dir $workspace_dir $filtname $maskname $cohortname ${SampleNames[*]}
            echo "SPLITREADER merging and filtering negative calls DP$depth submitted"
            
         else 
            echo -e "SPLITREADER $cohortname DONE!!\n"
            ALLDONE=1
        fi

    fi
fi

echo -e "\n\n\n\n"