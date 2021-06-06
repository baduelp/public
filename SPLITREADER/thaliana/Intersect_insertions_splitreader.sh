#!/bin/bash

in_dir=$1
out_dir=$2
cohortname=$3
famName=$4
shift 4
fqName=( "$@" )

indNB=${#fqName[@]}
echo -e "$cohortname\t$indNB"
anyinsertion=0

cd $out_dir/

if [ -e interbed.intersect.$famName.$cohortname.e ]
then
	rm -f *.e
	rm *.o*
	rm *.sh

fi

echo '#!/bin/sh' > tampon.intersect.$famName.$cohortname.sh 
echo -n "multiIntersectBed -i " >> tampon.intersect.$famName.$cohortname.sh #-header 

for (( fosmid=0; fosmid<$indNB; fosmid++ ))  
do

	filename=${fqName[$fosmid]}
	if [ -s $in_dir/$famName.$filename-insertion-sites.sort.bed ] #check .bed file exists
	then
		anyinsertion=1
		echo -n "$in_dir/$famName.$filename-insertion-sites.sort.bed " >> $out_dir/tampon.intersect.$famName.$cohortname.sh
	fi
done

cd $out_dir
if [ $anyinsertion -gt 0 ]
then
	echo -n "-names " >> tampon.intersect.$famName.$cohortname.sh
	for (( fosmid=0; fosmid<$indNB; fosmid++ ))  
	do
		filename=${fqName[$fosmid]}
		if [[ -s $in_dir/$famName.$filename-insertion-sites.sort.bed ]] #check .bed file exists
		then
			echo -n "$filename " >> tampon.intersect.$famName.$cohortname.sh
		fi
		
	done

	echo " > $out_dir/$famName.$cohortname-intersect.bed " >> tampon.intersect.$famName.$cohortname.sh

	command="condor_qsub -p 0 -l mem=5G -S /bin/bash -e interbed.intersect.$famName.$cohortname.e -N interbed.intersect.$famName.$cohortname tampon.intersect.$famName.$cohortname.sh"
	latest_id=$($command | awk ' { print $3 }')
	
	echo "Intersect bed $famName.$cohortname submitted: $latest_id"
else
	echo -en "no putative insertion detected\n"
	echo '' > $out_dir/$famName.$cohortname-intersect.bed

fi

