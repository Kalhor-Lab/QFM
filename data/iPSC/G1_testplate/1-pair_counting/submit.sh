#!/bin/bash


#sys=${1:-'osx'}

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
	sys="linux";
else
	sys='osx';

fi

for i in `ls ../0-raw_data_PB/*R1_001.fastq`; do
	prename="$(cut -d'/' -f3 <<<$i)"; name="$(cut -d'_' -f1 <<<$prename)"
        date=`date`
        echo processing $name,  $date.
	perl -w Count_PB.pl $i $sys
	# echo "${?}"
	if [ $? != 0 ]
	then
		echo "Failed to process $i"
		break
	fi
done
