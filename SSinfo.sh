#!/bin/bash

inpdb=$1
outfile=$2

if [ ! -f "$inpdb" ]
then
	exit
fi

#stride $inpdb | grep ASG | awk '{print $3}' | xargs| sed 's/ //g' > $outfile
#stride $inpdb | grep ASG | awk '{print $6}' | xargs| sed 's/ //g' >> $outfil

stride $inpdb > ${outfile}_1

while read line
do
	lb=$(echo $line | awk '{print $1}')
	if [[ $lb == "ASG" ]]
	then
		echo $line >> ${outfile}_2
	fi
done < ${outfile}_1

awk '{print $3}' ${outfile}_2 | xargs | sed 's/ //g' > $outfile
awk '{print $6}' ${outfile}_2 | xargs | sed 's/ //g' >> $outfile

rm ${outfile}_1
rm ${outfile}_2
