#!/bin/bash

FILE=$1
LINEAGE=$2
BASE_ROUTE=NEED_CUSTOMIZATION		#FILL WITH YOUR BASE ROUTE OF YOUR COMPUTER, i.e /home/computer_name/

clear 

echo "

Description: 

Script to split ONE multifasta file, into multiple singlefasta files. You need to provide the file you want to split and an identifier (LINEAGE) to be able to collapse all single fasta files from a single multifasta file into one folder. The Script ask you where you would want to store the output files, but it will be in the computer Hard Disk if you do not provide a different one. Use under your responsability.

@GCabot2021

"

mkdir $BASE_ROUTE/Documentos/
mkdir $BASE_ROUTE/Documentos/Multifasta_to_multiplefastas/
mkdir $BASE_ROUTE/Documentos/Multifasta_to_multiplefastas/$LINEAGE/

while read line
do

	if [[ $line = ">"* ]]; then

		SNAME=`echo $line | sed 's#>##g' | cut -d ' ' -f1`

		echo $SNAME

		echo $line > $BASE_ROUTE/Documentos/Multifasta_to_multiplefastas/$LINEAGE/$SNAME.fasta

	else

		echo $line >> $BASE_ROUTE/Documentos/Multifasta_to_multiplefastas/$LINEAGE/$SNAME.fasta

	fi

done < $FILE

#clear

#echo " Please CHECK your results are properly treated prior to work with them. You ARE the ultimate responsible for your analyses... @GCabot2021 "
