#!/bin/bash
file_path=${1}
len_threshold=${2}
output=$(sed "s/.fasta/_len_me_$len_threshold.fasta/" <<< $file_path)

awk -v n=$len_threshold '/^>/{ if(l>n) print b; b=$0;l=0;next }
           {l+=length;b=b ORS $0}END{if(l>n) print b }' $file_path > $output

echo Filtering results saved in $output
echo $(cat $file_path | grep '>' | wc -l) total contigs found
echo $(cat $output | grep '>' | wc -l) contigs are longer than $len_threshold

