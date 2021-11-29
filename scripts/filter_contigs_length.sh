#!/bin/bash
file_path=${1}
len_threshold=${2}
output_file=${3}

awk -v n=$len_threshold '/^>/{ if(l>n) print b; b=$0;l=0;next }
            {l+=length;b=b ORS $0}END{if(l>n) print b }' $file_path > $output_file

echo $(cat $file_path | grep '>' | wc -l) total contigs found
echo $(cat $output_file | grep '>' | wc -l) contigs are longer than $len_threshold

