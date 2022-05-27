#!/bin/bash
fasta=${1}
csplit -z $fasta '/>/' '{*}'
