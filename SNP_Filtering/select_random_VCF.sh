#!/bin/bash
# select random SNP from a VCF file.
[ $# -lt 2 ] && echo $0 [vcf] [number] && exit 1
vcf=$1
num=$2

if [ ${vcf:0-3} == '.gz' ];then
    # gzvcf.
    GREP='zgrep'
else
    GREP='grep'
fi

($GREP '^#' $vcf; $GREP -v '^#' $vcf | shuf -n $num | sort -k1,1V -k2,2n)
