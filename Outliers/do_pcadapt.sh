#!/bin/bash

for d in TH/ CH/ HZ/ LM/; 
do 
    cd $d; Rscript ../pcadapt.r; 
    cd ../; 
done

cat TH/pcadapt.outliers.id CH/pcadapt.outliers.id \
HZ/pcadapt.outliers.id LM/pcadapt.outliers.id \
| sort -V | uniq -c >outliers_pca.id

awk '$1==4 {print $2}' outliers_pca.id >all_4pop.outliers


