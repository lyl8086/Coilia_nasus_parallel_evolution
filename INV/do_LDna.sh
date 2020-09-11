#!/bin/bash
vcf="dj.vcf.gz"
for pop in popmap;
do
    out="Three_out"
    [ -f "$out/LDna.done" ] && continue # check
    rm -r $out
    :>$pop.r2.cmd
    :>$pop.LDna.cmd
    mkdir $out
    cut -f1 $pop >$out/sample.id
    for i in `seq 1 24`;
    do
        echo "vcftools --gzvcf $vcf --keep $out/sample.id --chr LG$i --geno-r2 --out $out/LG$i" >>$pop.r2.cmd
        echo "Rscript LDna.R $out/LG$i.geno.ld $out/LG$i.LDna $vcf $out/sample.id" >>$pop.LDna.cmd
    done
    ParaFly -c $pop.r2.cmd -CPU 24
    ParaFly -c $pop.LDna.cmd -CPU 24
    # get PCA
    for i in $out/*.recode.vcf;
    do
        b=${i%".recode.vcf"}
        vcf_pca.sh $i $b.PCA
    done
    touch $out/LDna.done
done

