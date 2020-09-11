#!/bin/bash
lostruct="/opt/bio/local_pca/templated/run_lostruct.R"
:>run.cmd
for pop in CJ+TH CJ+CH CJ+HZ CJ+LM;
do
    out="$pop.out"
    rm -r $out
    mkdir $out
    cut -f1 $pop >$out/sample.id
    for i in `seq 1 24`;
    do
        bcftools convert -r LG$i -S $out/sample.id -Oz -o $out/LG$i.vcf.gz dj.vcf.gz
        bcftools index $out/LG$i.vcf.gz
    done
    echo "$lostruct -i $out \
    -I popmap.tsv -t snp \
    -s 20 -o $out/res; Rscript -e \"templater::render_template('summarize_run.Rmd',output='$out/res/$pop.html',change.rootdir=TRUE)\"" >>run.cmd    
done
ParaFly -c run.cmd -CPU 4

    