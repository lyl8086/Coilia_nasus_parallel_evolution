# Scripts used in parallel evolution analysis of *Coilia nasus*
<u>**Note**: many of these scripts are highly customized, which can only be used in this study.</u>

## Folders:

- SNP_Filtering
  - **VCF_filter_multi.pl**:  do various filtering on VCF

    ```
    VCF_filter_multi.pl -h
    ```

  - **select_random_VCF.sh**: select random SNP from a given VCF file

    ```
    select_random_VCF.sh vcf_file number_to_retain 
    ```

  - **calc_ld_cov.pl**: calculate genome coverage based the given LD extent and SNP [positions](Examples/snp.index)

    ```
    calc_ld_cov.pl snp_index_file samtools_fai_file ld_extent_in_bp
    ```
    
  - **vcf_converter.pl**: convert VCF to [BayeScan](http://cmpg.unibe.ch/software/BayeScan/)/[Genepop](https://kimura.univ-montp2.fr/~rousset/Genepop.htm)/[Structure](http://web.stanford.edu/group/pritchardlab/structure.html)/[Arlequin](http://cmpg.unibe.ch/software/arlequin35/)/[BayEnv](https://bitbucket.org/tguenther/bayenv2_public/src/master/)/[BayesAss](https://github.com/brannala/BA3)

    format: 

    > bayescan  -> [BayeScan](http://cmpg.unibe.ch/software/BayeScan/)
    >
    > genepop   -> [Genepop](https://kimura.univ-montp2.fr/~rousset/Genepop.htm)
    >
    > str              -> [Structure](http://web.stanford.edu/group/pritchardlab/structure.html)
    >
    > arp             -> [Arlequin](http://cmpg.unibe.ch/software/arlequin35/)
    >
    > bayenv      -> [BayEnv](https://bitbucket.org/tguenther/bayenv2_public/src/master/)
    >
    > bayesass   -> [BayesAss](https://github.com/brannala/BA3)
    >
    > all               -> all the above

    ```
    vcf_converter.pl vcf popmap format prefix
    ```

- Allele_Frequency
  
  - **calc_FWA.pl**: calculate the frequency of Fresh Water favored Allele (FWA) given a VCF, a [popmap](Examples/popmap) and the name of reference population
  
    ```
    calc_FWA.pl vcf popmap reference_pop_name
    ```
  
  - **calc_maf_ref.pl**: calculate the frequency of major allele in the reference population, and the corresponding frequency on the rest populations
  
    ```
    calc_maf_ref.pl vcf popmap reference_pop_name
    ```
  
- PCA
  - **multi_VCF_PCA_XXX.R**: R script to do PCA plot on multiple VCF files in a folder based on [adegenet](https://github.com/thibautjombart/adegenet)/[GCTA](https://cnsgenomics.com/software/gcta/)/[SeqVarTools](https://github.com/smgogarten/SeqVarTools)/[SNPRelate](https://github.com/zhengxwen/SNPRelate)

    ```
    Rscript multi_VCF_PCA_XXX.R input_foler_containing_vcf_files
    ```

- Outliers

  - **do_pcadapt.sh**: call **pcadapt.r** to do [pcadapt](https://github.com/bcm-uga/pcadapt) on each population

    ```
    bash do_pcadapt.sh
    ```

  - **Fisher_exact_test.pl**: do Fisher's Exact test given a VCF file, a popmap and the reference population name to compare with

    ```
    Fisher_exact_test.pl vcf popmap reference_pop_name
    ```

- INV

  - **do_LDna.sh**: call **LDna.R** to do [LDna](https://github.com/petrikemppainen/LDna) analysis on each chromosome using the LD results of [VCFtools](https://github.com/vcftools/vcftools)

    ```
    bash do_LDna.sh
    ```

  - **LDna.R**: R script to do LDna analysis based on the r2 calculated by VCFtools

    ```
    Rscript LDna.R chr_xx.ld.gz out_prefix vcf.gz file_of_individuals_id number_of_threads
    ```

  - **do_localPCA.sh**: call [lostruct](https://github.com/petrelharp/local_pca) to do local PCA on each chromosome

    ```
    bash do_localPCA.sh
    ```

  - **plot_local_PCA.R**: R script to plot the results of local PCA 

    ```
    Rscript plot_local_PCA.R
    ```

  - **pca_and_boxplot_XXX.R**: R script to plot PCA and heterozygosity for multiple VCF in a folder

    ```
    Rscript pca_and_boxplot_XXX.R input_foler_containing_vcf_files
    ```

  - **plink_LDna.R**: R script to do LDna analysis based on the LD results of plink

    - Parameters need to be adjusted

    ```
    Rscript plink_LDna.R plin.ld.gz out_prefix
    ```

- Gene_Annotation

  - **topGO.R**: R script to do GO enrichment analysis using [topGO](https://www.bioconductor.org/packages/release/bioc/html/topGO.html) given a [**genes_to_go**](Examples/genes_to_go) and a **[genes_of_interest](Examples/genes_of_interest)** file 

    ```
    Rscript topGO.R
    ```

  - **GOA_genes_to_GO.pl**: extract the genes id to GO mapping file (genes_to_go) given the results of blastp (tsv) and [GOA](ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/) data base file

    ```
    GOA_genes_to_GO.pl blast_tsv_file GOA_db_file
    ```



------

**Questions should be addressed to [Yulong Li](mailto:liyulong12@mails.ucas.ac.cn)**

