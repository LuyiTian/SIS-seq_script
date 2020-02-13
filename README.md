# SIS-seq_script
This repository contains the scripts and data for for manuscript **SIS-seq, a single cell molecular ‘time machine’, identifies Bcor as a regulator of emergency dendritic cell fate**.

## data structure

`script` folder contains all code that been used in this study.

* [SIS-seq_data_analysis.Rmd](https://github.com/LuyiTian/SIS-seq_script/blob/master/scripts/SIS-seq_data_analysis.Rmd) R script to analyse SIS-seq data. Other random script for data exploration and plotting is in [scripts/misc/sis_seq](https://github.com/LuyiTian/SIS-seq_script/tree/master/scripts/misc/sis_seq). This part is related to Figure 1.

* [CRISPR_screen_analysis.R](https://github.com/LuyiTian/SIS-seq_script/blob/master/scripts/CRISPR_screen_analysis.R) R script to analyse CRISPR screening results. The genes used for CRISPR screen is from the SIS-seq data analysis. Other random script for data exploration and plotting is in [scripts/misc/crispr_screen](https://github.com/LuyiTian/SIS-seq_script/tree/master/scripts/misc/crispr_screen). This part is related to Figure 2.

* [BcorKO_scRNAseq_analysis.Rmd](https://github.com/LuyiTian/SIS-seq_script/blob/master/scripts/BcorKO_scRNAseq_analysis.Rmd) R script to analyse scRNAseq from Bcor knockout and wild type cells. This part is related to Figure 3.

`data` folder contains all sequencing data that been used in this study.

* [clone_split](https://github.com/LuyiTian/SIS-seq_script/tree/master/data/clone_split) This folder contains all SIS-seq data, including gene counting matrix and sister-well cell fate measures.


* [crispr_screen](https://github.com/LuyiTian/SIS-seq_script/tree/master/data/crispr_screen) This folder contains all CRISPR screen sequencing data. The R object include processed data using the `CRISPR_screen_analysis.R` script.

* [BcorKO_scRNAseq](https://github.com/LuyiTian/SIS-seq_script/tree/master/data/BcorKO_scRNAseq) This folder contains scRNAseq data of Bcor knockout and wild type cells. 
