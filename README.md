# README #

This repository includes code for processing RNA-Seq FASTQ files and clinical data from The Cancer Genome Atlas. In addition, we have included the codes for data analysis in the manuscript "RNA-Sequencing data for 7706 tumor samples across 20 cancer types from The Cancer Genome Atlas".  

### What is this repository for? ###

* We have used the 'Rsubread' R package to align and summarize reads at the gene level for 7706 TCGA RNA-Seq tumor samples. The R scripts can also be used to process new samples. We have also included the codes for compiling clinical data available for these tumors into a matrix format and matched the IDs for ease of matching phenotypes with clinical variables. 
* We have providied the codes and small datasets necessary for analyzing the analysis scenarios as described in the manuscript.

### How do I get set up? ###

* Clone the git repository to you local computing area using the url ttps://github.com/mumtahena/TCGA_RNASeq_clinical.git
* Obtain access to download the raw TCGA data or download any other raw FASTQ files that you want to process
* Place the raw files in fastq or fastq.gz formats in the FASTQ folder of the repository
* Install R packages "Rsubread", "limma", "edgeR" and "tools" for processing new samples. Install "stats" R package and its dependencies if you just want to run the analysis included in the manuscript.
* Process the data using our code for downstream analysis.

### Who do I talk to? ###

* Mumtahena Rahman. [moom.rahman@utah.edu](mailto:moom.rahman@utah.edu)
* Stephen R Piccolo. [https://piccolo.byu.edu](https://piccolo.byu.edu)
