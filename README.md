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
* To run the analysis in the manuscript, 
	* Download and install R packages "stats", "ROCR" and "pROC"
	* Download the data files from GEO Accession # GSE62820 (HER2 experimental data ) and #GSE62944 (TCGA reprocessed RNA-Seq data and clinical data).
	* Set the working directory to Analysis_datasets with all the downloaded files.
	*  To make the HER2 predictions on TCGA breast cancer samples, BinReg 2 algorithm in MatLab platform was used. Using our HER2 signature datasets as training samples and TCGA breast cancer datasets as test samples, the predictions were generated. We ran the BinReg 2 for 200 genes with 2 megatons and quantile normalization (-g 200 -m 2 -q) to minimize the batch effects between training and test samples. The original outputs from BinReg2 is located at Analysis_datasets/10_14_predictions_raw. The output predictions from each method are tabulated and are located in the Analysis_datasets folder for further analysis on the data.
	* To classify TCGA  lung adenocarcinoma and squamous carcinoma samples using "Random Forest" classifier with 10-fold cross validation, the R script at Code/Classify_luad_vs_lusc.R was used. The output predictions are located in the Analysis_datasets folder for further analysis on the data.
	* Run the TCGA_20_manuscript_analysis.Rmd file. Our results are stored as TCGA_20_manuscript_analysis.html file.


### Who do I talk to? ###

* Mumtahena Rahman. [moom.rahman@utah.edu](mailto:moom.rahman@utah.edu)
* Stephen R Piccolo. [https://piccolo.byu.edu](https://piccolo.byu.edu)
