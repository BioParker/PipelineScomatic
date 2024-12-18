# PipelineScomatic

Nextflow pipeline for performing somatic variant calling on single cell RNA-seq data using SComatic

This pipeline relies on the creation of a SComatic conda environment to handle all dependencies, as recommended by the SComatic github page (https://github.com/cortes-ciriano-lab/SComatic), and the cloning of the SComatic github repository onto your local working environment. Once set up, the paths to the SComatic conda environment and the SComatic repository can be set in the config file.

Initial inputs to the pipeline are a 2 column comma-separated file containing the paths to the sample .bam files in the first column and the corresponding sample names in the second, and a tab-separated file containing 3 columns and a row for every cell you want to analyse across all samples, where the columns are the cell barcode, cell type and sample ID respectively. The sample ID in column 2 of the .csv file MUST match the sample ID in column 3 of the .tsv for the same sample. 
