# PipelineScomatic

Nextflow pipeline for performing somatic variant calling on single cell RNA-seq data using SComatic

This pipeline relies on the creation of a SComatic conda environment to handle all dependencies, as recommended by the SComatic github page (https://github.com/cortes-ciriano-lab/SComatic), and the cloning of the SComatic github repository onto your local working environment. Once set up, the paths to the SComatic conda environment and the SComatic repository can be set in the config file.

SComatic requires cell type labels, and so is a tool designed to be used downstream of clustering and cell type annotation. Initial inputs to the pipeline are a 2 column comma-separated file containing the paths to the sample .bam files in the first column and the corresponding sample names in the second, and a tab-separated file containing 3 columns and a row for every cell you want to analyse across all samples, where the columns are the cell barcode, cell type and sample ID respectively. The sample ID in column 2 of the .csv file MUST match the sample ID in column 3 of the .tsv for the same sample.

## Parameters

All essential parameters can be set in config file. More task specific options may be added in future.

### Key files

  *	in_bams: Path to .csv file containing .bam file paths and sample IDs
  * bc2cell: Path to .tsv file containing the barcode, cell type and sample ID for each cell in the analysis

### Paths

  * scomatic: Path to cloned SComatic github respository
  * SCconda: Path to scomatic conda environment (setup detailed on SComatic github)

### Accessory files

  * ref_genome: Path to reference genome used for mapping
  * rna_editing: Path to file containing RNA editing genome positions, comes with SComatic (uses params.scomatic as initial path so should not need changing unless the file is manually relocated)
  * pons: Path to PoNs file, comes with SComatic, comes with SComatic (uses params.scomatic as initial path so should not need changing unless the file is manually relocated)

### Script options

  * chrom_opt: option for BaseCellCounter.py, which chromosomes to analyse. Default in config file is all.
