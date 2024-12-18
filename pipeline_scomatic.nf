#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */


//Output directories

params.cpdir="bams.dir"
params.bc2celldir="bc2cell.dir"
params.cellbamsdir="cellbams.dir"
params.basecountdir="basecount.dir"
params.mergecountsdir="mergecounts.dir"
params.basecall1dir="basecall1.dir"
params.basecall2dir="basecall2.dir"

/*
 * Processes
 */


//copy bams to nextflow directory and index
//possibly an unnecessary step but ensures that bam indexes are up to date 
process cp_idx {

        conda params.SCconda

	input:
        tuple path(input_bam), val(sample_id)
	    val cpdir

	output:
	    tuple path("${cpdir}/${sample_id}.bam"), path("${cpdir}/${sample_id}.bam.bai"), val("${sample_id}")

	script:
	"""
	mkdir '$cpdir';
        cp '$input_bam' '$cpdir'/'$sample_id'.bam;
	samtools index '$cpdir'/'$sample_id'.bam
	"""
}

//split barcode to cell type assignment table into sample specific tables
process mk_bc2cell {

	publishDir   params.bc2celldir, mode: 'symlink'

	input:
	   tuple path(indexed_bam), path(bam_index), val(sample_id)
	   path bc2celltbl

	output:
	   tuple path("${indexed_bam}"), path("${bam_index}"), val("${sample_id}"), path("${sample_id}_bc2cell.tsv")

	script:
	"""
	grep -e "Index" -e '$sample_id' '$bc2celltbl' | awk '{print \$1"\\t"\$2}' > '$sample_id'_bc2cell.tsv
	"""

}

//divide each bam file into cell type bams
process mk_cellbams {

    publishDir params.cellbamsdir, mode: 'symlink'
    
    conda params.SCconda

    input:
        tuple path(indexed_bam), path(bam_index), val(sample_id), path(bc2cell)
        path scomatic_path

    output:
	tuple path("${sample_id}.*.bam"), val("${sample_id}"), emit: out_bam
	path "${sample_id}.*.bam.bai", emit: out_index

    script:
    """
    python '$scomatic_path'/scripts/SplitBam/SplitBamCellTypes.py \
            --bam '$indexed_bam' \
            --meta '$bc2cell' \
            --id '$sample_id' \
            --max_nM 1 \
            --min_MQ 255 \
    """
}

//base count for each cell type bam. optional: true in case BaseCellCounter.py cannot output results for particular cell type
process base_count {

    publishDir "${params.basecountdir}/${sample_id}", mode: 'symlink'

    conda params.SCconda

    input:
        tuple path(cell_bam), val(sample_id), path(cell_bam_idx)
        path scomatic_path
        path ref_genome
        val chrom_opt

    output:
        tuple path("${sample_id}.*.tsv"), val("${sample_id}"), optional: true

    script:
    """
    python '$scomatic_path'/scripts/BaseCellCounter/BaseCellCounter.py \
	    --bam '$cell_bam' \
            --ref '$ref_genome' \
            --chrom '$chrom_opt'
    """
}

//merging base count matrices
process merge_counts {

    publishDir params.mergecountsdir, mode: 'copy'

    conda params.SCconda

    input:
        tuple path(counts_folder), val(sample_id)
	path scomatic_path

    output:
        tuple path("${sample_id}_aggrCounts.tsv"), val("${sample_id}")

    script:
    """
    python '$scomatic_path'/scripts/MergeCounts/MergeBaseCellCounts.py \
            --tsv_folder '$counts_folder' \
            --outfile '$sample_id'_aggrCounts.tsv \
    """
}

//first base calling step
process base_calling_1 {

    publishDir params.basecall1dir, mode: 'copy'

    conda params.SCconda

    input:
        tuple path(aggr_counts), val(sample_id)
        path scomatic_path
        path ref_genome

    output:
        tuple path("${sample_id}.calling.step1.tsv"), val("${sample_id}"), optional: true

    script:
    """
    python '$scomatic_path'/scripts/BaseCellCalling/BaseCellCalling.step1.py \
            --infile '$aggr_counts' \
            --outfile '$sample_id' \
	    --ref '$ref_genome'
    """
}

//second base calling step
process base_calling_2 {

    publishDir params.basecall2dir, mode: 'copy'

    conda params.SCconda

    input:
        tuple path(calls1), val(sample_id)
        path scomatic_path
        path rna_edits
	path pons

    output:
        tuple path("${sample_id}.calling.step2.tsv"), val("${sample_id}"), optional: true

    script:
    """
    python '$scomatic_path'/scripts/BaseCellCalling/BaseCellCalling.step2.py \
            --infile '$calls1' \
            --outfile '$sample_id' \
            --editing '$rna_edits' \
	    --pon '$pons'
    """
}

/*
 * Workflow
 */


workflow {

    //create input channel
    bampaths = Channel.fromPath(params.in_bams)
	       .splitCsv( header: true )
	       .map { row -> tuple(file(row.bam), row.sample) }
    
    //outdirs
    cpdir = params.cpdir
    def bcountdir = "${projectDir}/${params.basecountdir}"

    //main files
    bc2celltbl = file(params.bc2cell)

    //script path
    scomatic_path = params.scomatic

    //accessory files    
    ref_genome = file(params.ref_genome)
    edits = file(params.rna_editing)
    pons = file(params.pons)

    //script options
    chrom_opt = params.chrom_opt


    bampaths.view()
    
    cp_idx(bampaths,
	   cpdir)

    mk_bc2cell(cp_idx.out,
	       bc2celltbl)

    mk_cellbams(mk_bc2cell.out,
                scomatic_path)

    // transform cellbams output to get tuples of file-sample-index triplets
    def bamsam = mk_cellbams.out.out_bam.map { it -> it[0] + it[1] }
                            .flatMap{ it -> def sampl = it[-1]
                                            it[0..-2].collect { element -> [element, sampl] }[0..-2] }
			    .map { it -> it + "${it[0]}.bai" }			    	    	    
			    //.view()


    base_count(bamsam,
               scomatic_path,
               ref_genome,
               chrom_opt)

    //construct path to publishDir tsv folder for each sample
    def sample_tsv_path = base_count.out.map { it -> it[1] }
				    .flatten()
				    .unique()
				    .map { it -> ["${bcountdir}/${it}", it] }
				    .view()

    merge_counts(sample_tsv_path,
		 scomatic_path)

    base_calling_1(merge_counts.out,
		   scomatic_path,
		   ref_genome)

    base_calling_2(base_calling_1.out,
		   scomatic_path,
		   edits,
		   pons)
}

