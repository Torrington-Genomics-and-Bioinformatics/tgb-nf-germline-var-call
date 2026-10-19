process BWA_MEM {
	
	tag "$sample_id"
	publishDir "${params.outdir}/aligned", mode: 'symlink'

	input:
	tuple val(sample_id), path(fastq1), path(fastq2)
	path ref
	path ref_idx
	output:
	tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

	script:
	"""
	
	# Prepare RG info
	LB="${sample_id}_lib1"
	PL="ILLUMINA"
	FLOWCELL_LANE=\$(zcat ${fastq1} | head -1 | cut -d: -f3,4 --output-delimiter=".")

	# Run alignment while adding read groups 
	bwa mem \\
	-K 10000000 -Y \\
	-R "@RG\\tID:\${FLOWCELL_LANE}\\tSM:${sample_id}\\tLB:\${LB}\\tPL:\${PL}\\tPU:\${FLOWCELL_LANE}" \\
	${ref} \\
	${fastq1} ${fastq2} \\
	| samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
	
	# Index sorted bam
	samtools index ${sample_id}.sorted.bam



	"""




}
