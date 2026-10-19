process QC {
	
	tag "$sample_id"
	publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.fastp.*'
	
	input:
	tuple val(sample_id), path(fastq1), path(fastq2)
	
	output:
	tuple val(sample_id), path("${sample_id}_R1.qcd.fq.gz"), path("${sample_id}_R2.qcd.fq.gz"), emit: reads
	path("${sample_id}.fastp.{json,html}"), emit: report
	

	script:
	"""

	fastp \\
	-i ${sample_id}_R1.fastq.gz \\
	-I ${sample_id}_R2.fastq.gz \\
	-o ${sample_id}_R1.qcd.fq.gz \\
	-O ${sample_id}_R2.qcd.fq.gz \\
	--detect_adapter_for_pe \\
	-j ${sample_id}.fastp.json \\
	-h ${sample_id}.fastp.html
	
	"""

}
