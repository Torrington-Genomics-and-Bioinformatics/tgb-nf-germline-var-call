process QC {
	
	tag "$sample_id"
	publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*fastp.*'
	
	input:
	tuple val(sample_id), path(fastq1), path(fastq2), path(fastqs)
	
	output:
	tuple val(sample_id), path("${sample_id}_R1.qcd.fq.gz"), path("${sample_id}_R2.qcd.fq.gz"), path("${sample_id}_S.qcd.fq.gz"), emit: reads
	path("${sample_id}.fastp.{json,html}"), emit: report
	path("${sample_id}.singletons_fastp.{json,html}"), optional: true, emit: se_report
	

	script:
	// v3.2: unpaired reads (fastq_s) are trimmed single-end with the same fastp defaults; without them
	// an empty placeholder is passed on and the paired command is unchanged
	def se_cmd = fastqs.name == 'NO_SINGLETONS' ? ": > ${sample_id}_S.qcd.fq.gz" :
		"fastp -i ${fastqs} -o ${sample_id}_S.qcd.fq.gz -j ${sample_id}.singletons_fastp.json -h ${sample_id}.singletons_fastp.html"
	"""

	fastp \\
	-i ${sample_id}_R1.fastq.gz \\
	-I ${sample_id}_R2.fastq.gz \\
	-o ${sample_id}_R1.qcd.fq.gz \\
	-O ${sample_id}_R2.qcd.fq.gz \\
	--detect_adapter_for_pe \\
	-j ${sample_id}.fastp.json \\
	-h ${sample_id}.fastp.html

	${se_cmd}
	
	"""

}
