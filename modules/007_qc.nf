// v3.1-qc: alignment QC on the final (recalibrated) BAM, and one QC table + report per run.
// These processes only read pipeline outputs; they do not change the BAM, the gVCF or any call.

process SAMPLE_QC {

	tag "$sample_id"
	publishDir "${params.outdir}/qc", mode: 'copy'

	input:
	tuple val(sample_id), path(bam), path(bai), path(fmt)
	path ref
	path ref_index
	path ref_dict
	path intervals

	output:
	path("${sample_id}.{hs_metrics.txt,per_target_coverage.txt.gz,samtools.stats}"), emit: files

	script:
	"""
	# Coverage, on-target, fold-80 and dropout on the targets used for calling (bait = target)
	gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" CollectHsMetrics \\
	-I ${bam} -R ${ref} \\
	-BI ${intervals} -TI ${intervals} \\
	--PER_TARGET_COVERAGE ${sample_id}.per_target_coverage.txt \\
	-O ${sample_id}.hs_metrics.txt
	gzip ${sample_id}.per_target_coverage.txt

	samtools stats -@ ${task.cpus} --reference ${ref} ${bam} > ${sample_id}.samtools.stats
	"""
}

process BATCH_QC {

	publishDir "${params.outdir}/qc", mode: 'copy'

	input:
	path(files)

	output:
	path("batch_qc.tsv")
	path("batch_qc_report.html")
	path("multiqc_report.html")

	script:
	def title = file(params.outdir).name
	"""
	batch_qc.py . batch_qc.tsv batch_qc_report.html "${title}" "${params.pipeline_version}"
	multiqc -q -f --no-data-dir -n multiqc_report.html --title "${title}" .
	"""
}
