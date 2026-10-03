process GENOMICSDB_IMPORT {


	tag "interval:$interval.baseName"
	publishDir "${params.outdir}/genomicsdb", mode: 'symlink'	

	
	input:
	path sample_map
	path interval

	output:
	tuple val(interval.baseName), path("genomicsdb_${interval.baseName}")
	
	script:
	"""
	
	gatk GenomicsDBImport \\
	--sample-name-map ${sample_map} \\
	--genomicsdb-workspace-path genomicsdb_${interval.baseName} \\
	--merge-input-intervals \\
	--batch-size 50 \\
	-L ${interval} \\
	--reader-threads ${task.cpus}

	"""
	


}
