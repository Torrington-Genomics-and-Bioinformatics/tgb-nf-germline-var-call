process GENOMICSDB_IMPORT {


	tag "interval:${name}"	

	
	input:
	path sample_map
	tuple val(name), val(gendb)

	output:
	tuple val(name), val(gendb)
	
	script:
	"""
	
	DB="${gendb}"

	gatk GenomicsDBImport \\
	--sample-name-map ${sample_map} \\
	--genomicsdb-update-workspace-path "\$DB" \\
	--batch-size 50 \\
	--reader-threads ${task.cpus}

	"""
	


}
